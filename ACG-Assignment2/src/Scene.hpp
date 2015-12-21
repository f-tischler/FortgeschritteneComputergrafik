#ifndef Scene_H_included 
#define Scene_H_included

//#define USE_OMP

#include <memory>
#include <vector>

#ifdef USE_OMP
#include <omp.h>
#endif

#include "Sphere.hpp"
#include "TriangleMesh.hpp"
#include "ModelLoader.hpp"

struct IntersectionInfo
{
	bool hit;
	Vector normal;
	Vector hitpoint;
	SPGeometry geometry;
	size_t geometryId;
	double distance;
};

class Scene
{
	BoundingBox<Vector> mBB;
	std::vector<SPGeometry> mGeometries;

public:
	Scene()
	{

	}

	void addGeometry(SPGeometry geometry)
	{
		mBB.addPoint(geometry->getBB().min());
		mBB.addPoint(geometry->getBB().max());

		mGeometries.push_back(geometry);
	}

	bool Intersect(const Ray &ray, IntersectionInfo& info, bool culling = true)
	{
		if (!mBB.intersect(ray))
			return false;

		double T = DBL_MAX;

		info.hit = false;

		#ifdef USE_OMP
			#pragma omp parallel for
		#endif
		for (long i = 0u; i < (long)mGeometries.size(); i++)
		{
			Vector newNormal;
			double distance = mGeometries[i]->Intersect(ray, newNormal, culling);

			#ifdef USE_OMP
				#pragma omp critical
			#endif
			{
				if (distance > 0.0 && distance < T)
				{
					T = distance;
					info.hit = true;
					info.distance = distance;
					info.normal = newNormal;
					info.geometry = mGeometries[i];
					info.geometryId = i;
				}
			}
		}


		info.hitpoint = ray.org + ray.dir * info.distance;
		return info.hit;
	}

	const std::vector<SPGeometry>& geometries() { return mGeometries; }
};

typedef std::shared_ptr<Scene> SPScene;


#endif