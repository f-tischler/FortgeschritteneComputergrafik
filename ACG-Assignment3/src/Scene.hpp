#ifndef Scene_H
#define Scene_H

#include "Ray.hpp"
#include "Geometry.hpp"
#include <memory>
#include <vector>

struct IntersectionInfo
{
	Ray ray;
	Vector normal;
	Vector hitpoint;
	float distance;

	Geometry* geometry;
};

class Scene
{
public:
	using GeometryList = std::vector<std::unique_ptr<Geometry>>;

	Scene() { }
	
	bool Intersect(const Ray& ray, IntersectionInfo& info) const
	{
		info.ray = ray;
		info.normal = Vector(0, 0, 0);
		info.hitpoint = Vector(0, 0, 0);
		info.distance = std::numeric_limits<float>::max();
		info.geometry = nullptr;

		auto hit = false;
		for(const auto& geometry : _geometry)
		{
			float distance;
			Vector normal;

			if(geometry->Intersects(ray, normal, distance))
			{
				if(distance < info.distance)
				{
					info.normal = normal;
					info.distance = distance;
					info.geometry = geometry.get();
					info.hitpoint = ray.GetOrigin() + ray.GetDirection() * distance;

					hit = true;
				}
			}
		}

		return hit;
	};

	void AddGeometry(GeometryList::value_type&& geometry)
	{
		_geometry.push_back(std::move(geometry));
	};

	const auto& GetGeometry() const { return _geometry; }

private:
	
	GeometryList _geometry;
};

#endif // Scene_H
