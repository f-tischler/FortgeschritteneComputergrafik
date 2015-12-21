#ifndef TriangleMesh_h_included
#define TriangleMesh_h_included

#include <vector>
#include <algorithm>

#include "Consts.hpp"
#include "Geometry.hpp"
#include "BoundingBox.hpp"
#include "Matrix.hpp"

typedef Vector Vertex;

class TriangleMesh : public Geometry
{
public:
	TriangleMesh(Vector emission_,Vector color_, Refl_t refl_) :
		Geometry(emission_, color_, refl_), mDivCount(15)
	{
		mChanged = false;
	}

	void addVertex(const Vertex& vertex)
	{
		mVertices.push_back(vertex);
		mChanged = true;
	}

	virtual double Intersect(const Ray &ray, Vector& normal, bool culling)
	{
		assert(mVertices.size() % 3 == 0);
		double distance = DBL_MAX;
		bool found = false;

		if (!m_BB.intersect(ray))
			return 0.0;

		for (size_t i = 0; i < mBoundingBoxes.size(); i++)
		{
			if (mBoundingBoxes[i].intersect(ray))
			{
				auto start = i*mDivCount;
				auto end = std::min(mVertices.size(), (size_t)i*mDivCount + mDivCount);

				for (size_t j = start; j < end; j += 3)
				{
					auto d = 0.0;
					auto V1 = mVertices[j];
					auto V2 = mVertices[j + 1];
					auto V3 = mVertices[j + 2];

					auto edgeA = V2 - V1;
					auto edgeB = V3 - V1;
					auto edgeACrossEdgeB = edgeA.Cross(edgeB);

					auto currentNormal = edgeACrossEdgeB.Normalized();

					if (intersection(ray, d, V1, V2, V3, currentNormal, culling))
					{
						distance = std::min(distance, d);
						normal = currentNormal;
						found = true;
					}
				}
			}
		}

		if (!found)
			return 0.0f;
		
		return distance;
	}
    
    double toRad(double angleDeg) {
        return M_PI/180 * angleDeg;
    }

	void rotateX(double angleDeg)
	{
        double angleRad = toRad(angleDeg);
        
        double tmp[16] =
        {
            1.0,           0.0,            0.0, 0.0,
            0.0, cos(angleRad), -sin(angleRad), 0.0,
            0.0, sin(angleRad),  cos(angleRad), 0.0,
            0.0,           0.0,            0.0, 1.0
        };
        
		transformPoints(Matrix(tmp));
		mChanged = true;
	}

	void rotateY(float angleDeg)
	{
        double angleRad = toRad(angleDeg);
        
        double tmp[16] =
        {
            cos(angleRad) , 0.0, sin(angleRad), 0.0,
            0.0           , 1.0,           0.0, 0.0,
            -sin(angleRad), 0.0, cos(angleRad), 0.0,
            0.0           , 0.0,           0.0, 1.0
        };
        
		transformPoints(Matrix(tmp));
	}

	void rotateZ(float angleDeg)
	{
        double angleRad = toRad(angleDeg);
        
        double tmp[16] =
        {
            cos(angleRad), -sin(angleRad), 0.0, 0.0,
            sin(angleRad), cos(angleRad) , 0.0, 0.0,
            0.0          , 0.0           , 1.0, 0.0,
            0.0          , 0.0           , 0.0, 1.0
        };
        
		transformPoints(Matrix(tmp));
	}

	void translate(double x, double y, double z)
	{
        double tmp[16] =
        {
            1.0, 0.0, 0.0, x  ,
            0.0, 1.0, 0.0, y  ,
            0.0, 0.0, 1.0, z  ,
            0.0, 0.0, 0.0, 1.0
        };
		transformPoints(Matrix(tmp));
	}

	void scale(double x, double y, double z)
	{
        double tmp[16] =
        {
            x  , 0.0, 0.0, 0.0,
            0.0, y  , 0.0, 0.0,
            0.0, 0.0,   z, 0.0,
            0.0, 0.0, 0.0, 1.0
        };
        
		transformPoints(Matrix(tmp));
	}

	void scaleToWidth(double width)
	{
		if (mChanged)
		{
			buildCollision();
			mChanged = false;
		}

		auto newScale = width / m_BB.width();
		scale(newScale, newScale, newScale);
	}

	virtual size_t getType() const { return eGeometryType_TriangleMesh; }


	virtual SPGeometry clone()
	{
		auto mesh = std::make_shared<TriangleMesh>(emission, color, refl);
		mesh->mVertices = mVertices;
		mesh->mChanged = true;
		mesh->buildCollision();

		return mesh;
	}

	void buildCollision()
	{
		if (!mChanged)
			return;

		mChanged = false;

		mBoundingBoxes.clear();
		m_BB.reset();

		for (size_t div = 0; div < mVertices.size() / mDivCount; div++)
		{
			BoundingBox<Vector> bb;

			for (auto i = div*mDivCount; i < mVertices.size(); i++)
			{
				if (i > div*mDivCount + mDivCount)
					break;

				bb.addPoint(mVertices[i]);
				m_BB.addPoint(mVertices[i]);
			}

			mBoundingBoxes.push_back(bb);
		}



		if (mVertices.size() % mDivCount != 0)
		{
			BoundingBox<Vector> bb;

			auto start = mVertices.size() - mVertices.size() % mDivCount;
			auto end = mVertices.size();

			for (auto i = start; i < end; i++)
			{
				bb.addPoint(mVertices[i]);
				m_BB.addPoint(mVertices[i]);
			}

			mBoundingBoxes.push_back(bb);
		}
	}

private:
	std::vector<BoundingBox<Vector>> mBoundingBoxes;
	std::vector<Vertex> mVertices;
	bool mChanged;
	int mDivCount;

	void transformPoints(Matrix m)
	{
		for (auto it = mVertices.begin(); it != mVertices.end(); ++it)
		{
			auto& v = *it;
			v = m * v;
		}

		mChanged = true;
	}

	bool intersection(const Ray& ray, double& distance, const Vector& V1, const Vector& V2, const Vector& V3, const Vector& normal, bool culling) const
	{
		Vector e1, e2;  //Edge1, Edge2
		Vector P, Q, T;
		double det, inv_det, u, v;
		double t;

		auto D = ray.dir;
		auto O = ray.org;

#define SUB(a,b,c) a = (b) - (c)
#define CROSS(a,b,c) a = (b).Cross(c)
#define DOT(a,b) (a).Dot(b)
#define EPSILON 0.000001f


		// implementing the single/double sided feature
		if (culling && DOT(ray.dir, normal) > 0)
			return false; // back-facing surface 

		//Find vectors for two edges sharing V1
		e1 = V2 - V1;
		e2 = V3 - V1;

		//Begin calculating determinant - also used to calculate u parameter
		P = D.Cross(e2);

		//if determinant is near zero, ray lies in plane of triangle
		det = e1.Dot(P);

		//NOT CULLING
		if (det > -EPSILON && det < EPSILON)
			return false;


		inv_det = 1.f / det;
		//calculate distance from V1 to ray origin
		T = O - V1;

		//Calculate u parameter and test bound
		u = T.Dot(P) * inv_det;

		//The intersection lies outside of the triangle
		if (u < 0.f || u > 1.f)
			return false;

		//Prepare to test v parameter
		Q = T.Cross(e1);

		//Calculate V parameter and test bound
		v = D.Dot(Q) * inv_det;

		//The intersection lies outside of the triangle
		if (v < 0.f || u + v  > 1.f)
			return false;

		t = e2.Dot(Q) * inv_det;

		//ray intersection
		if (t > EPSILON)
		{
			distance = t;
			return true;
		}

		return false;
	}
};

#endif