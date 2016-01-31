#ifndef TriangleMesh_h_included
#define TriangleMesh_h_included

#include <vector>
#include <algorithm>

#include "glm/gtx/transform.hpp"
#include "Material.hpp"

#include "Geometry.hpp"
#include "BoundingBox.hpp"
#include "Common.hpp"

typedef Vector Vertex;

class TriangleMesh : public Geometry
{
public:
	explicit TriangleMesh(const Material& mat) :
		Geometry(mat, eGeometryType::eGeometryType_TriangleMesh), mDivCount(15)
	{
		mChanged = false;
	}

	void addVertex(const Vertex& vertex)
	{
		mVertices.push_back(vertex);
		mChanged = true;

		GetBoundingBox().addPoint(vertex);
	}

	virtual bool IntersectsGeometry(const Ray& ray, Vector& normal, float& distance)
	{
		assert(mVertices.size() % 3 == 0);
		distance = 9999999.f;
		bool found = false;

		for (size_t j = 0; j < mVertices.size(); j += 3)
		{
			auto d = 0.0f;
			auto V1 = mVertices[j];
			auto V2 = mVertices[j + 1];
			auto V3 = mVertices[j + 2];

			auto edgeA = V2 - V1;
			auto edgeB = V3 - V1;
			auto edgeACrossEdgeB = glm::cross(edgeA,edgeB);

			auto currentNormal = glm::normalize(edgeACrossEdgeB);

			if (intersection(ray, d, V1, V2, V3, currentNormal, true))
			{
				distance = std::min(distance, d);
				normal = currentNormal;
				found = true;
			}
		}
		
		return found;
	}

	void scale(const Vector& v)
	{
		transform(glm::scale(glm::mat4x4(1.0f),v));
	}

	void translate(const Vector& v)
	{
		GetBoundingBox().reset();

		for (auto& vert : mVertices)
		{
			vert += v;

			GetBoundingBox().addPoint(vert);
		}
	}

	void rotate(const Vector& v)
	{
		transform(glm::rotate(glm::mat4x4(1.0f), v.x*(PI / 180.f), Vector(1.f,0.f,0.f)));
		transform(glm::rotate(glm::mat4x4(1.0f), v.y*(PI / 180.f), Vector(0.f, 1.f, 0.f)));
		transform(glm::rotate(glm::mat4x4(1.0f), v.z*(PI / 180.f), Vector(0.f, 0.f, 1.f)));
	}

	void transform(glm::mat4x4 transform)
	{
		GetBoundingBox().reset();

		for (auto& v : mVertices)
		{
			glm::vec4 tmp(v.x,v.y,v.z,1.0f);
			tmp = tmp * transform;

			v.x = tmp.x;
			v.y = tmp.y;
			v.z = tmp.z;

			GetBoundingBox().addPoint(v);
		}
	}

private:
	std::vector<Vertex> mVertices;
	bool mChanged;
	int mDivCount;



	bool intersection(const Ray& ray, float& distance, const Vector& V1, const Vector& V2, const Vector& V3, const Vector& normal, bool culling) const
	{
		Vector e1, e2;  //Edge1, Edge2
		Vector P, Q, T;
		float det, inv_det, u, v;
		float t;

		auto D = ray.GetDirection();
		auto O = ray.GetOrigin();

#define SUB(a,b,c) a = (b) - (c)
#define CROSS(a,b,c) a = (b).Cross(c)
#define DOT(a,b) (a).Dot(b)
#define EPSILON 0.000001f


		// implementing the single/double sided feature
		if (culling && glm::dot(D, normal) > 0)
			return false; // back-facing surface 

		//Find vectors for two edges sharing V1
		e1 = V2 - V1;
		e2 = V3 - V1;

		//Begin calculating determinant - also used to calculate u parameter
		P = glm::cross(D,e2);

		//if determinant is near zero, ray lies in plane of triangle
		det = glm::dot(e1,P);

		//NOT CULLING
		if (det > -EPSILON && det < EPSILON)
			return false;


		inv_det = 1.f / det;
		//calculate distance from V1 to ray origin
		T = O - V1;

		//Calculate u parameter and test bound
		u = glm::dot(T,P) * inv_det;

		//The intersection lies outside of the triangle
		if (u < 0.f || u > 1.f)
			return false;

		//Prepare to test v parameter
		Q = glm::cross(T,e1);

		//Calculate V parameter and test bound
		v = glm::dot(D,Q) * inv_det;

		//The intersection lies outside of the triangle
		if (v < 0.f || u + v  > 1.f)
			return false;

		t = glm::dot(e2,Q) * inv_det;

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