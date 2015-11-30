#ifndef Triangle_H_included
#define Triangle_H_included

#include "Vector.hpp"
#include "Ray.hpp"

#include <cassert>
#include <vector>

class Triangle
{
public:
	Triangle() : _a(0.0) { }

	Triangle(const Vector p0_, const Vector p1_, const Vector p2_, const Color &emission_, const Color &color_)
		: _p0(p0_), _p1(p1_), _p2(p2_), _emission(emission_), _color(color_)
	{
		auto edgeA = _p1 - _p0;
		auto edgeB = _p2 - _p0;
        auto edgeACrossEdgeB = edgeA.Cross(edgeB);

		_a = edgeACrossEdgeB.Length() / 2.0;
		_normal = edgeACrossEdgeB.Normalized();
	}

	double area() const { return _a; }
	const Vector& normal() const { return _normal; }

	void addNeighbour(Triangle* neighbour) { _neighbors.push_back(neighbour); }

	bool intersectsSub(const Ray &ray, double& distance, unsigned int& subIndex) const
	{
		for (unsigned int i = 0; i < _subTriangles.size(); i++)
		{
			double tmp;
			if (_subTriangles[i].intersects(ray, tmp))
			{
				distance = tmp;
				subIndex = i;
				return true;
			}
		}

		return false;
	}

	bool intersects(const Ray &ray, double& distance) const
	{
		Vector e1, e2;  //Edge1, Edge2
		Vector P, Q, T;
		double det, inv_det, u, v;
		double t;

		auto V1 = _p0;
		auto V2 = _p1;
		auto V3 = _p2;
		auto D = ray.dir;
		auto O = ray.org;

		#define SUB(a,b,c) a = (b) - (c)
		#define CROSS(a,b,c) a = (b).Cross(c)
		#define DOT(a,b) (a).Dot(b)
		#define EPSILON 0.000001f


		// implementing the single/double sided feature
		if (DOT(ray.dir, _normal) > 0)
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

	size_t getSubTriangleCount() const { return _subTriangles.size(); }

	void split(Triangle(&out)[4]) const
	{
		auto midP1P0 = _p0 + (_p1 - _p0) / 2;
		auto midP2P0 = _p0 + (_p2 - _p0) / 2;
		auto midP2P1 = _p1 + (_p2 - _p1) / 2;

		out[0] = Triangle(midP1P0, midP2P1, midP2P0, _emission, Color(0, 0 ,0));
		out[1] = Triangle(_p0, midP1P0, midP2P0, _emission, Color(0, 0, 0));
		out[2] = Triangle(midP1P0, _p1, midP2P1, _emission, Color(0, 0, 0));
		out[3] = Triangle(midP2P1, _p2, midP2P0, _emission, Color(0, 0, 0));

		const auto epsilon = 0.00001f;
		const auto epsilonVec = Vector(epsilon, epsilon, epsilon);

		assert(Vector::AreEqual(normal(),
								out[0].normal(),
								epsilonVec) && "normals of subtriangles must match the parent's");

		assert(Vector::AreEqual(normal(),
								out[1].normal(),
								epsilonVec) && "normals of subtriangles must match the parent's");

		assert(Vector::AreEqual(normal(),
								out[2].normal(),
								epsilonVec) && "normals of subtriangles must match the parent's");

		assert(Vector::AreEqual(normal(),
								out[3].normal(),
								epsilonVec) && "normals of subtriangles must match the parent's");

// 
// 		//For testing interpolation
// 		out[0] = Triangle(midP1P0, midP2P1, midP2P0, _emission, _color*0.2);
// 		out[1] = Triangle(_p0, midP1P0, midP2P0, _emission, _color*0.4);
// 		out[2] = Triangle(midP1P0, _p1, midP2P1, _emission, _color*0.7);
// 		out[3] = Triangle(midP2P1, _p2, midP2P0, _emission, _color);
	}

	Vector point_inside() const
	{
		//static std::default_random_engine rnd(static_cast<unsigned int>(time(nullptr)));
		//static std::uniform_real_distribution<double> rng(0.0, 1.0);

		//auto e0 = rng(rnd);
		//auto e1 = rng(rnd);

		auto e0 = drand48();
		auto e1 = drand48();

		auto rootE0 = std::sqrt(e0);

		auto l0 = 1 - rootE0;
		auto l1 = e1 * rootE0;
		auto l2 = 1 - l0 - l1;

		return l0 * _p0 + l1 * _p1 + l2 * _p2;

		//return (_p0 + _p1 + _p2) / 3;
	}

	Vector center() const
	{
		return (_p0 + _p1 + _p2) / 3.0;
	}

	void init_patchs(const size_t divisions)
	{
		_subTriangles.clear();
		splitRec(divisions, *this);
	}

	const std::vector<Triangle*>& getNeighbours() const
	{
		return _neighbors;
	}

	void splitRec(size_t divisions, const Triangle& tri)
	{
		Triangle tmp[4];
		tri.split(tmp);

		for (auto i = 0; i < 4; i++)
		{
			const auto& t = tmp[i];

			if (divisions == 1)
				_subTriangles.push_back(t);
			else
				splitRec(divisions - 1, t);
		}
	}

	const Triangle& getSubTriangle(size_t index) const
	{
		return _subTriangles[index];
	}

	Triangle& getSubTriangle(size_t index)
	{
		return _subTriangles[index];
	}

	void setColor(const Color& color) { _color = color; }

	void getAdjacentTriangles(const Triangle& tri, std::vector<Triangle>& out) const
	{
		const auto epsilon = 0.00001f;
		const auto e = Vector(epsilon, epsilon, epsilon);

		for (auto it = _subTriangles.begin(); it != _subTriangles.end(); ++it)
		{
			const auto& t = *it;
			auto matches = 0;
			
			if(Vector::AreEqual(tri.getP1(), t.getP1(), e) ||
			   Vector::AreEqual(tri.getP1(), t.getP2(), e) ||
			   Vector::AreEqual(tri.getP1(), t.getP3(), e))
			{
				matches++;
			}

			if (Vector::AreEqual(tri.getP2(), t.getP1(), e) ||
				Vector::AreEqual(tri.getP2(), t.getP2(), e) ||
				Vector::AreEqual(tri.getP2(), t.getP3(), e))
			{
				matches++;
			}

			if (Vector::AreEqual(tri.getP3(), t.getP1(), e) ||
				Vector::AreEqual(tri.getP3(), t.getP2(), e) ||
				Vector::AreEqual(tri.getP3(), t.getP3(), e))
			{
				matches++;
			}

			if(matches == 2)
			{
				out.push_back(t);
			}
		}
	}

	const Color& getEmission() const { return _emission; }
	const Color& getColor() const { return _color; }
	const Vector& getNormal() const { return _normal; }
	const Vector& getP1() const { return _p0; }
	const Vector& getP2() const { return _p1; }
	const Vector& getP3() const { return _p2; }
    double getA() const { return _a; }

	const std::vector<Triangle>& getSubTriangles() const { return _subTriangles; }

private:
	Vector _p0, _p1, _p2, _normal;
	Color _emission, _color;
	std::vector<Triangle*> _neighbors;

	double _a;
	std::vector<Triangle> _subTriangles;
};

#endif