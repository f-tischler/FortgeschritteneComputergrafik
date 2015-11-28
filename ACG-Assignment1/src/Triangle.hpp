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

		_normal = edgeACrossEdgeB.Normalized();
		_a = edgeACrossEdgeB.Length();
	}

	double area() const { return _a; }
	Vector normal() const { return _normal; }

	bool intersectsSub(const Ray &ray, double& distance, unsigned int& subIndex) const
	{
		for (unsigned int i = 0; i < _subTriangles.size(); i++)
		{
			double tmp;
			if (_subTriangles[i].intersects(ray, tmp))
			{
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

	unsigned int getSubTriangleCount() { return _subTriangles.size(); }

	void split(Triangle(&out)[4]) const
	{
		auto midP1P0 = _p0 + (_p1 - _p0) / 2;
		auto midP2P0 = _p0 + (_p2 - _p0) / 2;
		auto midP2P1 = _p1 + (_p2 - _p1) / 2;

		out[0] = Triangle(midP1P0, midP2P1, midP2P0, _emission, _color);
		out[1] = Triangle(_p0, midP1P0, midP2P0, _emission, _color);
		out[2] = Triangle(midP1P0, _p1, midP2P1, _emission, _color);
		out[3] = Triangle(midP2P1, _p2, midP2P0, _emission, _color);

// 
// 		//For testing interpolation
// 		out[0] = Triangle(midP1P0, midP2P1, midP2P0, _emission, _color*0.2);
// 		out[1] = Triangle(_p0, midP1P0, midP2P0, _emission, _color*0.4);
// 		out[2] = Triangle(midP1P0, _p1, midP2P1, _emission, _color*0.7);
// 		out[3] = Triangle(midP2P1, _p2, midP2P0, _emission, _color);
	}

	Vector point_inside() const
	{
		return (_p0 + _p1 + _p2) / 3;
	}

	void init_patchs(const int divisions)
	{
		_subTriangles.clear();
		splitRec(divisions, *this);
	}

	void getNeighbours(unsigned int subIndex, Triangle* neighbours) const
	{
		int pos = subIndex % 4;
		int startIndex = subIndex - pos;

		neighbours[0] = getSubTriangle(startIndex);
		neighbours[1] = getSubTriangle(startIndex+1);
		neighbours[2] = getSubTriangle(startIndex+2);
		neighbours[3] = getSubTriangle(startIndex+3);
	}

	void splitRec(int divisions, const Triangle& tri)
	{
		Triangle tmp[4];
		tri.split(tmp);

		for (auto tri : tmp)
		{
			if (divisions == 1)
				_subTriangles.push_back(tri);
			else
				splitRec(divisions - 1, tri);
		}
	}

	const Triangle& getSubTriangle(int index) const
	{
		return _subTriangles[index];
	}

	Triangle& getSubTriangle(int index) 
	{
		return _subTriangles[index];
	}

	void setColor(const Color& color) { _color = color; }

	const Color& getEmission() const { return _emission; }
	const Color& getColor() const { return _color; }
	const Vector& getNormal() const { return _normal; }
	const Vector& getP1() const { return _p0; }
	const Vector& getP2() const { return _p1; }
	const Vector& getP3() const { return _p2; }

private:
	Vector _p0, _p1, _p2, _normal;
	Color _emission, _color;

	double _a;
	std::vector<Triangle> _subTriangles;
};

#endif