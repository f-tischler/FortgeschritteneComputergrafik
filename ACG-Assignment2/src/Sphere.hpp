#ifndef Sphere_h_included
#define Sphere_h_included

#include "Consts.hpp"
#include "Vector.hpp"
#include "Ray.hpp"
#include "Geometry.hpp"


class Sphere : public Geometry
{
public:
	double radius;
	Vector position;

public:

	Sphere(double radius_, Vector position_, Vector emission_, Vector color_, Refl_t refl_, double c = 0.0) :
		Geometry(emission_, color_, refl_, c),
		radius(radius_), 
		position(position_)
	{
		m_BB.addPoint(position_ + radius);
		m_BB.addPoint(position_ - radius);
	}

	virtual double Intersect(const Ray &ray, Vector& normal, bool culling)
	{
		/* Check for ray-sphere intersection by solving for t:
		t^2*d.d + 2*t*(o-p).d + (o-p).(o-p) - R^2 = 0 */
		Vector op = position - ray.org;
		double b = op.Dot(ray.dir);
		double radicant = b*b - op.Dot(op) + radius*radius;
		if (radicant < 0.0)
			return 0.0;      /* No intersection */
		else
			radicant = sqrt(radicant);

		double t;
		t = b - radicant;    /* Check smaller root first */
		if (t > 1e-4)
		{
			op = ray.org + ray.dir * t;
			normal = (op - position).Normalized();  /* Normal at intersection */
			return t;
		}

		t = b + radicant;
		if (t > 1e-4)          /* Check second root */
		{
			op = ray.org + ray.dir * t;
			normal = (op - position).Normalized();  /* Normal at intersection */
			return t;
		}

		return 0.0;          /* No intersection in ray direction */
	}


	virtual SPGeometry clone()
	{
		return SPGeometry(new Sphere(radius, position, emission, color, refl, c));
	}

	virtual size_t getType() const { return eGeometryType_Sphere; }
};

#endif