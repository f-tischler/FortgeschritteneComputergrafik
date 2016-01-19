#ifndef Sphere_H
#define Sphere_H

#include "Geometry.hpp"

class Sphere : public Geometry
{
public:
	explicit Sphere(const Material& material, const Vector& position, float radius)
		: Geometry(material, eGeometryType_Sphere), _radius(radius), _position(position)
	{
		GetBoundingBox().addPoint(_position + _radius);
		GetBoundingBox().addPoint(_position - _radius);
	}

	const auto& GetPosition() const { return _position; }
	auto GetRadius() const { return _radius; }

protected:
	bool IntersectsGeometry(const Ray& ray, Vector& normal, float& distance) override
	{
		/* Check for ray-sphere intersection by solving for t:
		t^2*d.d + 2*t*(o-p).d + (o-p).(o-p) - R^2 = 0 */
		auto op = _position - ray.GetOrigin();
		auto b = glm::dot(op, ray.GetDirection());
		auto radicant = b*b - glm::dot(op, op) + _radius * _radius;
		
		if (radicant < 0.0) return false;  /* No intersection */
		
		radicant = sqrt(radicant);

		auto t = b - radicant;    /* Check smaller root first */
		if (t > 1e-4)
		{
			op = ray.GetOrigin() + ray.GetDirection() * t;
			
			normal = glm::normalize(op - _position);  /* Normal at intersection */
			distance = t;

			return true;
		}

		t = b + radicant;
		if (t > 1e-4)          /* Check second root */
		{
			op = ray.GetOrigin() + ray.GetDirection() * t;

			normal = glm::normalize(op - _position);  /* Normal at intersection */
			distance = t;

			return true;
		}

		return false;          /* No intersection in ray direction */
	}

private:
	float _radius;
	Vector _position;
};

#endif // Sphere_H
