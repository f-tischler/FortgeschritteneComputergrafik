#ifndef Geometry_H
#define Geometry_H

#include "Material.hpp"
#include "BoundingBox.hpp"

enum eGeometryType
{
	eGeometryType_Sphere,
	eGeometryType_TriangleMesh,
};

class Geometry
{
public:
	explicit Geometry(const Material& material, eGeometryType type) 
		: _type(type), _material(material) { }

	virtual ~Geometry() = default;

	bool Intersects(const Ray& ray, Vector& normal, float& distance)
	{
		if (!_boundingBox.intersect(ray)) return false;

		return IntersectsGeometry(ray, normal, distance);
	}

	const auto& GetMaterial() const { return _material; }
	auto GetType() const { return _type; }

protected:
	virtual bool IntersectsGeometry(const Ray& r, Vector& normal, float& distance) = 0;

	auto& GetBoundingBox() { return _boundingBox; }

private:
	eGeometryType _type;
	Material _material;
	BoundingBox _boundingBox;
};

#endif // Geometry_H
