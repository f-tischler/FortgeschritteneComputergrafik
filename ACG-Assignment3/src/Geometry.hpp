#ifndef Geometry_H
#define Geometry_H

#include "Material.hpp"
#include "BoundingBox.hpp"
#include "kdtree.h"

enum eGeometryType
{
	eGeometryType_Sphere,
	eGeometryType_TriangleMesh,
};

class Geometry
{
public:
	explicit Geometry(const Material& material, eGeometryType type) 
		: _material(material), _type(type), _tree(kd_create(3)) { }

    ~Geometry() {
        kd_free(_tree);
    }

	bool Intersects(const Ray& ray, Vector& normal, float& distance)
	{
		if (!_boundingBox.intersect(ray)) return false;

		return IntersectsGeometry(ray, normal, distance);
	}

	const auto& GetMaterial() const { return _material; }
	auto GetType() const { return _type; }
    struct kdtree* GetTree() { return _tree; }

protected:
	virtual bool IntersectsGeometry(const Ray& r, Vector& normal, float& distance) = 0;

	auto& GetBoundingBox() { return _boundingBox; }

private:
	eGeometryType _type;
	Material _material;
	BoundingBox _boundingBox;
    struct kdtree* _tree;
};

#endif // Geometry_H
