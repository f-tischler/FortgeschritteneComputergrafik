#ifndef Geometry_H_included
#define Geometry_H_included

#include "Ray.hpp"
#include "BoundingBox.hpp"

class Geometry;
typedef std::shared_ptr<Geometry> SPGeometry;

class Geometry
{
public:
	Geometry(Vector emission_,
		Vector color_, Refl_t refl_, double c_ = 0.0) :
		emission(emission_),
		color(color_), refl(refl_), c(c_)
	{
	}

	virtual SPGeometry clone() = 0;

	virtual double Intersect(const Ray &ray, Vector& normal, bool culling) = 0;

	const BoundingBox<Vector>& getBB() const { return m_BB; }

	virtual size_t getType() const = 0;

	const Color& getEmmision() const { return emission; }
	const Color& getColor() const { return color; }
	Refl_t getReflectionType() const { return refl; }
    double getConstant() const { return c; }

	void setEmmision(Color v) { emission = v; }
	void setColor(Color v) { color = v; }
	void setReflectionType(Refl_t v) { refl = v; }
    void setConstant(double c_) { c = c_; }

	Color emission, color;
	Refl_t refl;
    double c; // glossiness/translucency factor

protected:
	BoundingBox<Vector> m_BB;
};



#endif