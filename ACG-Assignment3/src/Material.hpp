#ifndef Material_H
#define Material_H

#include "Common.hpp"

/*------------------------------------------------------------------
| Scene objects are spheres; material either perfectly diffuse,
| specular (mirror reflection) or transparent (refraction/reflection)
| (DIFFuse, SPECular, REFRactive, GLOSsy, TRANslucent)
------------------------------------------------------------------*/
enum eReflectionType { DIFF, SPEC, TRAN };

class Material
{
public:
	Material(eReflectionType reflectionType, 
			 const Color& color,
			 const Color& emission,
			 const double glossiness)
		: _emission(emission),
		  _color(color),
		  _refl(reflectionType),
		  _glossiness(glossiness) { }

	const auto& GetEmission() const { return _emission; }
	const auto& GetColor() const { return _color; }
	auto GetReflectionType() const { return _refl; }
	auto GetGlossiness() const { return _glossiness; }
	auto HasEmission() const { return _emission != Vector(0, 0, 0); }

	void SetEmmision(const Color& v) { _emission = v; }
	void SetColor(const Color& v) { _color = v; }
	void SetReflectionType(eReflectionType v) { _refl = v; }
	void SetGlossiness(double glossiness_) { _glossiness = glossiness_; }

private:
	Color _emission, _color;
	eReflectionType _refl;
	double _glossiness;
};

#endif // Material_H
