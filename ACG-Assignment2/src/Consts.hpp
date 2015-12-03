#ifndef Consts_H_included
#define Consts_H_included

/*------------------------------------------------------------------
| Scene objects are spheres; material either perfectly diffuse,
| specular (mirror reflection) or transparent (refraction/reflection)
| (DIFFuse, SPECular, REFRactive)
------------------------------------------------------------------*/
enum Refl_t { DIFF, SPEC, REFR };
enum eGeometryType
{
	eGeometryType_Sphere,
	eGeometryType_TriangleMesh
};


#endif



