#ifndef Consts_H_included
#define Consts_H_included

/*------------------------------------------------------------------
| Scene objects are spheres; material either perfectly diffuse,
| specular (mirror reflection) or transparent (refraction/reflection)
| (DIFFuse, SPECular, REFRactive, GLOSsy, TRANslucent)
------------------------------------------------------------------*/
enum Refl_t { DIFF, SPEC, REFR, GLOS, TRAN};
enum eGeometryType
{
	eGeometryType_Sphere,
	eGeometryType_TriangleMesh
};


#endif



