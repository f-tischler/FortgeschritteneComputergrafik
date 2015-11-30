#ifndef Ray_H_included
#define Ray_H_included

#include "Vector.hpp"

/*------------------------------------------------------------------
| Structure for rays (e.g. viewing ray, ray tracing)
------------------------------------------------------------------*/
struct Ray
{
	Vector org, dir;    /* Origin and direction */
	Ray(const Vector org_, const Vector &dir_) : org(org_), dir(dir_) {}
};


#endif