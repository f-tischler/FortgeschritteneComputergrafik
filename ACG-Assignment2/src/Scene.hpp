#ifndef Scene_H_included 
#define Scene_H_included

#include <vector>

#include "Sphere.hpp"


inline std::vector<Sphere> getScene()
{
	/******************************************************************
	* Hard-coded scene definition: the geometry is composed of spheres
	* (i.e. Cornell box walls are part of very large spheres).
	* These are defined by:
	* - radius, center
	* - emitted light (light sources), surface reflectivity (~color),
	*   material
	*******************************************************************/
	static Sphere spheres[] =
	{
		Sphere(1e5, Vector(1e5 + 1, 40.8, 81.6), Vector(), Vector(.75, .25, .25), DIFF), /* Left wall */
		Sphere(1e5, Vector(-1e5 + 99, 40.8, 81.6), Vector(), Vector(.25, .25, .75), DIFF), /* Rght wall */
		Sphere(1e5, Vector(50, 40.8, 1e5), Vector(), Vector(.75, .75, .75), DIFF), /* Back wall */
		Sphere(1e5, Vector(50, 40.8, -1e5 + 170), Vector(), Vector(), DIFF), /* Front wall */
		Sphere(1e5, Vector(50, 1e5, 81.6), Vector(), Vector(.75, .75, .75), DIFF), /* Floor */
		Sphere(1e5, Vector(50, -1e5 + 81.6, 81.6), Vector(), Vector(.75, .75, .75), DIFF), /* Ceiling */

		Sphere(16.5, Vector(27, 16.5, 47), Vector(), Vector(1, 1, 1)*.999, SPEC), /* Mirror sphere */
		Sphere(16.5, Vector(73, 16.5, 78), Vector(), Vector(1, 1, 1)*.999, REFR), /* Glas sphere */

		Sphere(1.5, Vector(50, 81.6 - 16.5, 81.6), Vector(4, 4, 4) * 100, Vector(), DIFF), /* Light */
	};


	return std::vector<Sphere>(spheres, spheres + sizeof(spheres)/sizeof(spheres[0]));
}


#endif