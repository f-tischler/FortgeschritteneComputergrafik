/******************************************************************
*
* PathTracing.cpp
*
* Description: This program demonstrates global illumination rendering
* based on the path tracing method. The intergral in the rendering
* equation is approximated via Monte-Carlo integration; explicit 
* direct lighting is included to improve quality; the rendered image 
* is saved in PPM format.
*
* The code is largely based on the software smallpt by Kevin Beason,
* released under the MIT License.
*
* Advanced Computer Graphics Proseminar WS 2015
* 
* Interactive Graphics and Simulation Group
* Institute of Computer Science
* University of Innsbruck
*
*******************************************************************/

#if defined(_WIN32)
	#define _CRT_SECURE_NO_WARNINGS
	#define _USE_MATH_DEFINES

	#define drand48() (((double)rand())/((double)RAND_MAX))
#endif


/* Standard includes */
#include <cmath>   
#include <cstdlib> 
#include <iostream>
#include <fstream>

using namespace std;

#include "Image.hpp"
#include "Ray.hpp"
#include "Vector.hpp"
#include "Renderer.hpp"
#include "Triangle.hpp"
#include "Scene.hpp"


inline void cornellBox(SPScene scene)
{
	scene->addGeometry(SPGeometry(new Sphere(1e5, Vector(1e5 + 1, 40.8, 81.6), Vector(), Vector(.75, .25, .25), DIFF)));	/* Left wall */
	scene->addGeometry(SPGeometry(new Sphere(1e5, Vector(-1e5 + 99, 40.8, 81.6), Vector(), Vector(.25, .25, .75), DIFF)));	/* Right wall */

	scene->addGeometry(SPGeometry(new Sphere(1e5, Vector(50, 40.8, 1e5), Vector(), Vector(.75, .75, .75), DIFF)));			/* Back wall */
	//scene->addGeometry(SPGeometry(new Sphere(1e5, Vector(50, 40.8, -1e5 + 170), Vector(), Vector(), DIFF)));				/* Front wall */

	scene->addGeometry(SPGeometry(new Sphere(1e5, Vector(50, 1e5, 81.6), Vector(), Vector(.75, .75, .75), DIFF)));			/* Floor */
	scene->addGeometry(SPGeometry(new Sphere(1e5, Vector(50, -1e5 + 81.6, 81.6), Vector(), Vector(.75, .75, .75), DIFF)));  /* Ceiling */
}

inline void lightSphere(SPScene scene)
{
	scene->addGeometry(SPGeometry(new Sphere(1.5, Vector(50, 81.6 - 16.5, 81.6), Vector(4, 4, 4) * 100, Vector(), DIFF))); /* Light */
}

inline void lightCube(SPScene scene)
{
	auto cubeMesh = loadObj("cube.obj", Vector(), Vector(), SPEC)[0];

	//Light
	cubeMesh->setColor(Vector());
	cubeMesh->setEmmision(Vector(4, 4, 4) * 20);
	cubeMesh->setReflectionType(DIFF);

	cubeMesh->scaleToWidth(5.0);
	cubeMesh->rotateY(45);
	cubeMesh->rotateX(45);
	cubeMesh->translate(40, 50, 80);
	cubeMesh->buildCollision();
	
	scene->addGeometry(cubeMesh); /* Light */
}

inline void spheres(SPScene scene)
{
	scene->addGeometry(SPGeometry(new Sphere(16.5, Vector(27, 50, 47), Vector(), Vector(1, 1, 1)*.999, GLOS, 0.5)));			/* Glossy sphere */
	scene->addGeometry(SPGeometry(new Sphere(16.5, Vector(73, 50, 47), Vector(), Vector(1, 1, 1)*.999, TRAN, 0.5)));			/* Translucent sphere */
}

inline void cornellSpheres(SPScene scene)
{
	scene->addGeometry(SPGeometry(new Sphere(16.5, Vector(27, 16.5, 47), Vector(), Vector(1, 1, 1)*.999, GLOS, 0.5)));			/* Glossy sphere */
	scene->addGeometry(SPGeometry(new Sphere(16.5, Vector(73, 16.5, 78), Vector(), Vector(1, 1, 1)*.999, TRAN, 0.5)));			/* Translucent sphere */
}

inline void cornellScene(SPScene scene)
{
	cornellBox(scene);
	spheres(scene);
	lightSphere(scene);

	//cornellSpheres(scene);
}

inline void finalScene(SPScene scene)
{
	cornellBox(scene);
	lightSphere(scene);
	spheres(scene);


	auto mesh = loadObj("cube.obj", Vector(), Vector(), GLOS)[0];

	//Mirrow
	mesh->setEmmision(Vector());
	mesh->setColor(Vector(1, 1, 1)*.999);
	mesh->setReflectionType(GLOS);

	mesh->scaleToWidth(20.0);
	mesh->rotateX(10);
	mesh->translate(75, 10, 85);
	mesh->buildCollision();
	scene->addGeometry(mesh);


	mesh = loadObj("cube.obj", Vector(), Vector(), REFR)[0];

	//Mirrow
	mesh->setEmmision(Vector());
	mesh->setColor(Vector(1, 1, 1)*.999);
	mesh->setReflectionType(REFR);

	mesh->scaleToWidth(15.0);
	mesh->rotateY(45);
	mesh->rotateX(45);
	mesh->translate(10, 10, 85);
	mesh->buildCollision();
	scene->addGeometry(mesh);
}


/******************************************************************
* Main routine: Computation of path tracing image (2x2 subpixels)
* Key parameters
* - Image dimensions: width, height 
* - Number of samples per subpixel (non-uniform filtering): samples 
* Rendered result saved as PPM image file
*******************************************************************/
int main(int argc, char* argv[]) 
{
	auto width = 640;
    auto height = 480;
	auto samples = 8;
	auto scene = std::make_shared<Scene>();
	auto apetureSize = 5.0;

	if (argc > 1)
		width = atoi(argv[1]);
	if (argc > 2)
		height = atoi(argv[2]);
	if (argc > 3)
		samples = atoi(argv[3]);
	if (argc > 4)
	{
		cornellBox(scene);
		lightSphere(scene);
		auto mesh = loadObj(argv[4], Vector(), Vector(), SPEC)[0];

		//Mirror
		mesh->setEmmision(Vector());
		mesh->setColor(Vector(1, 1, 1)*.999);
		mesh->setReflectionType(SPEC);

		mesh->scaleToWidth(30.0);
		mesh->translate(50, 10, 50);
		mesh->buildCollision();

		scene->addGeometry(mesh);
	}
	else
	{
		finalScene(scene);
	}

	const auto camPos = Vector(50.0, 52.0, 220);
	const auto lookAt = Vector(75, 10, 105);

	const Camera cam(camPos, lookAt);
	Renderer renderer(cam);

	// uncomment to visualize focus point
	//scene->addGeometry(std::make_shared<Sphere>(1.5, lookAt, Vector(), Vector(1, 0 , 0), DIFF));

	renderer.buildScene(scene);

	auto img = renderer.render(width, height, samples, apetureSize);

	auto filename = "image.ppm";
	img.Save(filename);

#if defined(_WIN32) || defined(_WIN64)
	system(filename);
#else
	system(string("open " + std::string(filename)).c_str());
#endif
}
