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


inline SPScene getScene()
{
	static SPScene scene = SPScene(new Scene());

	scene->addGeometry(SPGeometry(new Sphere(1e5, Vector(1e5 + 1, 40.8, 81.6), Vector(), Vector(.75, .25, .25), DIFF)));	/* Left wall */
	scene->addGeometry(SPGeometry(new Sphere(1e5, Vector(-1e5 + 99, 40.8, 81.6), Vector(), Vector(.25, .25, .75), DIFF)));	/* Right wall */

	scene->addGeometry(SPGeometry(new Sphere(1e5, Vector(50, 40.8, 1e5), Vector(), Vector(.75, .75, .75), DIFF)));			/* Back wall */
	scene->addGeometry(SPGeometry(new Sphere(1e5, Vector(50, 40.8, -1e5 + 170), Vector(), Vector(), DIFF)));				/* Front wall */

	scene->addGeometry(SPGeometry(new Sphere(1e5, Vector(50, 1e5, 81.6), Vector(), Vector(.75, .75, .75), DIFF)));			/* Floor */
	scene->addGeometry(SPGeometry(new Sphere(1e5, Vector(50, -1e5 + 81.6, 81.6), Vector(), Vector(.75, .75, .75), DIFF)));  /* Ceiling */

 	scene->addGeometry(SPGeometry(new Sphere(16.5, Vector(27, 16.5, 47), Vector(), Vector(1, 1, 1)*.999, SPEC)));			/* Mirror sphere */
	scene->addGeometry(SPGeometry(new Sphere(16.5, Vector(73, 16.5, 47), Vector(), Vector(1, 1, 1)*.999, REFR)));			/* Glas sphere */

//	scene->addGeometry(SPGeometry(new Sphere(1.5, Vector(50, 81.6 - 16.5, 81.6), Vector(4, 4, 4) * 100, Vector(), DIFF))); /* Light */

	auto cubeMesh = loadObj("cube.obj", Vector(), Vector(), SPEC)[0];
	auto rhino = loadObj("cube.obj", Vector(), Vector(), SPEC)[0];

	auto cube1 = cubeMesh;
	auto cube2 = rhino;
	auto cube3 = std::static_pointer_cast<TriangleMesh>(cubeMesh->clone());
	auto cube4 = std::static_pointer_cast<TriangleMesh>(cubeMesh->clone());

	//Mirrow
	cube1->setEmmision(Vector());
	cube1->setColor(Vector(1, 1, 1)*.999);
	cube1->setReflectionType(SPEC);

	cube1->scaleToWidth(10.0);
	cube1->rotateY(45);
	cube1->translate(20, 10, 80);
	cube1->buildCollision();

	//Diffuse
	cube2->setEmmision(Vector());
	cube2->setColor(Vector(.75, .25, .25));
	cube2->setReflectionType(DIFF);

	cube2->scaleToWidth(10.0);
	cube2->rotateY(15);
	cube2->rotateX(10);
	cube2->translate(50, 10, 80);
	cube2->buildCollision();


	//Glas
	cube3->setColor(Vector(1, 1, 1)*.999);
	cube3->setEmmision(Vector());
	cube3->setReflectionType(REFR);

	cube3->scaleToWidth(10.0);
	cube3->rotateY(45);
	cube3->rotateX(45);
	cube3->translate(80, 10, 80);
	cube3->buildCollision();

	//Light
	cube4->setColor(Vector());
	cube4->setEmmision(Vector(4,4,4) * 20);
	cube4->setReflectionType(DIFF);

	cube4->scaleToWidth(5.0);
	cube4->rotateY(45);
	cube4->rotateX(45);
	cube4->translate(40, 50, 80);
	cube4->buildCollision();


	scene->addGeometry(cube1);
	scene->addGeometry(cube2);
	scene->addGeometry(cube3);
	scene->addGeometry(cube4);

	return scene;
}

inline SPScene getScene2()
{
    static SPScene scene = SPScene(new Scene());
    
    scene->addGeometry(SPGeometry(new Sphere(1e5, Vector(1e5 + 1, 40.8, 81.6), Vector(), Vector(.75, .25, .25), DIFF)));	/* Left wall */
    scene->addGeometry(SPGeometry(new Sphere(1e5, Vector(-1e5 + 99, 40.8, 81.6), Vector(), Vector(.25, .25, .75), DIFF)));	/* Right wall */
    
    scene->addGeometry(SPGeometry(new Sphere(1e5, Vector(50, 40.8, 1e5), Vector(), Vector(.75, .75, .75), DIFF)));			/* Back wall */
    scene->addGeometry(SPGeometry(new Sphere(1e5, Vector(50, 40.8, -1e5 + 170), Vector(), Vector(), DIFF)));				/* Front wall */
    
    scene->addGeometry(SPGeometry(new Sphere(1e5, Vector(50, 1e5, 81.6), Vector(), Vector(.75, .75, .75), DIFF)));			/* Floor */
    scene->addGeometry(SPGeometry(new Sphere(1e5, Vector(50, -1e5 + 81.6, 81.6), Vector(), Vector(.75, .75, .75), DIFF)));  /* Ceiling */
    
    scene->addGeometry(SPGeometry(new Sphere(16.5, Vector(27, 16.5, 47), Vector(), Vector(1, 1, 1)*.999, SPEC)));			/* Mirror sphere */
    scene->addGeometry(SPGeometry(new Sphere(16.5, Vector(73, 16.5, 47), Vector(), Vector(1, 1, 1)*.999, REFR)));			/* Glas sphere */
    
    scene->addGeometry(SPGeometry(new Sphere(1.5, Vector(50, 81.6 - 16.5, 81.6), Vector(4, 4, 4) * 100, Vector(), DIFF))); /* Light */
    
    return scene;
}


/******************************************************************
* Main routine: Computation of path tracing image (2x2 subpixels)
* Key parameters
* - Image dimensions: width, height 
* - Number of samples per subpixel (non-uniform filtering): samples 
* Rendered result saved as PPM image file
*******************************************************************/
#include <chrono>

int main(int argc, char *argv[]) 
{
    int width = 320;
    int height = 240;
	int samples = 8;


	auto start = std::chrono::system_clock::now();


	Renderer renderer;
	renderer.buildScene(getScene());
	Image img = renderer.render(width, height, samples);

	img.Save(string("image.ppm"));


	auto end = std::chrono::system_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
	std::cout << "Time elapsed:" << elapsed.count() << '\n';

	system(string("image.ppm").c_str());
}
