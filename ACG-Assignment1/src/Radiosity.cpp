/******************************************************************
*
* Radiosity.cpp
*
* Description: This file demonstrates global illumination rendering
* based on the radiosity method. The geometry is divided into patches
* for which the form factors are determined employing Monte Carlo 
* integration. Radiosity values for the patches are computed with
* an iterative solver. The final image (i.e. radiance) is obtained 
* via tracing rays into the scene. Two output files are saved - 
* one with constant shading of patches and one with bicubic color
* interpolation.
*
* The code is extended from software by user Hole and Kevin Beason,
* released under the MIT License. 
*
* http://kagamin.net/hole/license.txt
* http://kagamin.net/hole/smallpt-license.txt
*
* Advanced Computer Graphics Proseminar WS 2015
* 
* Interactive Graphics and Simulation Group
* Institute of Computer Science
* University of Innsbruck
*
*******************************************************************/

#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS

#ifdef _WIN32
#define drand48() (((double)rand())/((double)RAND_MAX))
#endif

#include <cmath>

const double Over_M_PI = 1.0 / M_PI;

#include <iostream>
#include "Image.hpp"
#include "RendererTriangle.h"
#include "RendererRectangles.h"
#include "Scene.hpp"

int main(int argc, char **argv) 
{
	const auto width = 320;
	const auto height = 240;
	const auto samples = 4;
	const auto iterations = 40u;
	const auto divisions = 4u;
	const auto mcSamples = 10u;

	//RendererRectangles renderer(width, height, samples);
	RendererTriangles renderer(width, height, samples);
	renderer.buildScene(getScene2());

	Image img(width, height);
	Image img_interpolated(width, height);

	renderer.Render(img, img_interpolated, divisions, mcSamples, iterations);

	img.Save("image_patches.ppm");
	img_interpolated.Save("image_smooth.ppm");

	//system("image_patches.ppm");
	system("image_smooth.ppm");

	return 0;
}
