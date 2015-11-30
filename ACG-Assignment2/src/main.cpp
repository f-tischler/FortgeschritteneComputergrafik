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
#include "Scene.hpp"
#include "Renderer.hpp"

using namespace std;


#include "Image.hpp"
#include "Ray.hpp"
#include "Vector.hpp"
#include "Renderer.hpp"
#include "Triangle.hpp"




/******************************************************************
* Main routine: Computation of path tracing image (2x2 subpixels)
* Key parameters
* - Image dimensions: width, height 
* - Number of samples per subpixel (non-uniform filtering): samples 
* Rendered result saved as PPM image file
*******************************************************************/

int main(int argc, char *argv[]) 
{
    int width = 1024;
    int height = 768;
    int samples = 1;

	Renderer renderer;
	renderer.buildScene(getScene());
	Image img = renderer.render(width, height, samples);

	img.Save(string("image.ppm"));
	system(string("image.ppm").c_str());
}
