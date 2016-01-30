#include "KdPhotonMappingRadianceProvider.hpp"
#ifndef main_H
#define main_H

#include <cstdlib>

#include "Common.hpp"
#include "Image.hpp"
#include "Camera.hpp"
#include "Scene.hpp"
#include "Sphere.hpp"
#include "Raytracer.hpp"
#include "SimpleRadianceProvider.hpp"
#include "SimplePhotonMappingRadianceProvider.h"

int main(int argc, char* argv[])
{
	glutInit(&argc, argv);

	constexpr auto width = 400ul;
	constexpr auto height = 200ul;

	constexpr auto roomWidth  = 400.f;
	constexpr auto roomHeight = 300.f;

	const auto camPos = Vector(0, 130.0, 500);
	const auto lookAt = Vector(0, 40, 0);
	constexpr auto fov = 80.0f / 180.0f * PI;

	constexpr auto lightPower = 5000000.f;

	Scene scene;
	Image image(width, height);
	Camera camera(camPos, lookAt, fov);

	auto redDiffuseMat = Material(eReflectionType::DIFF, Color(1, 0, 0), Color(0, 0, 0), 0);
	auto blueDiffuseMat = Material(eReflectionType::DIFF, Color(0, 0, 1), Color(0, 0, 0), 0);
	auto greyDiffuseMat = Material(eReflectionType::DIFF, Color(0.5, 0.5, 0.5), Color(0, 0, 0), 0);
	auto greenDiffuseMat = Material(eReflectionType::DIFF, Color(0, 1, 0), Color(0, 0, 0), 0);
	auto whiteDiffuseMat = Material(eReflectionType::DIFF, Color(1, 1, 1), Color(0, 0, 0), 0);

	auto greySpecularMat = Material(eReflectionType::SPEC, Color(0.6, 0.6, 0.6), Color(0, 0, 0), 0);
	auto greyGlossyMat = Material(eReflectionType::SPEC, Color(0.6, 0.6, 0.6), Color(0, 0, 0), 0.8f);
	auto whiteTransMat = Material(eReflectionType::TRAN, Color(1, 1, 1)*.999f, Color(0), 0.0f);

	auto lightMat = Material(eReflectionType::DIFF, Color(255.f /255.f, 147.f / 255.f, 41.f / 255.f), Color(255.f / 255.f, 147.f / 255.f, 41.f / 255.f)*lightPower, 0.0f);


	//Back left
	scene.AddGeometry(std::make_unique<Sphere>(greyGlossyMat, Vector(-70, 45, 70), 40.0f));

	//back right
	scene.AddGeometry(std::make_unique<Sphere>(redDiffuseMat, Vector(75, 90, 70), 40.0f));

	//front left
	scene.AddGeometry(std::make_unique<Sphere>(greySpecularMat, Vector(-70, 0, 250), 60.0f));

	//front right
	scene.AddGeometry(std::make_unique<Sphere>(whiteTransMat, Vector(75, 50, 250), 40.0f));

	//light
	scene.AddGeometry(std::make_unique<Sphere>(lightMat, Vector(0, 190, 160), 8.0f));

	auto r = 50000.0f;

	scene.AddGeometry(std::make_unique<Sphere>(greyDiffuseMat, Vector(0, r + roomHeight, 0), r));
	scene.AddGeometry(std::make_unique<Sphere>(greyDiffuseMat, Vector(0,-r, 0), r));
	scene.AddGeometry(std::make_unique<Sphere>(greenDiffuseMat, Vector(-r - roomWidth/2, 0, 0), r));
	scene.AddGeometry(std::make_unique<Sphere>(blueDiffuseMat, Vector(r + roomWidth / 2, 0, 0), r));
	scene.AddGeometry(std::make_unique<Sphere>(whiteDiffuseMat, Vector(0, 0, -r - roomWidth / 2), r));

	auto radianceProvider = SimplePhotonMappingRadianceProvider(false);
	radianceProvider.CreatePhotonMap(scene);

	auto config = RaytracerConfiguration
	{
		2, // subSamplesPerPixel
		2, // unsigned int samplesPerSubSample;
		1, // unsigned int dofSamples;
		1, // float apetureSize;
	};

	auto raytracer = Raytracer<decltype(radianceProvider)>(image,	
														   camera, 
														   scene, 
														   radianceProvider, 
														   config);

	raytracer.Render();

    auto filename = "image.ppm";
	image.Save(filename);

#if defined(_WIN32) || defined(_WIN64)
    system(filename);
#else
    system(std::string("open " + std::string(filename)).c_str());
#endif

	return EXIT_SUCCESS;
}

#endif // main_H
