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
	constexpr auto width = 640ul;
	constexpr auto height = 480ul;

	const auto camPos = Vector(0, 120.0, 600);
	const auto lookAt = Vector(0, 10, 0);
	constexpr auto fov = 1.0472f;//60.0f / 180.0f * PI;

	Scene scene;
	Image image(width, height);
	Camera camera(camPos, lookAt, fov);

	auto redDiffuseMat = Material(eReflectionType::DIFF, Color(1, 0, 0), Color(0, 0, 0), 0);
	auto blueDiffuseMat = Material(eReflectionType::DIFF, Color(0, 0, 1), Color(0, 0, 0), 0);
	auto greyDiffuseMat = Material(eReflectionType::DIFF, Color(0.5, 0.5, 0.5), Color(0, 0, 0), 0);
	auto greenDiffuseMat = Material(eReflectionType::DIFF, Color(0, 1, 0), Color(0, 0, 0), 0);
	auto whiteDiffuseMat = Material(eReflectionType::DIFF, Color(1, 1, 1), Color(0, 0, 0), 0);

	auto whiteSpecularMat = Material(eReflectionType::SPEC, Color(1, 1, 1), Color(0, 0, 0), 0);
	auto greyGlossyMat = Material(eReflectionType::GLOS, Color(0.6, 0.6, 0.6), Color(0, 0, 0), 0.5f);
	auto whiteTransMat = Material(eReflectionType::TRAN, Color(1,1,1), Color(0, 0, 0), 0.0f);

	auto lightMat = Material(eReflectionType::DIFF, Color(1, 1, 1), Color(500000, 500000, 500000), 0.0f);


	//Back left
	scene.AddGeometry(std::make_unique<Sphere>(greyGlossyMat, Vector(-70, 45, 70), 40.0f));

	//back right
	scene.AddGeometry(std::make_unique<Sphere>(whiteSpecularMat, Vector(75, 45, 70), 40.0f));

	//front left
	scene.AddGeometry(std::make_unique<Sphere>(whiteTransMat, Vector(-70, 45, 190), 70.0f));

	//front right
	scene.AddGeometry(std::make_unique<Sphere>(blueDiffuseMat, Vector(75, 35, 190), 20.0f));

	//light
	scene.AddGeometry(std::make_unique<Sphere>(lightMat, Vector(0, 200, 60), 8.0f));

	scene.AddGeometry(std::make_unique<Sphere>(whiteDiffuseMat, Vector(0, -5000, 0), 5000.0f));
	scene.AddGeometry(std::make_unique<Sphere>(whiteDiffuseMat, Vector(0, -5500, 0), 5000.0f));
	scene.AddGeometry(std::make_unique<Sphere>(greenDiffuseMat, Vector(-5150, 0, 0), 5000.0f));
	scene.AddGeometry(std::make_unique<Sphere>(blueDiffuseMat, Vector(5150, 0, 0), 5000.0f));
	scene.AddGeometry(std::make_unique<Sphere>(greyDiffuseMat, Vector(0, 0, -5050), 5000.0f));

	auto radianceProvider = SimplePhotonMappingRadianceProvider(false);
	radianceProvider.CreatePhotonMap(scene);

	auto config = RaytracerConfiguration
	{
		1, // subSamplesPerPixel
		1, // unsigned int samplesPerSubSample;
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
