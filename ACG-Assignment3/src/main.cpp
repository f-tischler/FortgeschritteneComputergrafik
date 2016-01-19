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
	constexpr auto height = 320ul;

	const auto camPos = Vector(50.0, 200.0, 400);
	const auto lookAt = Vector(0, 0, 0);
	constexpr auto fov = 60.0f / 180.0f * PI;

	Scene scene;
	Image image(width, height);
	Camera camera(camPos, lookAt, fov);

	auto redDiffuseMat = Material(eReflectionType::DIFF, Color(1, 0, 0), Color(0, 0, 0), 0);
	auto blueDiffuseMat = Material(eReflectionType::DIFF, Color(0, 0, 1), Color(0, 0, 0), 0);
	auto greyDiffuseMat = Material(eReflectionType::DIFF, Color(0.5, 0.5, 0.5), Color(0, 0, 0), 0);
	auto greenDiffuseMat = Material(eReflectionType::DIFF, Color(0, 1, 0), Color(0, 0, 0), 0);

	auto lightMat = Material(eReflectionType::DIFF, Color(1, 1, 1), Color(5, 5, 5), 5.0f);

	scene.AddGeometry(std::make_unique<Sphere>(redDiffuseMat, Vector(0, 0, 0), 50.0f));
	scene.AddGeometry(std::make_unique<Sphere>(blueDiffuseMat, Vector(0, 0, 100), 20.0f));
	scene.AddGeometry(std::make_unique<Sphere>(lightMat, Vector(60, 80, 10), 5.0f));

	scene.AddGeometry(std::make_unique<Sphere>(greyDiffuseMat, Vector(0, -5000, 0), 5000.0f));
	scene.AddGeometry(std::make_unique<Sphere>(greenDiffuseMat, Vector(-5100, 0, 0), 5000.0f));
	scene.AddGeometry(std::make_unique<Sphere>(greyDiffuseMat, Vector(0, 0, -5100), 5000.0f));

	auto radianceProvider = SimplePhotonMappingRadianceProvider(false);
	radianceProvider.CreatePhotonMap(scene);

	auto config = RaytracerConfiguration
	{
		2, // subSamplesPerPixel
		2, // unsigned int samplesPerSubSample;
		2, // unsigned int dofSamples;
		1, // float apetureSize;
	};

	auto raytracer = Raytracer<decltype(radianceProvider)>(image,	
														   camera, 
														   scene, 
														   radianceProvider, 
														   config);

	raytracer.Render();

	image.Save("image.ppm");

	system("image.ppm");

	return EXIT_SUCCESS;
}

#endif // main_H
