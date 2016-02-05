
#ifndef main_H
#define main_H

#include <cstdlib>

#include "Common.hpp"
#include "Image.hpp"
#include "Camera.hpp"
#include "Scene.hpp"
#include "Sphere.hpp"
#include "Raycaster.hpp"
#include "NanoFlannKdPMRP.hpp"
#include "SimpleRadianceProvider.hpp"
#include "SimplePhotonMappingRadianceProvider.h"
#include "NanoFlannPhotonMap.h"
#include "PathTracingWithPhotonMappingRP.hpp"
#include "ModelLoader.hpp"

int main(int argc, char* argv[])
{
	glutInit(&argc, argv);

	constexpr auto width = 640ul;
	constexpr auto height = 480ul;

	constexpr auto roomWidth  = 400.f;
	constexpr auto roomHeight = 300.f;

	const auto camPos = Vector(0, 130.0, 500);
	const auto lookAt = Vector(0, 145, 190);
	constexpr auto fov = 80.0f / 180.0f * PI;

	constexpr auto lightPower = 1000000.f;

	Scene scene;
	Image image(width, height);
	Camera camera(camPos, lookAt, fov);

	auto redDiffuseMat = Material(eReflectionType::DIFF, Color(1, 0, 0), Color(0, 0, 0), 0);
	auto blueDiffuseMat = Material(eReflectionType::DIFF, Color(0, 0, 1), Color(0, 0, 0), 0);
	auto greyDiffuseMat = Material(eReflectionType::DIFF, Color(0.1, 0.1, 0.1), Color(0, 0, 0), 0);
	auto brownDiffuseMat = Material(eReflectionType::DIFF, Color(0.823529, 0.411765, 0.117647), Color(0, 0, 0), 0);
	auto greenDiffuseMat = Material(eReflectionType::DIFF, Color(0, 1, 0), Color(0, 0, 0), 0);
	auto whiteDiffuseMat = Material(eReflectionType::DIFF, Color(1, 1, 1), Color(0, 0, 0), 0);

	auto greySpecularMat = Material(eReflectionType::SPEC, Color(0.6, 0.6, 0.6), Color(0, 0, 0), 0);
	auto greyGlossyMat = Material(eReflectionType::SPEC, Color(0.6, 0.6, 0.6), Color(0, 0, 0), 0.8f);
	auto whiteTransMat = Material(eReflectionType::TRAN, Color(1, 1, 1)*.999f, Color(0), 0.0f);

	auto lightColor = Color(255.f / 255.f / 2, 147.f / 255.f / 2, 41.f / 255.f / 2);
	auto lightMat = Material(eReflectionType::DIFF, lightColor, lightColor*lightPower, 0.0f);

	//Back left
	scene.AddGeometry(std::make_unique<Sphere>(greyGlossyMat, Vector(-70, 45, 70), 40.0f));

	//back right
	scene.AddGeometry(std::make_unique<Sphere>(whiteDiffuseMat, Vector(95, 90, 70), 40.0f));

	//front left
	scene.AddGeometry(std::make_unique<Sphere>(greySpecularMat, Vector(-100, 0, 280), 50.0f));

	// camera focus point
	//scene.AddGeometry(std::make_unique<Sphere>(whiteDiffuseMat, lookAt, 5.0f));

	//front right
	//scene.AddGeometry(std::make_unique<Sphere>(whiteTransMat, Vector(75, 60, 250), 40.0f));

	auto diamond = loadObj("diamond.obj", whiteTransMat);
	diamond->scale(Vector(30.0f, 30.0f, 30.0f));
	//diamond->rotate(Vector(20.f, 0.f, 20.f));
	diamond->translate(Vector(0, 140, 160));

	scene.AddGeometry(std::move(diamond));

	auto lightHousing = loadObj("light_housing.obj", redDiffuseMat);
	lightHousing->scale(Vector(300, 300, 300));
	lightHousing->translate(Vector(0, 210, 160));
	scene.AddGeometry(std::move(lightHousing));

	//light
	scene.AddGeometry(std::make_unique<Sphere>(lightMat, Vector(0, 210, 160), 8.0f));

	auto r = 50000.0f;

	scene.AddGeometry(std::make_unique<Sphere>(brownDiffuseMat, Vector(0, r + roomHeight, 0), r));
	scene.AddGeometry(std::make_unique<Sphere>(brownDiffuseMat, Vector(0,-r, 0), r));
	scene.AddGeometry(std::make_unique<Sphere>(greenDiffuseMat, Vector(-r - roomWidth/2, 0, 0), r));
	scene.AddGeometry(std::make_unique<Sphere>(blueDiffuseMat, Vector(r + roomWidth / 2, 0, 0), r));
	scene.AddGeometry(std::make_unique<Sphere>(greyDiffuseMat, Vector(0, 0, -r - roomWidth / 2), r));

	// add extra scope to for the radiance provider so that its memory is freed after rendering
	{
		auto radianceProvider = PathTracingWithPhotonMappingRP<NanoFlannPhotonMap>(scene, false, true);
		radianceProvider.CreatePhotonMap(scene);

		auto config = RaycasterConfiguration
		{
			2, // sub samples per pixel
			2, // samples per sub sample;
			1, // dof samples;
			4, // apeture size;
		};

		auto raycaster = Raycaster<decltype(radianceProvider)>(image, camera, radianceProvider, config);

		raycaster.Render();
	}

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
