#ifndef Raytracer_H
#define Raytracer_H

#include "Image.hpp"
#include "Camera.hpp"
#include "Scene.hpp"
#include <iostream>
#include <random>

struct TracingInfo
{
	std::vector<Ray> spawnedRays;
};

struct RaytracerConfiguration
{
	unsigned int subSamplesPerPixel;
	unsigned int samplesPerSubSample;
	unsigned int dofSamples;
	float apetureSize;
};

template<class RadianceProviderType>
class Raytracer
{
public:


	explicit Raytracer(Image& image, 
					   Camera& camera, 
					   Scene& scene, 
					   RadianceProviderType& radianceProvider, 
					   const RaytracerConfiguration& configuration) : _image(image),
															 _camera(camera),
															 _scene(scene),
															 _radianceProvider(radianceProvider),
															 _config(configuration)
	{
		
	}

	void Render() 
	{
		const auto fov = _camera.GetFov();
		const auto aspect = _image.GetAspect();

		const auto u = glm::normalize(glm::cross(_camera.GetRay().GetDirection(), Vector(0, 1, 0)));
		const auto v = glm::normalize(glm::cross(u, _camera.GetRay().GetDirection()));

		const auto focalDistance = glm::length(_camera.GetLookAt() - _camera.GetPosition());
		const auto viewPlaneHalfWidth = std::tan(fov / 2.0f) * focalDistance;
		const auto viewPlaneHalfHeight = viewPlaneHalfWidth / aspect;

		const auto viewPlaneBottomLeft = _camera.GetLookAt() - v * viewPlaneHalfHeight - u * viewPlaneHalfWidth;

		const auto xIncVector = (2.0f * u * viewPlaneHalfWidth) / static_cast<float>(_image.GetWidth());
		const auto yIncVector = (2.0f * v * viewPlaneHalfHeight) / static_cast<float>(_image.GetHeight());

		const auto apeture = u * _config.apetureSize + v * _config.apetureSize;

		/* Loop over image rows */
		for (auto y = 0u; y < _image.GetHeight(); y++)
		{
			std::cout << "\rRendering (" << _config.subSamplesPerPixel * _config.samplesPerSubSample << " spp) " << (100.0 * y / (_image.GetHeight() - 1)) << "%     ";
			srand(y * y * y);

			/* Loop over row pixels */
			for (auto x = 0u; x < _image.GetWidth(); x++)
			{
				_image.SetColor(x, y, Color());

				/* 2x2 subsampling per pixel */
				for (auto sy = 0u; sy < _config.subSamplesPerPixel; sy++)
				{
					for (auto sx = 0u; sx < _config.subSamplesPerPixel; sx++)
					{
						auto accumulated_radiance = Color();

						/* Compute radiance at subpixel using multiple samples */
						for (auto s = 0u; s < _config.samplesPerSubSample; s++)
						{
							accumulated_radiance += Sample(x, y,
								sx, sy,
								viewPlaneBottomLeft,
								xIncVector,
								yIncVector,
								apeture) / static_cast<float>(_config.samplesPerSubSample);
						}

						accumulated_radiance = glm::clamp(accumulated_radiance, 0.0f, 1.0f) 
							/ static_cast<float>(_config.subSamplesPerPixel * _config.subSamplesPerPixel);

						_image.AddColor(x, y, accumulated_radiance);
					}
				}
			}
		}

		std::cout << std::endl;
	}

	const auto& GetImage() const { return _image; }
	const auto& GetCamera() const { return _camera; }
	const auto& GetScene() const { return _scene; }
	
private:
	Image& _image;
	Camera& _camera;
	Scene& _scene;
	RaytracerConfiguration _config;

	RadianceProviderType& _radianceProvider;
	Color _clearColor;

	std::default_random_engine _rnd;
	std::uniform_real_distribution<float> _rng;

	auto GetRandom() { return _rng(_rnd); }

	Color Sample(int x, int y, int, int, const Vector& viewPlaneBottomLeft, const Vector& xIncVector, const Vector& yIncVector, const Vector& apeture)
	{
		const auto r1 = 2.0f * GetRandom();
		const auto r2 = 2.0f * GetRandom();

		/* Transform uniform into non-uniform filter samples */
		auto sdx = (r1 < 1.0f) ? sqrt(r1) - 1.0f : 1.0f - sqrt(2.0f - r1);
		auto sdy = (r2 < 1.0f) ? sqrt(r2) - 1.0f : 1.0f - sqrt(2.0f - r2);

		auto radiance = Color();

		auto dofSamples = 2ul;
		for (auto i = 0ul; i < dofSamples; i++)
		{
			const auto r3 = 2.0f * (GetRandom() - 0.5f);
			const auto r4 = 2.0f * (GetRandom() - 0.5f);

			const auto blurVec = Vector(r3 * apeture.x, r4 * apeture.y, 0);

			const auto viewPlanePoint = viewPlaneBottomLeft + static_cast<float>(sdx + x) * xIncVector
															+ static_cast<float>(sdy + y) * yIncVector;

			const auto eyePoint = _camera.GetPosition() + blurVec;
			const auto dir = glm::normalize(viewPlanePoint - eyePoint);

			radiance += GetRadiance(Ray(eyePoint, dir));
		}

		return radiance / static_cast<float>(dofSamples);
	}

	Color GetRadiance(const Ray& ray) const
	{
		IntersectionInfo intersectionInfo;
		TracingInfo tracingInfo;

		auto radiance = _clearColor;

		if(_scene.Intersect(ray, intersectionInfo))
		{
			radiance = _radianceProvider.GetRadiance(intersectionInfo, tracingInfo);
			if(!tracingInfo.spawnedRays.empty())
			{
				for (const auto& subRay : tracingInfo.spawnedRays)
					radiance += GetRadiance(subRay);
			}
		}

		return radiance;
	}
};

#endif // Raytracer_H
