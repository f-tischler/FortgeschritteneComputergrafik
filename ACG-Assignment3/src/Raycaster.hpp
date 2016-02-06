#ifndef Raycaster_H
#define Raycaster_H

#include "Image.hpp"
#include "Camera.hpp"
#include "Scene.hpp"
#include "Material.hpp"

#include <iostream>
#include <random>
#include <future>

#include "LiveImage.hpp"

#if defined(_WIN32)
	#define _CRT_SECURE_NO_WARNINGS
	#define _USE_MATH_DEFINES

	#define drand48() (((float)rand())/((float)RAND_MAX))
#endif

class TracingInfo
{
public:
	using ColorAdaptorFunc = std::function<Color(const Color&)>;
	
	struct SubRay
	{
		Ray ray;
		ColorAdaptorFunc adaptor;
	};

	TracingInfo()
		: _adaptor([](const Color& color) { return color; }),
		  _radiance(Color(0, 0, 0)), 
		  _depth(0) { }

	explicit TracingInfo(const ColorAdaptorFunc& adaptor)
		: _adaptor(adaptor),
		_radiance(Color(0, 0, 0)),
		_depth(0) { }

	void SpawnRay(const Ray& ray)
	{
		_subRays.push_back({ ray, [](const auto& c) { return c; } });
	}

	void SpawnRay(const Ray& ray, const ColorAdaptorFunc& adaptor)
	{
		_subRays.push_back({ ray, adaptor });
	}

	void SetRadiance(const Color& radiance) { _radiance = radiance; }
	void AddRadiance(const Color& radiance) { _radiance += radiance; }

	const auto& GetSubRays() const { return _subRays; }
	const auto& GetRadiance() const { return _radiance; }
	
	auto GetAdaptedRadiance() const { return _adaptor(_radiance); }

	auto GetSubRayTracingInfo()
	{
		auto info = *this;

		info._subRays.clear();
		info._depth++;

		return info;
	}

private:
	ColorAdaptorFunc _adaptor;
	std::vector<SubRay> _subRays;
	Color _radiance;
	size_t _depth;
};

struct RaycasterConfiguration
{
	unsigned int subSamplesPerPixel;
	unsigned int samplesPerSubSample;
	unsigned int dofSamples;
	float apetureSize;
};

template<class RadianceProviderType>
class Raycaster
{
public:
	explicit Raycaster(Image& image, 
					   Camera& camera,  
					   RadianceProviderType& radianceProvider, 
					   const RaycasterConfiguration& configuration) 
		: _image(image),
		  _camera(camera),
		  _config(configuration),
		  _radianceProvider(radianceProvider)
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

		const auto numThreads = std::thread::hardware_concurrency();
		const auto numXperThread = static_cast<Image::SizeType>(_image.GetWidth() / numThreads + (_image.GetWidth() % numThreads == 0 ? 0 : 1));
		//const auto numYperThread = static_cast<Image::SizeType>(_image.GetHeight() / numThreads + (_image.GetHeight() % numThreads == 0 ? 0 : 1));

		LiveImage liveImage(_image.GetWidth(), _image.GetHeight());
		liveImage.Show();

		std::atomic<Image::SizeType> fragmentsDone(0);
		std::vector<std::future<void>> tasks;

		for (auto tidx = 0ul; tidx < numThreads; tidx++)
		{
			const auto startX = tidx * numXperThread;
			const auto startY = 0u;

			auto endX = startX + numXperThread;
			if (endX >= _image.GetWidth()) endX = _image.GetWidth();

			auto endY = _image.GetHeight();

			tasks.push_back(std::async(/*std::launch::deferred,*/ [&, this, startY, startX, endY, endX, tidx]()
			{
				for (auto y = startY; y < endY; y++)
				{
					for (auto x = startX; x < endX; x++)
					{
						std::vector<std::future<Color>> futures;

						for (auto sy = 0u; sy < _config.subSamplesPerPixel; sy++)
						{
							for (auto sx = 0u; sx < _config.subSamplesPerPixel; sx++)
							{
								/* Compute radiance at subpixel using multiple samples */
								for (auto s = 0u; s < _config.samplesPerSubSample; s++)
								{
									futures.push_back(std::async(/*std::launch::deferred, */[&, this]()
									{
										return Sample(x, y,
											sx, sy,
											viewPlaneBottomLeft,
											xIncVector,
											yIncVector,
											apeture) / static_cast<float>(_config.samplesPerSubSample);
									}));
								}
							}
						}

						_image.SetColor(x, y, Color());

						auto accumulated_radiance = Color();
						for (auto& future : futures)
						{
							accumulated_radiance += glm::clamp(future.get(), 0.0f, 1.0f)
								/ static_cast<float>(_config.subSamplesPerPixel * _config.subSamplesPerPixel);
						}

						_image.AddColor(x, y, accumulated_radiance);

						liveImage.Set(x, y, accumulated_radiance);
					}
					
					fragmentsDone += endX - startX;
					if (tidx == 0)
					{
						std::cout << "\rRendering (" << _config.subSamplesPerPixel * _config.samplesPerSubSample * _config.samplesPerSubSample << " spp) " << (100.0 * fragmentsDone / (_image.GetHeight() * _image.GetWidth())) << "%     ";
					}
				}
			}));
		}

		auto allReady = false;

		while(!allReady)
		{
			allReady = true;

			for(auto& task : tasks)
			{
				if(task.wait_for(std::chrono::milliseconds(30)) != std::future_status::ready)
				{
					allReady = false;
				}
			}

			liveImage.Update();
		}

		std::cout << std::endl;
	}

	const auto& GetImage() const { return _image; }
	const auto& GetCamera() const { return _camera; }
	
private:
	Image& _image;
	Camera& _camera;
	RaycasterConfiguration _config;
	RadianceProviderType& _radianceProvider;

	auto GetRandom() const
	{
		static std::default_random_engine _rnd;
		static std::uniform_real_distribution<float> _rng;
		return _rng(_rnd);
	}

	Color Sample(Image::SizeType x, Image::SizeType y, int, int, const Vector& viewPlaneBottomLeft, const Vector& xIncVector, const Vector& yIncVector, const Vector& apeture) const
	{
		const auto r1 = 2.0f * GetRandom();
		const auto r2 = 2.0f * GetRandom();

		/* Transform uniform into non-uniform filter samples */
		auto sdx = (r1 < 1.0f) ? sqrt(r1) - 1.0f : 1.0f - sqrt(2.0f - r1);
		auto sdy = (r2 < 1.0f) ? sqrt(r2) - 1.0f : 1.0f - sqrt(2.0f - r2);

		auto radiance = Color();
		auto dofSamples = _config.dofSamples;

		if(dofSamples > 0)
		{
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
			
			radiance /= static_cast<float>(dofSamples);
		}
		else
		{
			const auto viewPlanePoint = viewPlaneBottomLeft
				+ static_cast<float>(sdx + x) * xIncVector
				+ static_cast<float>(sdy + y) * yIncVector;

			const auto eyePoint = _camera.GetPosition();
			const auto dir = glm::normalize(viewPlanePoint - eyePoint);

			radiance = GetRadiance(Ray(eyePoint, dir));
		}

		return radiance;
	}

	Color GetRadiance(const Ray& ray) const
	{
		return _radianceProvider.GetRadiance(ray);
	}
};

#endif // Raytracer_H
