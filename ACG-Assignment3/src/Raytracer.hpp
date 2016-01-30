#ifndef Raytracer_H
#define Raytracer_H

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

	#define drand48() (((double)rand())/((double)RAND_MAX))
#endif


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

		const auto numThreads = std::thread::hardware_concurrency();
		const auto numXperThread = static_cast<Image::SizeType>(_image.GetWidth() / numThreads + (_image.GetWidth() % numThreads == 0 ? 0 : 1));
		//const auto numYperThread = static_cast<Image::SizeType>(_image.GetHeight() / numThreads + (_image.GetHeight() % numThreads == 0 ? 0 : 1));

		LiveImage liveImage(_image.GetWidth(), _image.GetHeight());
		liveImage.show();

		std::atomic<Image::SizeType> fragmentsDone = 0;
		std::vector<std::future<void>> tasks;

		for (auto tidx = 0ul; tidx < numThreads; tidx++)
		{
			const auto startX = tidx * numXperThread;
			const auto startY = 0;

			auto endX = startX + numXperThread;
			if (endX >= _image.GetWidth()) endX = _image.GetWidth();

			auto endY = _image.GetHeight();

			tasks.push_back(std::async(std::launch::deferred, [&, this, startY, startX, endY, endX, tidx]()
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
									futures.push_back(std::async(std::launch::deferred, [&, this]()
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

						liveImage.set(x, y, accumulated_radiance);
					}
					
					fragmentsDone += endX - startX;
					if (tidx == 0)
					{
						std::cout << "\rRendering (" << _config.subSamplesPerPixel * _config.samplesPerSubSample * _config.samplesPerSubSample << " spp) " << (100.0 * fragmentsDone / (_image.GetHeight() * _image.GetWidth())) << "%     ";
					}

					liveImage.update();
				}
			}));
		}

		for(auto& task : tasks)
		{
			task.wait();
		}

		/* Loop over image rows
		for (auto y = 0u; y < _image.GetHeight(); y++)
		{
			std::cout << "\rRendering (" << _config.subSamplesPerPixel * _config.samplesPerSubSample * _config.samplesPerSubSample << " spp) " << (100.0 * y / (_image.GetHeight() - 1)) << "%     ";

			for (auto x = 0u; x < _image.GetWidth(); x++)
			{
				std::vector<std::future<Color>> futures;
				
				for (auto sy = 0u; sy < _config.subSamplesPerPixel; sy++)
				{
					for (auto sx = 0u; sx < _config.subSamplesPerPixel; sx++)
					{
						for (auto s = 0u; s < _config.samplesPerSubSample; s++)
						{	
							futures.push_back(std::async([&, this] ()
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
				for(auto& future : futures)
				{
					accumulated_radiance += glm::clamp(future.get(), 0.0f, 1.0f)
						/ static_cast<float>(_config.subSamplesPerPixel * _config.subSamplesPerPixel);
				}
				
				_image.AddColor(x, y, accumulated_radiance);
			}
		} */

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

	std::uniform_real_distribution<float> _rng;

	auto GetRandom() const { static std::default_random_engine _rnd; return _rng(_rnd); }

	Color Sample(Image::SizeType x, Image::SizeType y, int, int, const Vector& viewPlaneBottomLeft, const Vector& xIncVector, const Vector& yIncVector, const Vector& apeture) const
	{
		const auto r1 = 2.0f * GetRandom();
		const auto r2 = 2.0f * GetRandom();

		/* Transform uniform into non-uniform filter samples */
		auto sdx = (r1 < 1.0f) ? sqrt(r1) - 1.0f : 1.0f - sqrt(2.0f - r1);
		auto sdy = (r2 < 1.0f) ? sqrt(r2) - 1.0f : 1.0f - sqrt(2.0f - r2);

		auto radiance = Color();

		auto dofSamples = _config.dofSamples;
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

	Color GetRadiance(const Ray& ray, int depth = 0) const
	{
		IntersectionInfo intersectionInfo;
		TracingInfo tracingInfo;

		if (_scene.Intersect(ray, intersectionInfo))
		{
			auto radience = Radiance(ray, 0, 1);
			return radience;
		}


		return _clearColor;
	}


	


	Color Radiance(const Ray &ray, int depth, float E) const
	{
		IntersectionInfo intersectionInfo;
		TracingInfo tracingInfo;

		depth++;

		if (!_scene.Intersect(ray, intersectionInfo))   /* No intersection with scene */
			return _clearColor;

		//const auto &obj = *(mIntersectionInfo.geometry.get());
		const auto &obj = *intersectionInfo.geometry;


		auto glossiness = obj.GetMaterial().GetGlossiness();
		auto emission = obj.GetMaterial().GetEmission() * E;
		auto col = obj.GetMaterial().GetColor();

		/* Maximum RGB reflectivity for Russian Roulette */
		float p = glm::max(col.x,glm::max(col.y, col.z));

		if (depth > 5 || !p)   /* After 5 bounces or if max reflectivity is zero */
		{
			if (drand48() < p*0.9)            /* Russian Roulette */
				col = col * (1.0f / p);        /* Scale estimator to remain unbiased */
			else
				return emission;  /* No further bounces, only return potential emission */
		}


		if (obj.GetMaterial().GetReflectionType() == DIFF)
		{
			return emission + _radianceProvider.GetRadiance(intersectionInfo, tracingInfo);
		}
		else if(obj.GetMaterial().GetReflectionType() == SPEC)
		{
			auto reflt = glm::reflect(ray.GetDirection(), intersectionInfo.normal);
			varyVector(reflt, glossiness);

			return emission + col * Radiance(Ray(intersectionInfo.hitpoint, reflt), depth, 1);
		}
		else if (obj.GetMaterial().GetReflectionType() == TRAN)
		{

			auto hitpoint = intersectionInfo.hitpoint;
			auto normal = intersectionInfo.normal;
			auto origDir = ray.GetDirection();


			Ray perfectReflective(hitpoint, glm::reflect(origDir, normal));





			Vector nl = normal;

			/* Obtain flipped normal, if object hit from inside */
			if (glm::dot(normal, origDir) >= 0.f)
				nl = nl*-1.0f;

			/* Otherwise object transparent, i.e. assumed dielectric glass material */

			bool into = glm::dot(normal, nl) > 0;       /* Bool for checking if ray from outside going in */
			float nc = 1.f;                        /* Index of refraction of air (approximately) */
			float nt = 1.5f;                      /* Index of refraction of glass (approximately) */
			float nnt;

			if (into)      /* Set ratio depending on hit from inside or outside */
				nnt = nc / nt;
			else
				nnt = nt / nc;

			float ddn = glm::dot(origDir, nl);
			float cos2t = 1.f - nnt * nnt * (1.f - ddn*ddn);

			/* Check for total internal reflection, if so only reflect */
			if (cos2t < 0)
				return emission + col * (Radiance(perfectReflective, depth, 1));

			/* Otherwise reflection and/or refraction occurs */
			Vector tdir;

			/* Determine transmitted ray direction for refraction */
			if (into)
				tdir = (origDir * nnt - normal * (ddn * nnt + sqrt(cos2t)));
			else
				tdir = (origDir * nnt + normal * (ddn * nnt + sqrt(cos2t)));

			tdir = glm::normalize(tdir);

			/* Determine R0 for Schlick´s approximation */
			float a = nt - nc;
			float b = nt + nc;
			float R0 = a*a / (b*b);

			/* Cosine of correct angle depending on outside/inside */
			float c;
			if (into)
				c = 1 + ddn;
			else
				c = 1 - glm::dot(tdir, normal);

			/* Compute Schlick´s approximation of Fresnel equation */
			float Re = R0 + (1.f - R0) *c*c*c*c*c;   /* Reflectance */
			float Tr = 1.f - Re;                     /* Transmittance */

													 /* Probability for selecting reflectance or transmittance */
			float P = .25f + .5f * Re;
			float RP = Re / P;         /* Scaling factors for unbiased estimator */
			float TP = Tr / (1.f - P);

			if (depth < 3)   /* Initially both reflection and trasmission */
				return emission + col*(Radiance(perfectReflective, depth, 1) * Re +
					Radiance(Ray(hitpoint, tdir), depth, 1) * Tr);
			else             /* Russian Roulette */
				if (drand48() < P)
					return emission + col*(Radiance(perfectReflective, depth, 1) * RP);
				else
					return emission + col*(Radiance(Ray(hitpoint, tdir), depth, 1) * TP);
		}
		
		assert("Unknown material");
		return Color();
	}


	// myFunction2
	// input:  x should be a random number between 0 and 1
	//         glossiness should be a value between 0 and 1
	//         small glossiness means really blurry reflection (evenly distributed randomness)
	//         large glossiness means really sharp reflection (unevenly distributed randomness)
	// output: a value between 0 and 1
	float myFunction2(float x, float glossiness) const
	{
		if (glossiness < 0) 
			glossiness = 0;

		if (glossiness >= 1.0) 
			glossiness = 1.0 - 1.0e-6;

		float c = -1.0f / glm::log(glossiness);

		return glm::exp(c * x - c);
	}

	void varyVector(Vector &vector, float glossiness) const
	{
		if (glossiness == 0.0f)
			return;

		/* Compute random reflection vector on half hemisphere */
		float oldtheta = glm::acos(vector.z);
		float oldphi = glm::atan(vector.y / vector.x);
		float deltatheta = PI * (2.0f * drand48() - 1);
		float deltaphi = (PI / 4.0f) * (2.0f * drand48() - 1.0f);

		float theta = oldtheta + deltatheta;
		float phi = oldphi + deltaphi;

		/* Random reflection vector d */
		auto d = Vector(glm::sin(theta) * glm::cos(phi),
			glm::sin(theta) * glm::sin(phi),
			glm::cos(theta));

		vector += myFunction2(drand48(), glossiness) * d; // myFunction2 works better than myFunction1, since the vector will still point in the same general direction
		vector = glm::normalize(vector);
	}
};

#endif // Raytracer_H
