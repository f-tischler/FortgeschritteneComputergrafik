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
			auto photonRadiance = 0.0f;// _radianceProvider.GetRadiance(intersectionInfo, tracingInfo, depth);
			auto radience = Radiance(ray, 0, 1);
			return radience + photonRadiance;
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


		auto col = obj.GetMaterial().GetColor();

		/* Maximum RGB reflectivity for Russian Roulette */
		float p = glm::max(col.x,glm::max(col.y, col.z));

		if(depth > 1)


		if (depth > 5 || !p)   /* After 5 bounces or if max reflectivity is zero */
		{
			if (drand48() < p*0.9)            /* Russian Roulette */
				col = col * (1.0f / p);        /* Scale estimator to remain unbiased */
			else
				return obj.GetMaterial().GetEmission() * E;  /* No further bounces, only return potential emission */
		}


		if (obj.GetMaterial().GetReflectionType() == DIFF)
		{
			/* Obtain flipped normal, if object hit from inside */
			auto nl = intersectionInfo.normal;
			if (glm::dot(intersectionInfo.normal, ray.GetDirection()) >= 0)
				nl = nl*-1.0f;

			/* Compute random reflection vector on hemisphere */
			float r1  = 2.0 * PI * drand48();
			float r2  = drand48();
			float r2s = sqrt(r2);

			/* Set up local orthogonal coordinate system u,v,w on surface */
			Vector w = nl;
			Vector u;

			if (fabs(w.x) > .1)
				u = Vector(0.0, 1.0, 0.0);
			else
				u = glm::normalize(glm::cross(Vector(1.0, 0.0, 0.0), w));

			Vector v = glm::cross(w,u);

			/* Random reflection vector d */
			Vector d =  u * glm::cos(r1) * r2s +
						v * glm::sin(r1) * r2s +
					    w * glm::sqrt(1.0f - r2);

			d = glm::normalize(d);

			/* Explicit computation of direct lighting */
			Vector e;
			for (auto i = 0u; i < _scene.GetGeometry().size(); i++)
			{
				const auto& geom = *(_scene.GetGeometry()[i].get());
				if (!geom.GetMaterial().HasEmission())
					continue;

				auto hitpoint = intersectionInfo.hitpoint;
				auto r = geom.GetBoundingBox().diameter() / 2;
				auto pos = geom.GetBoundingBox().position();

				/* Set up local orthogonal coordinate system su,sv,sw towards light source */
				Vector sw = pos - hitpoint;
				Vector su;

				if (fabs(sw.x) > 0.1)
					su = Vector(0.0, 1.0, 0.0);
				else
					su = Vector(1.0, 0.0, 0.0);

				su = glm::cross(su, w);
				su = glm::normalize(su);
				auto sv = glm::cross(sw, su);

				/* Create random sample direction l towards spherical light source */
				auto hit_Light = hitpoint - pos;

				float cos_a_max = sqrt(1.0f - r * r / glm::dot(hit_Light, hit_Light));
				float eps1 = drand48();
				float eps2 = drand48();
				float cos_a = 1.0f - eps1 + eps1 * cos_a_max;
				float sin_a = sqrt(1.0f - cos_a * cos_a);
				float phi = 2.0f*PI * eps2;

				Vector l = su * glm::cos(phi) * sin_a +
					sv * glm::sin(phi) * sin_a +
					sw * cos_a;

				l = glm::normalize(l);

				/* Shoot shadow ray, check if intersection is with light source */
				IntersectionInfo info;
				if (_scene.Intersect(Ray(hitpoint, l), info) && 
					info.geometry == _scene.GetGeometry()[i].get())
				{
					auto omega = 2 * PI * (1 - cos_a_max);
					e = e + col * (geom.GetMaterial().GetEmission() * glm::dot(l, nl) * omega) * 1.0f / PI;
				}
			}

			/* Return potential light emission, direct lighting, and indirect lighting (via
			recursive call for Monte-Carlo integration */

			col *= Radiance(Ray(intersectionInfo.hitpoint, d), depth, 0);

			return obj.GetMaterial().GetEmission() * E + e + col;
		}
		else if(obj.GetMaterial().GetReflectionType() == eReflectionType::SPEC)
		{
			auto reflt = glm::reflect(ray.GetDirection(), intersectionInfo.normal);
			col *= Radiance(Ray(intersectionInfo.hitpoint, reflt), depth, 1);
			return col;
		}
		else if (obj.GetMaterial().GetReflectionType() == GLOS)
		{
			auto L = ray.GetDirection();
			auto N = intersectionInfo.normal;
			auto L_prime = L - N * 2.0f * (glm::dot(N,L));

			varyVector(L_prime, obj.GetMaterial().GetGlossiness());

			col *= Radiance(Ray(intersectionInfo.hitpoint, L_prime), depth, 1);

			return obj.GetMaterial().GetEmission() + col;
		}
		else if (obj.GetMaterial().GetReflectionType() == TRAN)
		{
			auto isGlossy = obj.GetMaterial().GetGlossiness() > 0.0f;


			//perfectly reflection
			auto prefectReflected = glm::reflect(ray.GetDirection(), intersectionInfo.normal);

			//glossy material
			//randomize it
			if (isGlossy)
				varyVector(prefectReflected, obj.GetMaterial().GetGlossiness());

			Ray reflRay(intersectionInfo.hitpoint, prefectReflected);

			/* Obtain flipped normal, if object hit from inside */
			auto nl = intersectionInfo.normal;
			if (glm::dot(intersectionInfo.normal, ray.GetDirection()) >= 0)
				nl = nl*-1.0f;

			auto into = glm::dot(intersectionInfo.normal, nl) > 0;       /* Bool for checking if ray from outside going in */
			auto nc = 1.f;                        /* Index of refraction of air (approximately) */
			auto nt = 1.5f;                      /* Index of refraction of glass (approximately) */
			auto nnt = 0.0f;

			if (into)      /* Set ratio depending on hit from inside or outside */
				nnt = nc / nt;
			else
				nnt = nt / nc;

			auto ddn = glm::dot(ray.GetDirection(), nl);
			auto cos2t = 1.f - nnt * nnt * (1.f - ddn*ddn);


			/* Check for total internal reflection, if so only reflect */
			if (cos2t < 0)
				return obj.GetMaterial().GetEmission() + col * Radiance(reflRay, depth, 1);

			/* Otherwise reflection and/or refraction occurs */
			Vector tdir;

			/* Determine transmitted ray direction for refraction */
			if (into)
				tdir = glm::normalize(ray.GetDirection() * nnt - intersectionInfo.normal * (ddn * nnt + sqrt(cos2t)));
			else
				tdir = glm::normalize(ray.GetDirection() * nnt + intersectionInfo.normal * (ddn * nnt + sqrt(cos2t)));


			//glossy material
			//randomize it
			if (isGlossy)
				varyVector(tdir, obj.GetMaterial().GetGlossiness());


			/* Determine R0 for Schlick´s approximation */
			auto a = nt - nc;
			auto b = nt + nc;
			auto R0 = a*a / (b*b);

			/* Cosine of correct angle depending on outside/inside */
			auto c = 0.0f;
			if (into)
				c = 1 + ddn;
			else
				c = 1 - glm::dot(tdir, intersectionInfo.normal);

			/* Compute Schlick´s approximation of Fresnel equation */
			auto Re = R0 + (1 - R0) *c*c*c*c*c;   /* Reflectance */
			auto Tr = 1 - Re;                     /* Transmittance */

													/* Probability for selecting reflectance or transmittance */
			auto P = .25f + .5f * Re;
			auto RP = Re / P;         /* Scaling factors for unbiased estimator */
			auto TP = Tr / (1 - P);

			if (depth < 5)   /* Initially both reflection and trasmission */
				return obj.GetMaterial().GetEmission() + col * 
				(Radiance(reflRay, depth, 1) * Re + Radiance(Ray(intersectionInfo.hitpoint, tdir), depth, 1) * Tr);
			else             /* Russian Roulette */
				if (drand48() < P)
					return obj.GetMaterial().GetEmission() + col * (Radiance(reflRay, depth, 1) * RP);
				else
					return obj.GetMaterial().GetEmission() + col * (Radiance(Ray(intersectionInfo.hitpoint, tdir), depth, 1) * TP);
		}

		return col;
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

	void varyVector(Vector &vector, double glossiness) const
	{
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
