#ifndef SimplePhotonMappingRadianceProvider_H
#define SimplePhotonMappingRadianceProvider_H

#include <map>
#include <random>
#include <cmath>
#include "glm/gtx/norm.hpp"

#include "Scene.hpp"
#include "Sphere.hpp"
#include <list>


struct TracingInfo;

class SimplePhotonMappingRadianceProvider
{
public:
	SimplePhotonMappingRadianceProvider(bool debug = false)
		: _debug(debug) { }

	void CreatePhotonMap(const Scene& scene)
	{
		const auto& objects = scene.GetGeometry();

		for(auto& object : objects)
		{
			const auto& material = object->GetMaterial();

			const auto& emission = material.GetEmission();
			
			if(object->GetType() == eGeometryType_Sphere && 
			   material.GetEmission() != Vector(0, 0, 0))
			{
				for (auto i = 0ul; i < photonCount; i++)
				{
					const auto sphere = static_cast<Sphere*>(object.get());

					Photon photon;
					photon.radiance = material.GetEmission() / static_cast<float>(photonCount);
					photon.direction = glm::normalize(Vector(GetRandom(), GetRandom(), GetRandom()));
					photon.position = sphere->GetPosition() + photon.direction * (sphere->GetRadius() + 0.00001f);
					photon.depth = 0;

					SendPhoton(scene, photon);
				}	
			}
		}
	}

	Color GetRadiance(const IntersectionInfo& intersectionInfo, TracingInfo& tracingInfo) const
	{
		if (_debug) return DisplayPhoton(intersectionInfo, tracingInfo);

		//Is light
		if (intersectionInfo.geometry->GetMaterial().HasEmission())
			return intersectionInfo.geometry->GetMaterial().GetEmission();







		//Has photon

		auto it = _photonMap.find(reinterpret_cast<size_t>(intersectionInfo.geometry));
		if (it == _photonMap.end()) 
			return Color(0, 0, 0);



		auto photons = it->second;

		for(auto& photon : photons)
		{
			photon.distance =  glm::length(intersectionInfo.hitpoint - photon.position);
		}

		std::sort(photons.begin(), photons.end(), 
			[&intersectionInfo] (const Photon& p1, const Photon& p2)
		{
			return p1.distance < p2.distance;
		});
		
		constexpr auto maxPhotonGathered = 100;
		auto maxDist = 0.0f;
		if (maxPhotonGathered < photons.size())
		{
			maxDist = photons[maxPhotonGathered].distance;
		}
		else
		{
			maxDist = photons.back().distance;
		}

		auto n = 0;	
		constexpr auto k = 1;
		auto color = Color(0, 0, 0);

		for (auto& photon : photons)
		{
			auto weight = std::max(0.0f, -glm::dot(intersectionInfo.normal, photon.direction));

			weight *= 1.0f - photon.distance / (k * maxDist);/// exposure;
			color += photon.radiance * weight;

			if (++n > maxPhotonGathered) break;
		}

		return (color / (maxDist * maxDist * PI)) * intersectionInfo.geometry->GetMaterial().GetColor();
	}

	Color DisplayPhoton(const IntersectionInfo& intersectionInfo, TracingInfo& tracingInfo) const
	{
		if (_photonMap.find(reinterpret_cast<size_t>(intersectionInfo.geometry)) == _photonMap.end()) return Color(0, 0, 0);
		
		auto photons = _photonMap.find(reinterpret_cast<size_t>(intersectionInfo.geometry))->second;
		for (auto& photon : photons)
		{
			if(glm::length2(photon.position - intersectionInfo.hitpoint) < debugEpsilon)
			{
				return photon.radiance;
			}
		}

		return Color(0, 0, 0);
	}

private:
	static constexpr auto maxDepth = 3;
	static constexpr auto photonCount = 5000;
	static constexpr auto debugEpsilon = 0.01f;

	struct Photon
	{
		Color radiance;
		Vector position;
		Vector direction;
		float distance;
		int depth;
	};

	float GetRandom() { return _rng(_rnd); }

	void SendPhoton(const Scene& scene, Photon& photon, int depth = 1)
	{
		if (depth >= maxDepth) 
			return;

		Ray ray(photon.position, photon.direction);
		IntersectionInfo info;

		if(scene.Intersect(ray, info))
		{
			const auto& obj = *info.geometry;

			auto glossiness = obj.GetMaterial().GetGlossiness();
			auto emission = obj.GetMaterial().GetEmission();
			auto col = obj.GetMaterial().GetColor();

			/* Maximum RGB reflectivity for Russian Roulette */
			float p = glm::max(col.x, glm::max(col.y, col.z));

			if (depth > 5 || !p)   /* After 5 bounces or if max reflectivity is zero */
			{
				if (drand48() < p*0.9)            /* Russian Roulette */
					col = col * (1.f / p);        /* Scale estimator to remain unbiased */
				else
					return;  /* No further bounces, only return potential emission */
			}

			switch(obj.GetMaterial().GetReflectionType())
			{
				case DIFF:
				{
					photon.position = info.hitpoint;
					_photonMap[reinterpret_cast<size_t>(info.geometry)].push_back(photon);

					photon.radiance *= info.geometry->GetMaterial().GetColor();
					photon.direction = glm::reflect(photon.direction, info.normal);
					photon.radiance *= 0.8f;
					photon.depth++;

					SendPhoton(scene, photon, ++depth);
					break;
				}
				case SPEC:
				{
					photon.position = info.hitpoint;
					_photonMap[reinterpret_cast<size_t>(info.geometry)].push_back(photon);


					photon.radiance *= col;
					photon.direction = glm::reflect(photon.direction, info.normal);
					photon.depth++;

					varyVector(photon.direction, glossiness);

					SendPhoton(scene, photon, ++depth);
					break;
				}
				case TRAN:
				{
					auto hitpoint = info.hitpoint;
					auto normal = info.normal;
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
					{
						photon.position = info.hitpoint;
						_photonMap[reinterpret_cast<size_t>(info.geometry)].push_back(photon);


						photon.radiance *= col;
						photon.direction = glm::reflect(photon.direction, info.normal);
						photon.depth++;

						varyVector(photon.direction, glossiness);

						SendPhoton(scene, photon, ++depth);
						break;
					}

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

					//leave photon
					photon.position = info.hitpoint;
					_photonMap[reinterpret_cast<size_t>(info.geometry)].push_back(photon);

					if (depth < 4)   /* Initially both reflection and trasmission */
					{
						photon.radiance *= col;
						photon.depth++;

						photon.direction = perfectReflective.GetDirection();
						SendPhoton(scene, photon, ++depth);

						photon.direction = tdir;
						SendPhoton(scene, photon, ++depth);
					}
					else             /* Russian Roulette */
					{
						if (drand48() < P)
						{
							photon.radiance *= col;
							photon.depth++;

							photon.direction = perfectReflective.GetDirection();
							SendPhoton(scene, photon, ++depth);
							photon.radiance *= RP;
						}
						else
						{
							photon.radiance *= col;
							photon.depth++;

							photon.direction = tdir;
							SendPhoton(scene, photon, ++depth);
							photon.radiance *= TP;
						}
					}

					break;
				}
			}
		}
	}

	// myFunction2
	// input:  x should be a random number between 0 and 1
	//         glossiness should be a value between 0 and 1
	//         small glossiness means really blurry reflection (evenly distributed randomness)
	//         large glossiness means really sharp reflection (unevenly distributed randomness)
	// output: a value between 0 and 1
	float myFunction2(float x, float glossiness) const
	{
		if (glossiness < 0.f)
			glossiness = 0.f;

		if (glossiness >= 1.0f)
			glossiness = 1.0f - 1.0e-6f;

		auto c = -1.0f / glm::log(glossiness);

		return glm::exp(c * x - c);
	}

	void varyVector(Vector &vector, float glossiness) const
	{
		if (glossiness == 0.0f)
			return;

		/* Compute random reflection vector on half hemisphere */
		float oldtheta = glm::acos(vector.z);
		float oldphi = glm::atan(vector.y / vector.x);
		float deltatheta = PI * (2.0f * drand48() - 1.f);
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

	using PhotonMap = std::map<size_t, std::vector<Photon>>;
	PhotonMap _photonMap;

	std::default_random_engine _rnd;
	std::uniform_real_distribution<float> _rng = std::uniform_real_distribution<float>(-1.0f, 1.0f);
	bool _debug;
};

#endif // SimplePhotonMappingRadianceProvider_H
