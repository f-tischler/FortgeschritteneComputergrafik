#ifndef SimplePhotonMappingRadianceProvider_H
#define SimplePhotonMappingRadianceProvider_H

#include <map>
#include <random>
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

		if (intersectionInfo.geometry->GetMaterial().GetEmission() != Vector(0, 0, 0))
			return intersectionInfo.geometry->GetMaterial().GetEmission();

		auto it = _photonMap.find(reinterpret_cast<size_t>(intersectionInfo.geometry));
		if (it == _photonMap.end()) return Color(0, 0, 0);

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
	static constexpr auto photonCount = 100000;
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
		if (depth >= maxDepth) return;

		Ray ray(photon.position, photon.direction);
		IntersectionInfo info;

		if(scene.Intersect(ray, info))
		{
			photon.position = info.hitpoint;
			_photonMap[reinterpret_cast<size_t>(info.geometry)].push_back(photon);

			photon.radiance *= info.geometry->GetMaterial().GetColor();
			photon.direction = glm::reflect(photon.direction, info.normal);
			photon.radiance *= 0.8f;
			photon.depth++;

			SendPhoton(scene, photon, ++depth);
		}
	}

	using PhotonMap = std::map<size_t, std::vector<Photon>>;
	PhotonMap _photonMap;

	std::default_random_engine _rnd;
	std::uniform_real_distribution<float> _rng = std::uniform_real_distribution<float>(-1.0f, 1.0f);
	bool _debug;
};

#endif // SimplePhotonMappingRadianceProvider_H
