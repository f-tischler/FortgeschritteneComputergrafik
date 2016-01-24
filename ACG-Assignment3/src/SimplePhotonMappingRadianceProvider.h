#ifndef SimplePhotonMappingRadianceProvider_H
#define SimplePhotonMappingRadianceProvider_H

#include <map>
#include <random>
#include "glm/gtx/norm.hpp"

#include "Scene.hpp"
#include "Sphere.hpp"
#include <list>
#include "kdtree.h"


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
					photon.radiance = material.GetEmission();
					photon.direction = glm::normalize(Vector(GetRandom(), GetRandom(), GetRandom()));
					photon.position = sphere->GetPosition() + photon.direction * (sphere->GetRadius() + 0.00001f);
					
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

		if (_photonMap.find(reinterpret_cast<size_t>(intersectionInfo.geometry)) == _photonMap.end()) return Color(0, 0, 0);

		auto color = Color(0, 0, 0);

        struct kdres* set = kd_nearest_range3f(intersectionInfo.geometry->GetTree(), intersectionInfo.hitpoint.x, intersectionInfo.hitpoint.y, intersectionInfo.hitpoint.z, 100.0f);
        
        std::vector<Photon> photons;
        while(!kd_res_end(set)) {
            Photon* ptr = static_cast<Photon*>(kd_res_item_data(set));
            photons.push_back(*ptr);
            kd_res_next(set);
        }

		/*std::sort(photons.begin(), photons.end(),
			[&intersectionInfo] (const Photon& p1, const Photon& p2)
		{
			auto distP1 = glm::length2(intersectionInfo.hitpoint - p1.position);
			auto distP2 = glm::length2(intersectionInfo.hitpoint - p2.position);

			return distP1 < distP2;
		});*/

		auto n = 0;
		
		constexpr auto maxPhotonGathered = 50;
        
		for (auto& photon : photons)
		{
			auto weight = std::max(0.0f, -glm::dot(intersectionInfo.normal, photon.direction));
			//weight *= (1.0f - glm::sqrt(distance));/// exposure;
			color += photon.radiance * weight;

			if (++n > maxPhotonGathered) break;
		}

		auto maxDistSq = 0.0f;
        if(n == 0) {
            return Color(0,0,0);
        } else if(n < photons.size())
		{
			maxDistSq = glm::length2(intersectionInfo.hitpoint - photons[n].position);
		}
		else
		{
			maxDistSq = glm::length2(intersectionInfo.hitpoint - photons.back().position);
		}
        
        kd_res_free(set);
        
		return color / (maxDistSq * PI);
	}

	Color DisplayPhoton(const IntersectionInfo& intersectionInfo, TracingInfo& tracingInfo) const
	{
		if (_photonMap.find(reinterpret_cast<size_t>(intersectionInfo.geometry)) == _photonMap.end()) return Color(0, 0, 0);
		
        struct kdres* set = kd_nearest_range3f(intersectionInfo.geometry->GetTree(), intersectionInfo.hitpoint.x, intersectionInfo.hitpoint.y, intersectionInfo.hitpoint.z, 100.0f);
        
        std::vector<Photon> photons;
        while(!kd_res_end(set)) {
            Photon* ptr = static_cast<Photon*>(kd_res_item_data(set));
            photons.push_back(*ptr);
            kd_res_next(set);
        }
        
		for (auto& photon : photons)
		{
			if(glm::length2(photon.position - intersectionInfo.hitpoint) < debugEpsilon)
			{
				return Color(1, 1, 1);
			}
		}

		return Color(0, 0, 0);
	}

private:
	static constexpr auto maxDepth = 4;
	static constexpr auto photonCount = 1000;
	static constexpr auto debugEpsilon = 0.1f;

	struct Photon
	{
		Color radiance;
		Vector position;
		Vector direction;
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
			photon.radiance *= info.geometry->GetMaterial().GetColor();
			_photonMap[reinterpret_cast<size_t>(info.geometry)].push_back(photon);
            kd_insert3f(info.geometry->GetTree(), photon.position.x, photon.position.y, photon.position.z, &photon);

			photon.direction = glm::reflect(photon.direction, info.normal);
			photon.radiance *= 1.0f / glm::sqrt(depth);

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
