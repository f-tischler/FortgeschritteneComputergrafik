#ifndef KdPhotonMappingRadianceProvider_H
#define KdPhotonMappingRadianceProvider_H

#include <map>
#include <random>
#include "glm/gtx/norm.hpp"

#include "Scene.hpp"
#include "Sphere.hpp"
#include "nanoflann.hpp"

struct TracingInfo;

struct Photon
{
	Color radiance;
	Vector position;
	Vector direction;
};

class Photons
{
public:
	Photons() : _index(3, *this, nanoflann::KDTreeSingleIndexAdaptorParams(10)) { }

	size_t kdtree_get_point_count() const { return _data.size(); }

	auto kdtree_distance(const float* p1, const size_t idx_p2, size_t) const
	{
		const auto d0 = p1[0] - _data[idx_p2].position.x;
		const auto d1 = p1[1] - _data[idx_p2].position.y;
		const auto d2 = p1[2] - _data[idx_p2].position.z;
		return d0*d0 + d1*d1 + d2*d2;
	}

	auto kdtree_get_pt(const size_t idx, int dim) const
	{
		if (dim == 0) return _data[idx].position.x;
		if (dim == 1) return _data[idx].position.y;

		return _data[idx].position.z;
	}

	template <class BBOX>
	bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }

	void Add(const Photon& photon)
	{
		_data.push_back(photon);
	}

	auto GetNeighbours(const Vector& pos, const size_t n) const
	{
		const float query_pt[3] = { pos.x, pos.y, pos.z };
		const auto num_results = n;
		size_t ret_index;
		float out_dist_sqr;

		nanoflann::KNNResultSet<float> resultSet(num_results);

		resultSet.init(&ret_index, &out_dist_sqr);
		_index.findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));

		return resultSet.;
	}

	const auto& GetPhotons() const { return _data; }

	void BuildIndex()
	{
		_index.buildIndex();
	}

private:
	using MyKdTreeType = 
		nanoflann::KDTreeSingleIndexAdaptor<
		nanoflann::L2_Simple_Adaptor<float, Photons>,
		Photons,
		3>;
	
	MyKdTreeType _index;

	std::vector<Photon> _data;
};

class KdPhotonMappingRadianceProvider
{
public:
	explicit KdPhotonMappingRadianceProvider(bool debug = false)
		: _debug(debug) { }

	void CreatePhotonMap(const Scene& scene)
	{
		const auto& objects = scene.GetGeometry();

		for (auto& object : objects)
		{
			const auto& material = object->GetMaterial();

			const auto& emission = material.GetEmission();

			if (object->GetType() == eGeometryType_Sphere &&
				material.GetEmission() != Vector(0, 0, 0))
			{
				for (auto i = 0ul; i < photonCount; i++)
				{
					const auto sphere = static_cast<Sphere*>(object.get());

					Photon photon;
					photon.radiance = material.GetEmission() / static_cast<float>(photonCount);
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

		auto it = _photonMap.find(reinterpret_cast<size_t>(intersectionInfo.geometry));
		if (it == _photonMap.end()) return Color(0, 0, 0);

		auto color = Color(0, 0, 0);

		constexpr auto radius = 100.0f;
		auto set = kd_nearest_range3f(it->second, intersectionInfo.hitpoint.x, intersectionInfo.hitpoint.y, intersectionInfo.hitpoint.z, radius);

		std::vector<Photon> photons;
		while (!kd_res_end(set)) {
			auto ptr = static_cast<Photon*>(kd_res_item_data(set));
			photons.push_back(*ptr);
			kd_res_next(set);
		}

		auto n = 0;

		constexpr auto maxPhotonGathered = 50;

		for (auto& photon : photons)
		{
			auto weight = std::max(0.0f, -glm::dot(intersectionInfo.normal, photon.direction));
			//weight *= (1.0f - glm::sqrt(distance));/// exposure;
			color += photon.radiance * weight;

			if (++n > maxPhotonGathered) break;
		}

		kd_res_free(set);

		auto maxDistSq = 0.0f;
		if (n == 0) {
			return Color(0, 0, 0);
		}
		else if (n < photons.size())
		{
			maxDistSq = glm::length2(intersectionInfo.hitpoint - photons[n].position);
		}
		else
		{
			maxDistSq = glm::length2(intersectionInfo.hitpoint - photons.back().position);
		}

		maxDistSq = radius * radius;



		return color / (maxDistSq * PI);
	}

	Color DisplayPhoton(const IntersectionInfo& intersectionInfo, TracingInfo& tracingInfo) const
	{
		auto it = _photonMap.find(reinterpret_cast<size_t>(intersectionInfo.geometry));
		if (it == _photonMap.end()) return Color(0, 0, 0);

		auto set = it->second.GetNeighbours(intersectionInfo.hitpoint, 50);

		std::vector<Photon> photons;
		for (auto i = 0ul; i < set.size; i++)
		{
			photons.push_back(it->second.GetPhotons()[set.indices[i]]);
		}

		for (auto& photon : photons)
		{
			if (glm::length2(photon.position - intersectionInfo.hitpoint) < debugEpsilon)
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

	float GetRandom() { return _rng(_rnd); }

	void SendPhoton(const Scene& scene, Photon& photon, int depth = 1)
	{
		if (depth >= maxDepth) return;

		Ray ray(photon.position, photon.direction);
		IntersectionInfo info;

		if (scene.Intersect(ray, info))
		{
			photon.position = info.hitpoint;

			_photonMap[reinterpret_cast<size_t>(info.geometry)].Add(photon);

			photon.direction = glm::reflect(photon.direction, info.normal);
			photon.radiance *= info.geometry->GetMaterial().GetColor() * 0.8f;

			SendPhoton(scene, photon, ++depth);
		}
	}

	using PhotonMap = std::map<size_t, Photons>;
	PhotonMap _photonMap;

	std::vector<Photon> _data;

	std::default_random_engine _rnd;
	std::uniform_real_distribution<float> _rng = std::uniform_real_distribution<float>(-1.0f, 1.0f);
	bool _debug;
};

#endif // KdPhotonMappingRadianceProvider_H

