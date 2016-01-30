#ifndef NanoFlannKdPMRP_H
#define NanoFlannKdPMRP_H

#include <map>
#include <random>
#include "glm/gtx/norm.hpp"

#include "Scene.hpp"
#include "Sphere.hpp"
#include "nanoflann.hpp"
#include <list>

struct TracingInfo;

struct Photon
{
	Color radiance;
	Vector position;
	Vector direction;
};

struct NeighbourSet
{
	std::vector<size_t> indices;
	std::vector<float> distances;
	float worstDistance;
};

class Photons
{
public:
	Photons() : _data(0), _index(3, *this, nanoflann::KDTreeSingleIndexAdaptorParams(10)) { }

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
		std::lock_guard<std::mutex> lock(_dataMtx);
		_data.push_back(photon);
	}

	auto GetNeighbours(const Vector& pos, const size_t n) const
	{
		const float query_pt[3] = { pos.x, pos.y, pos.z };
		const auto num_results = n;

		NeighbourSet set;
		set.indices = std::vector<size_t>(n, -1);
		set.distances = std::vector<float>(n, -1);

		nanoflann::KNNResultSet<float> resultSet(num_results);

		resultSet.init(set.indices.data(), set.distances.data());
		_index.findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));

		set.worstDistance = resultSet.worstDist();

		return set;
	}

	const auto& GetPhotons() const { return _data; }

	void BuildIndex()
	{
		if(!_data.empty())
			_index.buildIndex();
	}

	bool Empty() const { return _data.empty(); }

private:
	using MyKdTreeType = 
		nanoflann::KDTreeSingleIndexAdaptor<
		nanoflann::L2_Simple_Adaptor<float, Photons>,
		Photons,
		3>;

	std::vector<Photon> _data;
	MyKdTreeType _index;
	std::mutex _dataMtx;
};

class NanoFlannKdPMRP
{
public:
	explicit NanoFlannKdPMRP(bool debug = false)
		: _debug(debug) { }

	void CreatePhotonMap(const Scene& scene)
	{
		const auto& objects = scene.GetGeometry();

		for (auto& object : objects)
			_photonMap.try_emplace(reinterpret_cast<size_t>(object.get()));

		std::list<std::future<void>> tasks;
		auto threads = std::thread::hardware_concurrency();
		
		for (auto& object : objects)
		{
			const auto& material = object->GetMaterial();

			const auto& emission = material.GetEmission();

			if (object->GetType() == eGeometryType_Sphere &&
				material.GetEmission() != Vector(0, 0, 0))
			{
				for (auto t = 0ul; t < threads; t++)
				{						
					tasks.push_back(std::async([&]
					{
						for (auto i = 0ul; i < photonCount / threads; i++)
						{
							const auto sphere = static_cast<Sphere*>(object.get());

							Photon photon;
							photon.radiance = material.GetEmission() / static_cast<float>(photonCount);
							photon.direction = glm::normalize(Vector(GetRandom(), GetRandom(), GetRandom()));
							photon.position = sphere->GetPosition() + photon.direction * (sphere->GetRadius() + 0.00001f);

							SendPhoton(scene, photon);
						}
					}));
				}

				for (auto& t : tasks)
					t.wait();

				tasks.clear();
			}
		}

		for(auto& p : _photonMap)
		{
			tasks.push_back(std::async([&]()
			{
				p.second.BuildIndex();
			}));

			if(tasks.size() >= threads*2)
			{
				for (auto& t : tasks)
					t.wait();

				tasks.clear();
			}
		}

		for(auto& task : tasks)
		{
			task.wait();
		}
	}

	Color GetRadiance(const IntersectionInfo& intersectionInfo, TracingInfo& tracingInfo) const
	{
		if (_debug) return DisplayPhoton(intersectionInfo, tracingInfo);

		if (intersectionInfo.geometry->GetMaterial().GetEmission() != Vector(0, 0, 0))
			return intersectionInfo.geometry->GetMaterial().GetEmission();

		auto it = _photonMap.find(reinterpret_cast<size_t>(intersectionInfo.geometry));
		if (it == _photonMap.end()) return Color(0, 0, 0);

		if (it->second.Empty()) return Color(0, 0, 0);

		auto color = Color(0, 0, 0);
		
		auto neighbours = it->second.GetNeighbours(intersectionInfo.hitpoint, 200);

		float worstDistance = glm::sqrt(neighbours.worstDistance);
		float k = 1;

		auto normalizationCoeff = 1 - 2 / (3 * k);

		for (auto& idx : neighbours.indices)
		{
			if (idx < 0) continue;

			auto& photon = it->second.GetPhotons()[idx];
			auto weight = std::max(0.0f, -glm::dot(intersectionInfo.normal, photon.direction));
			weight *= 1 - glm::length(intersectionInfo.hitpoint - photon.position) / (worstDistance * k);/// exposure;
			color += photon.radiance * weight;
		}

		return color / (neighbours.worstDistance * PI * normalizationCoeff) * intersectionInfo.geometry->GetMaterial().GetColor();
	}

	Color DisplayPhoton(const IntersectionInfo& intersectionInfo, TracingInfo& tracingInfo) const
	{
		if (intersectionInfo.geometry->GetMaterial().GetEmission() != Vector(0, 0, 0))
			return intersectionInfo.geometry->GetMaterial().GetEmission();

		auto it = _photonMap.find(reinterpret_cast<size_t>(intersectionInfo.geometry));
		if (it == _photonMap.end()) return Color(0, 0, 0);

		if (it->second.Empty()) return Color(0, 0, 0);

		auto color = Color(0, 0, 0);

		auto neighbour = it->second.GetNeighbours(intersectionInfo.hitpoint, 1);
		if (neighbour.worstDistance < debugEpsilon && neighbour.indices[0] >= 0)
			return it->second.GetPhotons()[neighbour.indices[0]].radiance;

		return Color(0, 0, 0);
	}

private:
	static constexpr auto maxDepth = 5;
	static constexpr auto photonCount = 5000000;
	static constexpr auto debugEpsilon = 0.05f;

	float GetRandom() { return _rng(_rnd); }

	void SendPhoton(const Scene& scene, Photon& photon, int depth = 1)
	{
		Ray ray(photon.position, photon.direction);
		IntersectionInfo info;

		if (scene.Intersect(ray, info))
		{
			const auto& obj = *info.geometry;

			auto glossiness = obj.GetMaterial().GetGlossiness();
			auto emission = obj.GetMaterial().GetEmission();
			auto col = obj.GetMaterial().GetColor();

			/* Maximum RGB reflectivity for Russian Roulette */
			float p = glm::max(col.x, glm::max(col.y, col.z));

			if (depth > maxDepth || !p)   /* After 5 bounces or if max reflectivity is zero */
			{
				if (drand48() < p*0.9)            /* Russian Roulette */
					col = col * (1.f / p);        /* Scale estimator to remain unbiased */
				else
					return;  /* No further bounces, only return potential emission */
			}

			photon.position = info.hitpoint;
			photon.radiance *= 0.8f;

			if(depth >= 1)
			{
				auto& photons = _photonMap[reinterpret_cast<size_t>(info.geometry)];
				photons.Add(photon);
			}

			switch (obj.GetMaterial().GetReflectionType())
			{
				case DIFF:
				{
					photon.radiance *= info.geometry->GetMaterial().GetColor();
					photon.direction = glm::reflect(photon.direction, info.normal);
					photon.radiance *= 0.8f;

					SendPhoton(scene, photon, ++depth);
					break;
				}
				case SPEC:
				{
					photon.radiance *= col;
					photon.direction = glm::reflect(photon.direction, info.normal);

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
						photon.radiance *= col;
						photon.direction = glm::reflect(photon.direction, info.normal);

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

					if (depth < maxDepth)   /* Initially both reflection and trasmission */
					{
						photon.radiance *= col;
						photon.direction = perfectReflective.GetDirection();

						auto newPhoton = photon;
						SendPhoton(scene, newPhoton, ++depth);

						photon.direction = tdir;
						auto newPhoton2 = photon;
						SendPhoton(scene, newPhoton2, ++depth);
					}
					else             /* Russian Roulette */
					{
						if (drand48() < P)
						{
							photon.radiance *= col;
							photon.direction = perfectReflective.GetDirection();

							SendPhoton(scene, photon, ++depth);

							photon.radiance *= RP;
						}
						else
						{
							photon.radiance *= col;
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


	using PhotonMap = std::map<size_t, Photons>;
	PhotonMap _photonMap;

	std::vector<Photon> _data;

	std::default_random_engine _rnd;
	std::uniform_real_distribution<float> _rng = std::uniform_real_distribution<float>(-1.0f, 1.0f);
	bool _debug;
};

#endif // NanoFlannKdPMRP_H

