#ifndef PathTracingWithPhotonMappingRP_H
#define PathTracingWithPhotonMappingRP_H

#include <map>
#include <random>
#include "glm/gtx/norm.hpp"

#include "Scene.hpp"
#include "Sphere.hpp"
#include "nanoflann.hpp"
#include <list>
#include "Raycaster.hpp"

template<class PhotonMapType>
class PathTracingWithPhotonMappingRP
{
public:
	explicit PathTracingWithPhotonMappingRP(const Scene& scene, bool debug = false, bool softShadows = false)
		: _debug(debug), _scene(scene), _softShadows(softShadows ? 1.0f : 0.0f)
	{ }

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

							typename PhotonMapType::PhotonType photon;
							photon.SetRadiance(emission / static_cast<float>(photonCount));
							photon.SetDirection(glm::normalize(Vector(GetRandom(), GetRandom(), GetRandom())));
							photon.SetPosition(sphere->GetPosition() + photon.GetDirection() * (sphere->GetRadius() + 0.00001f));

							SendPhoton(scene, photon);
						}
					}));
				}

				for (auto& t : tasks)
					t.wait();

				tasks.clear();
			}
		}

		for (auto& p : _photonMap)
		{
			tasks.push_back(std::async([&]()
			{
				p.second.BuildIndex();
			}));

			if (tasks.size() >= threads * 2)
			{
				for (auto& t : tasks)
					t.wait();

				tasks.clear();
			}
		}

		for (auto& task : tasks)
		{
			task.wait();
		}
	}

	Color GetRadiance(const Ray& ray) const
	{
		return Radiance(ray, 0, 1);
	}

private:
	static constexpr auto minDepth = 5;
	static constexpr auto photonCount = 1000000;
	static constexpr auto debugEpsilon = 0.05f;

	float GetRandom() const { static std::default_random_engine rnd; return _rng(rnd); }

	void SendPhoton(const Scene& scene, typename PhotonMapType::PhotonType photon, int depth = 1)
	{
		Ray ray(photon.GetPosition(), photon.GetDirection());
		IntersectionInfo info;

		if (!scene.Intersect(ray, info)) return;

		const auto& obj = *info.geometry;

		auto color = obj.GetMaterial().GetColor();
		auto glossiness = obj.GetMaterial().GetGlossiness();
		auto emission = obj.GetMaterial().GetEmission();

		/* Maximum RGB reflectivity for Russian Roulette */
		float p = glm::max(photon.GetRadiance().x, glm::max(photon.GetRadiance().y, photon.GetRadiance().z));

		if (depth > minDepth || p < 0.0001f)   /* After 5 bounces or if max reflectivity is zero */
		{
			if (drand48() < p*0.9)            /* Russian Roulette */
				photon.SetRadiance(photon.GetRadiance() * (1.f / p));        /* Scale estimator to remain unbiased */
			else
				return;  /* No further bounces */
		}

		photon.SetPosition(info.hitpoint);

		switch (obj.GetMaterial().GetReflectionType())
		{
			case DIFF:
			{
				if (depth >= 2)
				{
					_photonMap[reinterpret_cast<size_t>(info.geometry)].Add(photon);
				}

				auto nl = info.normal;

				/* Obtain flipped normal, if object hit from inside */
				if (glm::dot(info.normal, info.ray.GetDirection()) >= 0.f)
					nl = nl*-1.0f;

				/* Compute random reflection vector on hemisphere */
				auto r1 = 2.0f * PI * drand48();
				auto r2 = drand48();
				auto r2s = sqrt(r2);

				/* Set up local orthogonal coordinate system u,v,w on surface */
				auto w = nl;
				Vector u;

				if (fabs(w.x) > .1f)
					u = Vector(0.0f, 1.0f, 0.0f);
				else
					u = glm::normalize(glm::cross(Vector(1.0f, 0.0f, 0.0f), w));

				auto v = glm::cross(w, u);

				/* Random reflection vector d */
				auto d = u * cos(r1) * r2s +
					v * sin(r1) * r2s +
					w * sqrt(1.f - r2);

				photon.SetDirection(d);
				photon.SetRadiance(photon.GetRadiance() * color * std::max(0.0f, -glm::dot(info.normal, -photon.GetDirection())) * 0.7f);
				
				SendPhoton(scene, photon, depth + 1);

				break;
			}
			case SPEC:
			{
				auto direction = glm::reflect(photon.GetDirection(), info.normal);
				varyVector(direction, glossiness);
				photon.SetDirection(direction);

				SendPhoton(scene, photon, depth + 1);
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
					auto direction = glm::reflect(photon.GetDirection(), info.normal);
					varyVector(direction, glossiness);

					photon.SetDirection(direction);
					SendPhoton(scene, photon, depth + 1);

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

				/* Determine R0 for Schlick큦 approximation */
				float a = nt - nc;
				float b = nt + nc;
				float R0 = a*a / (b*b);

				/* Cosine of correct angle depending on outside/inside */
				float c;
				if (into)
					c = 1 + ddn;
				else
					c = 1 - glm::dot(tdir, normal);

				/* Compute Schlick큦 approximation of Fresnel equation */
				float Re = R0 + (1.f - R0) *c*c*c*c*c;   /* Reflectance */
				float Tr = 1.f - Re;                     /* Transmittance */

															/* Probability for selecting reflectance or transmittance */
				float P = .25f + .5f * Re;
				float RP = Re / P;         /* Scaling factors for unbiased estimator */
				float TP = Tr / (1.f - P);

				if (depth <= minDepth)   /* Initially both reflection and trasmission */
				{
					photon.SetDirection(perfectReflective.GetDirection());
					SendPhoton(scene, photon, depth + 1);

					photon.SetDirection(tdir);
					SendPhoton(scene, photon, depth + 1);
				}
				else             /* Russian Roulette */
				{
					if (drand48() < P)
					{
						photon.SetDirection(perfectReflective.GetDirection());

						SendPhoton(scene, photon, depth + 1);

						photon.SetRadiance(photon.GetRadiance() * RP);
					}
					else
					{

						photon.SetDirection(tdir);

						SendPhoton(scene, photon, depth + 1);

						photon.SetRadiance(photon.GetRadiance() * TP);
					}
				}

				break;
			}
		}
	}

	Color DisplayPhoton(const IntersectionInfo& intersectionInfo) const
	{
		if (intersectionInfo.geometry->GetMaterial().GetEmission() != Vector(0, 0, 0))
			return intersectionInfo.geometry->GetMaterial().GetEmission();

		auto it = _photonMap.find(reinterpret_cast<size_t>(intersectionInfo.geometry));
		if (it == _photonMap.end()) return Color(0, 0, 0);

		if (it->second.Empty()) return Color(0, 0, 0);

		auto color = Color(0, 0, 0);

		auto neighbour = it->second.GetNeighbours(intersectionInfo.hitpoint, 1);
		if (neighbour.GetWorstDistanceSq() < debugEpsilon && neighbour.GetIndices()[0] >= 0)
			return it->second.GetPhotons()[neighbour.GetIndices()[0]].GetRadiance();

		return Color(0, 0, 0);
	}

	Color GetPhotonMapRadianceEstimate(const IntersectionInfo& intersectionInfo) const
	{
		if (_debug) return DisplayPhoton(intersectionInfo);

		if (intersectionInfo.geometry->GetMaterial().GetEmission() != Vector(0, 0, 0))
			return intersectionInfo.geometry->GetMaterial().GetEmission();

		auto it = _photonMap.find(reinterpret_cast<size_t>(intersectionInfo.geometry));
		if (it == _photonMap.end()) return Color(0, 0, 0);

		if (it->second.Empty()) return Color(0, 0, 0);

		auto color = Color(0, 0, 0);

		auto neighbours = it->second.GetNeighbours(intersectionInfo.hitpoint, 200);

		float worstDistance = glm::sqrt(neighbours.GetWorstDistanceSq());
		float k = 1;

		auto normalizationCoeff = 1 - 2 / (3 * k);

		for (auto& idx : neighbours.GetIndices())
		{
			if (idx < 0) continue;

			auto& photon = it->second.GetPhotons()[idx];
			auto weight = std::max(0.0f, -glm::dot(intersectionInfo.normal, photon.GetDirection()));
			weight *= 1 - glm::length(intersectionInfo.hitpoint - photon.GetPosition()) / (worstDistance * k);/// exposure;
			color += photon.GetRadiance() * weight;
		}

		return color / (neighbours.GetWorstDistanceSq() * PI * normalizationCoeff);
	}

	Color Radiance(const Ray& ray, int depth, float E) const
	{
		depth++;

		IntersectionInfo intersectionInfo;
		if (!_scene.Intersect(ray, intersectionInfo)) return Color(0, 0, 0);

		const auto &obj = *intersectionInfo.geometry;

		auto glossiness = obj.GetMaterial().GetGlossiness();
		auto emission = obj.GetMaterial().GetEmission() * E;
		auto col = obj.GetMaterial().GetColor();

		/* Maximum RGB reflectivity for Russian Roulette */
		float p = glm::max(col.x, glm::max(col.y, col.z));

		if (depth > 5 || !p)   /* After 5 bounces or if max reflectivity is zero */
		{
			if (drand48() < p*0.9)            /* Russian Roulette */
				col = col * (1.0f / p);        /* Scale estimator to remain unbiased */
			else
			{
				return emission;  /* No further bounces, only return potential emission */
			}
		}

		if (obj.GetMaterial().GetReflectionType() == DIFF)
		{
			Color directLighting(0, 0 ,0);
			for(const auto& geom : _scene.GetGeometry())
			{
				if(geom->GetMaterial().HasEmission() && geom->GetType() == eGeometryType_Sphere)
				{
					auto sphere = static_cast<Sphere*>(geom.get());

					auto randomVec = glm::normalize(Vector(GetRandom(), GetRandom(), GetRandom())) * _softShadows;
					auto origin = sphere->GetPosition() + sphere->GetRadius() * randomVec;
					auto dir = glm::normalize(intersectionInfo.hitpoint - origin);
					
					Ray shadowRay(origin, dir);

					IntersectionInfo shadowInfo;
					if(_scene.Intersect(shadowRay, shadowInfo))
					{
						while(shadowInfo.geometry == geom.get())
						{
							shadowRay = Ray(shadowInfo.hitpoint +  0.0001f * dir, dir);
							if (!_scene.Intersect(shadowRay, shadowInfo)) break;
						}

						if (shadowInfo.geometry == intersectionInfo.geometry)
						{
							auto weight = std::max(0.0f, -glm::dot(intersectionInfo.normal, shadowRay.GetDirection()));
							weight *= weight;

							//directly lighted
							auto lightColor = geom->GetMaterial().GetColor();

							directLighting += col * glm::clamp(geom->GetMaterial().GetEmission(), Color(0, 0, 0), Color(1, 1, 1)) * weight;
						}
					}
				}
			}

			auto indirectLighting = col * GetPhotonMapRadianceEstimate(intersectionInfo);

			return emission + directLighting + indirectLighting;
		}

		if (obj.GetMaterial().GetReflectionType() == SPEC)
		{
			auto reflt = glm::reflect(intersectionInfo.ray.GetDirection(), intersectionInfo.normal);
			varyVector(reflt, glossiness);

			return emission + col * Radiance(Ray(intersectionInfo.hitpoint, reflt), depth, 1);
		}

		if (obj.GetMaterial().GetReflectionType() == TRAN)
		{
			auto hitpoint = intersectionInfo.hitpoint;
			auto normal = intersectionInfo.normal;
			auto origDir = intersectionInfo.ray.GetDirection();

			Ray perfectReflective(hitpoint, glm::reflect(origDir, normal));

			auto nl = normal;

			/* Obtain flipped normal, if object hit from inside */
			if (glm::dot(normal, origDir) >= 0.f)
				nl = nl*-1.0f;

			/* Otherwise object transparent, i.e. assumed dielectric glass material */

			auto into = glm::dot(normal, nl) > 0;       /* Bool for checking if ray from outside going in */
			auto nc = 1.f;                        /* Index of refraction of air (approximately) */
			auto nt = 1.5f;                      /* Index of refraction of glass (approximately) */
			float nnt;

			if (into)      /* Set ratio depending on hit from inside or outside */
				nnt = nc / nt;
			else
				nnt = nt / nc;

			auto ddn = glm::dot(origDir, nl);
			auto cos2t = 1.f - nnt * nnt * (1.f - ddn*ddn);

			/* Check for total internal reflection, if so only reflect */
			if (cos2t < 0)
			{
				return emission + col * (Radiance(perfectReflective, depth, 1));
			}

			/* Otherwise reflection and/or refraction occurs */
			Vector tdir;

			/* Determine transmitted ray direction for refraction */
			if (into)
				tdir = (origDir * nnt - normal * (ddn * nnt + sqrt(cos2t)));
			else
				tdir = (origDir * nnt + normal * (ddn * nnt + sqrt(cos2t)));

			tdir = glm::normalize(tdir);

			/* Determine R0 for Schlick큦 approximation */
			auto a = nt - nc;
			auto b = nt + nc;
			auto R0 = a*a / (b*b);

			/* Cosine of correct angle depending on outside/inside */
			float c;
			if (into)
				c = 1 + ddn;
			else
				c = 1 - glm::dot(tdir, normal);

			/* Compute Schlick큦 approximation of Fresnel equation */
			auto Re = R0 + (1.f - R0) *c*c*c*c*c;   /* Reflectance */
			auto Tr = 1.f - Re;                     /* Transmittance */

													/* Probability for selecting reflectance or transmittance */
			auto P = .25f + .5f * Re;
			auto RP = Re / P;         /* Scaling factors for unbiased estimator */
			auto TP = Tr / (1.f - P);

			if (depth < 3)   /* Initially both reflection and trasmission */
			{
				return emission + col*(Radiance(perfectReflective, depth, 1) * Re +
					Radiance(Ray(hitpoint, tdir), depth, 1) * Tr);
			}

			/* Russian Roulette */
			if (drand48() < P)
			{
				return emission + col*(Radiance(perfectReflective, depth, 1) * RP);
			}

			return emission + col*(Radiance(Ray(hitpoint, tdir), depth, 1) * TP);
		}

		assert("Unknown material");

		return Color(0, 0, 0);
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

	std::map<size_t, PhotonMapType> _photonMap;

	
	std::uniform_real_distribution<float> _rng = std::uniform_real_distribution<float>(-1.0f, 1.0f);
	bool _debug;
	const Scene& _scene;
	float _softShadows;
};

#endif // PathTracingWithPhotonMappingRP

