#ifndef NanoFlannPhotonMap_H
#define NanoFlannPhotonMap_H

#include <mutex>
#include <vector>
#include "nanoflann.hpp"
#include "Common.hpp"

class Photon
{
public:
	void SetPosition(const Vector& position) { _position = position; }
	void SetRadiance(const Color& radiance) { _radiance = radiance; }
	void SetDirection(const Vector& direction) { _direction = direction; }

	const auto& GetPosition() const { return _position; }
	const auto& GetRadiance() const { return _radiance; }
	const auto& GetDirection() const { return _direction; }

private:
	Color _radiance;
	Vector _position;
	Vector _direction;
};

class NeighbourSet
{
public:
	explicit NeighbourSet(size_t n)
		: _indices(n, -1),
		  _distances(n, -1),
		  _worstDistanceSq(-1)
	{
	}

	const auto& GetIndices() const { return _indices; }
	const auto& GetDistances() const { return _distances; }
	
	auto& GetIndices() { return _indices; }
	auto& GetDistances() { return _distances; }

	auto GetWorstDistanceSq() { return _worstDistanceSq; }
	void SetWorstDistanceSq(float distanceSq) { _worstDistanceSq = distanceSq; }

private:
	std::vector<size_t> _indices;
	std::vector<float> _distances;
	float _worstDistanceSq;
};

class NanoFlannPhotonMap
{
public:
	using PhotonType = Photon;

	NanoFlannPhotonMap() : _data(0), _index(3, *this, nanoflann::KDTreeSingleIndexAdaptorParams(10)) { }

	size_t kdtree_get_point_count() const { return _data.size(); }

	auto kdtree_distance(const float* p1, const size_t idx_p2, size_t) const
	{
		const auto d0 = p1[0] - _data[idx_p2].GetPosition().x;
		const auto d1 = p1[1] - _data[idx_p2].GetPosition().y;
		const auto d2 = p1[2] - _data[idx_p2].GetPosition().z;
		return d0*d0 + d1*d1 + d2*d2;
	}

	auto kdtree_get_pt(const size_t idx, int dim) const
	{
		if (dim == 0) return _data[idx].GetPosition().x;
		if (dim == 1) return _data[idx].GetPosition().y;

		return _data[idx].GetPosition().z;
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
		const float queryPos[3] = { pos.x, pos.y, pos.z };

		NeighbourSet set(n);
		nanoflann::KNNResultSet<float> resultSet(n);

		resultSet.init(set.GetIndices().data(), set.GetDistances().data());
		_index.findNeighbors(resultSet, &queryPos[0], nanoflann::SearchParams(10));

		set.SetWorstDistanceSq(resultSet.worstDist());

		return set;
	}

	const auto& GetPhotons() const { return _data; }

	void BuildIndex()
	{
		if (!_data.empty())
			_index.buildIndex();
	}

	bool Empty() const { return _data.empty(); }

private:
	using MyKdTreeType =
		nanoflann::KDTreeSingleIndexAdaptor<
		nanoflann::L2_Simple_Adaptor<float, NanoFlannPhotonMap>,
		NanoFlannPhotonMap,
		3>;

	std::vector<Photon> _data;
	MyKdTreeType _index;
	std::mutex _dataMtx;
};

#endif // NanoFlannPhotonMap_H
