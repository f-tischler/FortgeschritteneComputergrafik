#ifndef BoundingBox_H
#define BoundingBox_H

#include "Common.hpp"
#include "Ray.hpp"

#include <algorithm>

class BoundingBox
{
public:
	BoundingBox()
	{
		mMin.x = FLT_MAX;
		mMin.y = FLT_MAX;
		mMin.z = FLT_MAX;

		mMax.x = FLT_MIN;
		mMax.y = FLT_MIN;
		mMax.z = FLT_MIN;
	}

	void reset()
	{
		mMin.x = std::numeric_limits<float>::max();
		mMin.y = std::numeric_limits<float>::max();
		mMin.z = std::numeric_limits<float>::max();

		mMax.x = std::numeric_limits<float>::min();
		mMax.y = std::numeric_limits<float>::min();
		mMax.z = std::numeric_limits<float>::min();
	}

	void addPoint(Vector v)
	{
		mMax.x = std::max(mMax.x, v.x);
		mMax.y = std::max(mMax.y, v.y);
		mMax.z = std::max(mMax.z, v.z);

		mMin.x = std::min(mMin.x, v.x);
		mMin.y = std::min(mMin.y, v.y);
		mMin.z = std::min(mMin.z, v.z);
	}

	bool intersect(const Ray &r) const
	{
		double t0, t1;
		t0 = 0.0f;
		t1 = DBL_MAX;

		double tmin, tmax, tymin, tymax, tzmin, tzmax;
		if (r.GetDirection().x >= 0)
		{
			tmin = (mMin.x - r.GetOrigin().x) / r.GetDirection().x;
			tmax = (mMax.x - r.GetOrigin().x) / r.GetDirection().x;
		}
		else
		{
			tmin = (mMax.x - r.GetOrigin().x) / r.GetDirection().x;
			tmax = (mMin.x - r.GetOrigin().x) / r.GetDirection().x;
		}

		if (r.GetDirection().y >= 0)
		{
			tymin = (mMin.y - r.GetOrigin().y) / r.GetDirection().y;
			tymax = (mMax.y - r.GetOrigin().y) / r.GetDirection().y;
		}
		else
		{
			tymin = (mMax.y - r.GetOrigin().y) / r.GetDirection().y;
			tymax = (mMin.y - r.GetOrigin().y) / r.GetDirection().y;
		}

		if ((tmin > tymax) || (tymin > tmax))
			return false;

		if (tymin > tmin)
			tmin = tymin;

		if (tymax < tmax)
			tmax = tymax;

		if (r.GetDirection().z >= 0) {

			tzmin = (mMin.z - r.GetOrigin().z) / r.GetDirection().z;
			tzmax = (mMax.z - r.GetOrigin().z) / r.GetDirection().z;
		}
		else
		{
			tzmin = (mMax.z - r.GetOrigin().z) / r.GetDirection().z;
			tzmax = (mMin.z - r.GetOrigin().z) / r.GetDirection().z;
		}

		if ((tmin > tzmax) || (tzmin > tmax))
			return false;

		if (tzmin > tmin)
			tmin = tzmin;

		if (tzmax < tmax)
			tmax = tzmax;

		return tmin < t1 && tmax > t0;
	}

	const auto& min() const { return mMin; }
	const auto& max() const { return mMax; }

	auto width()  const { return mMax.x - mMin.x; }
	auto height() const { return mMax.y - mMin.y; }
	auto depth()  const { return mMax.z - mMin.z; }
	auto diameter() const { return glm::length(mMax - mMin); }
	auto position() const { return mMin + (mMax - mMin) / 2.0f; }

private:
	Vector mMin;
	Vector mMax;
};

#endif // BoundingBox_H
