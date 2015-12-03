#ifndef BoundingBox_H_included
#define BoundingBox_H_included

#include "Ray.hpp"

template<class vec3>
class BoundingBox
{
	vec3 mMin;
	vec3 mMax;

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
		mMin.x = FLT_MAX;
		mMin.y = FLT_MAX;
		mMin.z = FLT_MAX;

		mMax.x = FLT_MIN;
		mMax.y = FLT_MIN;
		mMax.z = FLT_MIN;
	}

	void addPoint(vec3 v)
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
		if (r.dir.x >= 0)
		{
			tmin = (mMin.x - r.org.x) / r.dir.x;
			tmax = (mMax.x - r.org.x) / r.dir.x;
		}
		else 
		{
			tmin = (mMax.x - r.org.x) / r.dir.x;
			tmax = (mMin.x - r.org.x) / r.dir.x;
		}

		if (r.dir.y >= 0)
		{
			tymin = (mMin.y - r.org.y) / r.dir.y;
			tymax = (mMax.y - r.org.y) / r.dir.y;
		}
		else
		{
			tymin = (mMax.y - r.org.y) / r.dir.y;
			tymax = (mMin.y - r.org.y) / r.dir.y;
		}

		if ((tmin > tymax) || (tymin > tmax))
			return false;
		
		if (tymin > tmin)
			tmin = tymin;

		if (tymax < tmax)
			tmax = tymax;

		if (r.dir.z >= 0) {

			tzmin = (mMin.z - r.org.z) / r.dir.z;
			tzmax = (mMax.z - r.org.z) / r.dir.z;
		}
		else 
		{
			tzmin = (mMax.z - r.org.z) / r.dir.z;
			tzmax = (mMin.z - r.org.z) / r.dir.z;
		}

		if ((tmin > tzmax) || (tzmin > tmax))
			return false;

		if (tzmin > tmin)
			tmin = tzmin;

		if (tzmax < tmax)
			tmax = tzmax;

		return ((tmin < t1) && (tmax > t0));
	}
	
	const vec3& min() const { return mMin; }
	const vec3& max() const { return mMax; }

	double width()  const   { return mMax.x - mMin.x; }
	double height() const   { return mMax.y - mMin.y; }
	double depth()  const   { return mMax.z - mMin.z; }
	double diameter() const { return (mMax - mMin).Length(); }
	vec3 position() const   { return mMin + (mMax - mMin) / 2; }
};


#endif