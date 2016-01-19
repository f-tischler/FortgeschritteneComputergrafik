#ifndef SimpleRadianceProvider_H
#define SimpleRadianceProvider_H

#include "Raytracer.hpp"

class SimpleRadianceProvider
{
public:
	Color GetRadiance(const IntersectionInfo& intersectionInfo, TracingInfo& tracingInfo) const
	{
		auto color = intersectionInfo.geometry->GetMaterial().GetColor();
		
		tracingInfo.spawnedRays.clear();

		return color;
	}
};

#endif // SimpleRadianceProvider_H

