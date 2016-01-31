#ifndef SimpleRadianceProvider_H
#define SimpleRadianceProvider_H

class SimpleRadianceProvider
{
public:
	explicit SimpleRadianceProvider(const Scene& scene) : _scene(scene) {}

	Color GetRadiance(const Ray& ray) const
	{
		IntersectionInfo info;
		if(_scene.Intersect(ray, info))
		{
			return info.geometry->GetMaterial().GetColor();
		}
	
		return Color(0, 0, 0);
	}

private:
	const Scene& _scene;
};

#endif // SimpleRadianceProvider_H

