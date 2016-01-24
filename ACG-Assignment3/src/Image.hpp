#ifndef Image_H_included
#define Image_H_included

#ifndef _MSVC
#include <stdio.h>
#endif

#include <string>
#include <vector>

#include "glm/glm.hpp"

using Color = glm::vec3;

/*------------------------------------------------------------------
| Struct holds pixels/colors of rendered image
------------------------------------------------------------------*/
class Image
{
public:
	using StorageType = std::vector<Color>;
	using SizeType = size_t;

	Image(SizeType w, SizeType h) 
		: _width(w), _height(h)
	{
		_pixels = StorageType(_width * _height, glm::vec3(0, 0, 0));
	}

	const Color& GetColor(SizeType x, SizeType y) const
	{
		auto image_index = (_height - y - 1) * _width + x;
		return _pixels[image_index];
	}

	void SetColor(SizeType x, SizeType y, const Color& color)
	{
		auto image_index = (_height - y - 1) * _width + x;
		_pixels[image_index] = color;
	}

	void AddColor(SizeType x, SizeType y, const Color& color)
	{
		auto image_index = (_height - y - 1) * _width + x;
		_pixels[image_index] = _pixels[image_index] + color;
	}

	auto GetAspect() const { return _width / static_cast<float>(_height); }
	auto GetWidth() const { return _width; }
	auto GetHeight() const { return _height; }


	void Save(const std::string &filename)
	{
		/* Save image in PPM format */
		FILE *f;

#ifdef _MSC_VER 
		fopen_s(&f, filename.c_str(), "wb");
#else
        f = fopen(filename.c_str(), "wb");
#endif
		
		fprintf(f, "P3\n%d %d\n%d\n", static_cast<unsigned>(_width), static_cast<unsigned>(_height), 255);

		for (auto i = 0ul; i < _width * _height; i++)
		{
			fprintf(f, "%d %d %d ", 
				ToInteger(_pixels[i].x),
				ToInteger(_pixels[i].y),
				ToInteger(_pixels[i].z));
		}

		fclose(f);
	}

private:
	SizeType _width, _height;
	StorageType _pixels;

	int ToInteger(double x) const
	{
		/* Clamp to [0,1] */
		if (x < 0.0)
			x = 0.0;

		if (x > 1.0)
			x = 1.0;

		/* Apply gamma correction and convert to integer */
		return int(pow(x, 1 / 2.2) * 255 + .5);
	}
};

#endif