#ifndef Renderer_H_included
#define Renderer_H_included

#include "Image.hpp"
#include "Rectangle.hpp"

class Renderer
{
protected:
	unsigned int mWidth;
	unsigned int mHeight;
	unsigned int mSamples;
	Color mClearColor;

public:
	virtual ~Renderer()
	{
	}

	Renderer(unsigned int width, unsigned int height, unsigned int samples) :
		mWidth(width),
		mHeight(height),
		mSamples(samples),
		mClearColor(0,1,1)
	{
	}

	virtual void buildScene(const std::vector<Rectangle>& rects) = 0;

	virtual void Render(Image& normal, Image& smoothed, size_t divisions, size_t mcSamples, size_t iterations) = 0;
};

#endif