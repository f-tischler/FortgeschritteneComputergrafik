#ifndef Renderer_H_included
#define Renderer_H_included

#include "Image.hpp"

class Renderer
{
protected:
	unsigned int mWidth;
	unsigned int mHeight;
	unsigned int mSamples;
	Color mClearColor;

public:
	Renderer(unsigned int width, unsigned int height, unsigned int samples) :
		mWidth(width),
		mHeight(height),
		mSamples(samples),
		mClearColor(0,0,0)
	{
	}

	virtual void buildScene(const std::vector<Rectangle>& rects) = 0;

	virtual void Render(Image& normal, Image& smoothed) = 0;
};

#endif