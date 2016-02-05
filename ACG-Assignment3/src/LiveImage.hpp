#ifndef LiveImage_H_included
#define LiveImage_H_included

#include "Common.hpp"

#include "opengl/include/GL/freeglut.h"
#include <GL/gl.h>

class LiveImage;

class LiveImage
{
public:
	LiveImage(size_t w, size_t h)
		: _width(w), _height(h), _pixels(w * h, Color(1.f, 0.f, 1.f)) { }

	void Show() const
	{
		glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
		glutInitWindowSize(static_cast<int>(_width), static_cast<int>(_height));
		glutCreateWindow("Photon-Mapper");
		glEnable(GL_DEPTH_TEST);
		glClearColor(0.0, 0.0, 0.0, 1.0);
	}

	void Update() const
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glDrawPixels(static_cast<int>(_width), static_cast<int>(_height), GL_RGB, GL_FLOAT, _pixels.data());
		glutSwapBuffers();
		glutMainLoopEvent();
	}

	void Set(size_t x, size_t y, const Color& c)
	{
		_pixels[y*_width + x] = c;
	}

private:
	size_t _width;
	size_t _height;

	std::vector<Color> _pixels;
};


#endif 
