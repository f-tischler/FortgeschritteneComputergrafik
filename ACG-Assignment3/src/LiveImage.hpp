#ifndef LiveImage_H_included
#define LiveImage_H_included

#include <mutex>
#include "Common.hpp"

#include "opengl/include/GL/freeglut.h"
#include <GL/gl.h>

class LiveImage;

class LiveImage
{
	size_t w;
	size_t h;

	std::mutex mGuard;
	std::vector<Color> mPixels;

public:
	LiveImage(size_t w, size_t h) : w(w), h(h)
	{
		mPixels.resize(w*h);

		for (size_t i = 0; i < w*h; i++)
			mPixels[i] = Color(1.f,0.f,1.f);
	}

	void show()
	{
		glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
		glutInitWindowSize(static_cast<int>(w), static_cast<int>(h));
		glutCreateWindow("Photon-Mapper");
		glEnable(GL_DEPTH_TEST);
		glClearColor(0.0, 0.0, 0.0, 1.0);
	}

	void update()
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		{
			std::lock_guard<std::mutex> lock(mGuard);
			glDrawPixels(static_cast<int>(w), static_cast<int>(h), GL_RGB, GL_FLOAT, mPixels.data());
		}
		glutSwapBuffers();

		glutMainLoopEvent();
	}

	void set(size_t x, size_t y, Color c)
	{
		std::lock_guard<std::mutex> lock(mGuard);
		mPixels[y*w + x] = c;
	}
};


#endif 
