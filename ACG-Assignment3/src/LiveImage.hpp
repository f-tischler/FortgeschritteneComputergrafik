#ifndef LiveImage_H_included
#define LiveImage_H_included


#include "Common.hpp"

#include "opengl/include/GL/freeglut.h"
#include <GL/gl.h>


#include <cmath>

unsigned int window_width = 512, window_height = 512;
const int size=window_width*window_height;

struct rgbf {float r; float g; float b;};
//WBL 9 May 2007 Based on
//http://www.codeguru.com/cpp/w-d/dislog/commondialogs/article.php/c1861/
//Common.h
void toRGBf(const float h, const float s, const float v,
	    rgbf* rgb)
{
  /*
RGBType rgb;
	if(!h  && !s)
	{
		rgb.r = rgb.g = rgb.b = v;
	}
  */
  //rgbf* rgb = (rgbf*) out;
double min,max,delta,hue;
	
	max = v;
	delta = max * s;
	min = max - delta;

	hue = h;
	if(h > 300 || h <= 60)
	{
		rgb->r = max;
		if(h > 300)
		{
			rgb->g = min;
			hue = (hue - 360.0)/60.0;
			rgb->b = ((hue * delta - min) * -1);
		}
		else
		{
			rgb->b = min;
			hue = hue / 60.0;
			rgb->g = (hue * delta + min);
		}
	}
	else
	if(h > 60 && h < 180)
	{
		rgb->g = max;
		if(h < 120)
		{
			rgb->b = min;
			hue = (hue/60.0 - 2.0 ) * delta;
			rgb->r = min - hue;
		}
		else
		{
			rgb->r = min;
			hue = (hue/60 - 2.0) * delta;
			rgb->b = (min + hue);
		}
	}
	else
	{
		rgb->b = max;
		if(h < 240)
		{
			rgb->r = min;
			hue = (hue/60.0 - 4.0 ) * delta;
			rgb->g = (min - hue);
		}
		else
		{
			rgb->g = min;
			hue = (hue/60 - 4.0) * delta;
			rgb->r = (min + hue);
		}
	}
}


//Convert a wide range of data values into nice colours 
void colour(const float data, float* out) {
  //convert data to angle
  const float a = atan2(data,1)/(2*atan2(1,1)); // -1 .. +1
  const float angle = (1+a)*180; //red=0 at -1,+1

  const float saturation = 1;

  const float h = (data<-1||data>1)? 1 : fabs(data);

  toRGBf(angle,saturation,h,(rgbf*)out);
}


#include <mutex>

class LiveImage;
std::vector<LiveImage*> images;

class Renderer
{
	Renderer()
	{

	}
};


class LiveImage
{
	int w;
	int h;

	std::mutex mGuard;
	std::vector<Color> mPixels;

public:
	LiveImage(int w, int h) : w(w), h(h)
	{
		mPixels.resize(w*h);

		for (size_t i = 0; i < w*h; i++)
			mPixels[i] = Color(1.f,0.f,1.f);

		images.push_back(this);
	}

	void draw()
	{
		std::lock_guard<std::mutex> lock(mGuard);
		glDrawPixels(w, h, GL_RGB, GL_FLOAT, mPixels.data());
	}

	void show()
	{
		glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
		glutInitWindowSize(w, h);
		glutCreateWindow("OpenGL glDrawPixels demo");

		glutDisplayFunc([]
		{
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

			for (auto img : images)
				img->draw();

			glutSwapBuffers();
		});


		glEnable(GL_DEPTH_TEST);
		glClearColor(0.0, 0.0, 0.0, 1.0);
	}

	void draw(Image& img)
	{

	}

	void set(int x, int y, Color c)
	{
		static long t = 0;
		if (t++ % 500 == 0)
		{
			glutMainLoopEvent();
			glutPostRedisplay();
		}

		std::lock_guard<std::mutex> lock(mGuard);
		mPixels[y*w + x] = c;
	}
};


#endif 
