=========================================================
Assignment 1

Group: Bernhard Fritz, Mathias Hölzl, Florian Tischler
---------------------------------------------------------
How to build:
- cmake:
	call cmake in the same directory in which the CMakeLists.txt file is located

	> cmake ./

	this will generate a make file Makefile in the same directory

	> make 
	
	use make to build the program

	> ./Assignment1 

	to start the program

- or alternatively use the provided make file (this make file was built using the zid-gpl)
--------------------------------------------------------
How to use:

All configurable parameters can be found in the main function (Radiosity.cpp)
	- divisions must be a number between 1 and 4
	- there is another interpolation method commented out in the radiance function 
	  (RendererTriangle.h). This method uses the approach you suggested to us in the
	  last proseminar. Maybe you could take a look at it if we failed to implement it
	  correctly or if this method just does not apply to this issue.
	- interestingly we achieved a better rendering result when using fixed sample points
	  instead of random samples. See the form factor calculation function in RendererTriangle.h
	  for further information
--------------------------------------------------------
Comparison:

Because of the non trivial task of interpolating the sub triangle colors we noticed that for
eg. the shadow cast by the cuboid is less noticable than in the rectangle based rendering.

--------------------------------------------------------
Notes: 

The provided images where produced with 
	- width = 320
	- height = 240
	- division = 2
	- mcSamples = 20
	- samples = 2
	- iterations = 40

