#ifndef Scene_H_included
#define Scene_H_included

#include <vector>
#include "Rectangle.hpp"

inline std::vector<Rectangle> getScene3() {
    static Rectangle rects[] = {
        /* Cornell Box walls */
        //Rectangle(Vector(0.0, 0.0, 0.0), Vector(100.0, 0.0, 0.0), Vector(0.0, 80.0, 0.0),
        //          Vector(), Color(0.75, 0.75, 0.75)), /* Back */
        Rectangle(Vector(0.0, 0.0, 170.0), Vector(100.0, 0.0, 0.0), Vector(0.0, 0.0, -170.0),
                  Vector(), Color(0.75, 0.75, 0.75)), /* Bottom */
        // Rectangle(Vector(0.0, 80.0, 0.0), Vector(100.0, 0.0, 0.0), Vector(0.0, 0.0, 170.0),
        //          Vector(), Color(0.75, 0.75, 0.75)), /* Top */
        //Rectangle(Vector(0.0, 0.0, 170.0), Vector(0.0, 0.0, -170.0), Vector(0.0, 80.0, 0.0),
        //          Vector(), Color(0.75, 0.25, 0.25)), /* Left */
        //Rectangle(Vector(100.0, 0.0, 0.0), Vector(0.0, 0.0, 170.0), Vector(0.0, 80.0, 0.0),
        //          Vector(), Color(0.25, 0.25, 0.75)), /* Right */
        //Rectangle(Vector(100.0, 0.0, 170.0), Vector(-100.0, 0.0, 0.0), Vector(0.0, -80.0, 0.0),
        //          Vector(), Color(0, 1, 0)), /* Front (not visible) */
        
        /* Area light source on top */
        Rectangle(Vector(40.0, 79.99, 65.0), Vector(20.0, 0.0, 0.0), Vector(0.0, 0.0, 20.0),
                  Vector(12, 12, 12), Color(0.75, 0.75, 0.75))
    };
    return std::vector<Rectangle>(rects, rects + sizeof(rects) / sizeof(rects[0]));
}


inline std::vector<Rectangle> getScene2()
{
	static Rectangle rects[] =
	{
		/* Cornell Box walls */
		Rectangle(Vector(0.0, 0.0, 0.0), Vector(100.0, 0.0, 0.0), Vector(0.0, 80.0, 0.0),
		Vector(), Color(0.75, 0.75, 0.75)), /* Back */
		Rectangle(Vector(0.0, 0.0, 170.0), Vector(100.0, 0.0, 0.0), Vector(0.0, 0.0, -170.0),
		Vector(), Color(0.75, 0.75, 0.75)), /* Bottom */
		Rectangle(Vector(0.0, 80.0, 0.0), Vector(100.0, 0.0, 0.0), Vector(0.0, 0.0, 170.0),
		Vector(), Color(0.75, 0.75, 0.75)), /* Top */
		Rectangle(Vector(0.0, 0.0, 170.0), Vector(0.0, 0.0, -170.0), Vector(0.0, 80.0, 0.0),
		Vector(), Color(0.75, 0.25, 0.25)), /* Left */
		Rectangle(Vector(100.0, 0.0, 0.0), Vector(0.0, 0.0, 170.0), Vector(0.0, 80.0, 0.0),
		Vector(), Color(0.25, 0.25, 0.75)), /* Right */
		Rectangle(Vector(100.0, 0.0, 170.0), Vector(-100.0, 0.0, 0.0), Vector(0.0, -80.0, 0.0),
		Vector(), Color(0, 1, 0)), /* Front (not visible) */

		/* Area light source on top */
		Rectangle(Vector(40.0, 79.99, 65.0), Vector(20.0, 0.0, 0.0), Vector(0.0, 0.0, 20.0),
		Vector(12, 12, 12), Color(0.75, 0.75, 0.75)),

		/* Cuboid in room */
		Rectangle(Vector(30.0, 0.0, 100.0), Vector(0.0, 0.0, -20.0), Vector(0.0, 40.0, 0.0),
		Vector(), Color(0.75, 0.75, 0.75)), /* Right */
		Rectangle(Vector(10.0, 0.0, 80.0), Vector(0.0, 0.0, 20.0), Vector(0.0, 40.0, 0.0),
		Vector(), Color(0.75, 0.75, 0.75)), /* Left */
		Rectangle(Vector(10.0, 0.0, 100.0), Vector(20.0, 0.0, 0.0), Vector(0.0, 40.0, 0.0),
		Vector(), Color(0.75, 0.75, 0.75)), /* Front */
		Rectangle(Vector(30.0, 0.0, 80.0), Vector(-20.0, 0.0, 0.0), Vector(0.0, -40.0, 0.0),
		Vector(), Color(0.75, 0.75, 0.75)), /* Back */
		Rectangle(Vector(10.0, 40.0, 100.0), Vector(20.0, 0.0, 0.0), Vector(0.0, 0.0, -20.0),
		Vector(), Color(0.75, 0.75, 0.75)), /* Top */
	};

	return std::vector<Rectangle>(rects, rects + sizeof(rects) / sizeof(rects[0]));
}

inline std::vector<Rectangle> getScene()
{
	static Rectangle rects[] =
	{
		/* Cornell Box walls */
		Rectangle(Vector(0.0, 0.0, 0.0), Vector(100.0, 0.0, 0.0), Vector(0.0, 80.0, 0.0),
		Vector(), Color(0.75, 0, 0)), /* Back */
		Rectangle(Vector(0.0, 0.0, 170.0), Vector(100.0, 0.0, 0.0), Vector(0.0, 0.0, -170.0),
		Vector(), Color(0, 0.75, 0)), /* Bottom */
		Rectangle(Vector(0.0, 80.0, 0.0), Vector(100.0, 0.0, 0.0), Vector(0.0, 0.0, 170.0),
		Vector(), Color(0, 0.75, 0)), /* Top */
		Rectangle(Vector(0.0, 0.0, 170.0), Vector(0.0, 0.0, -170.0), Vector(0.0, 80.0, 0.0),
		Vector(), Color(0, 0, 0.25)), /* Left */
		Rectangle(Vector(100.0, 0.0, 0.0), Vector(0.0, 0.0, 170.0), Vector(0.0, 80.0, 0.0),
		Vector(), Color(0, 0, 0.25)), /* Right */
		Rectangle(Vector(100.0, 0.0, 170.0), Vector(-100.0, 0.0, 0.0), Vector(0.0, -80.0, 0.0),
		Vector(), Color(0, 1, 0)), /* Front (not visible) */

		/* Area light source on top */
		Rectangle(Vector(40.0, 79.99, 65.0), Vector(20.0, 0.0, 0.0), Vector(0.0, 0.0, 20.0),
		Vector(12, 12, 12), Color(0.75, 0.75, 0.75)),

		/* Cuboid in room */
		Rectangle(Vector(30.0, 0.0, 100.0), Vector(0.0, 0.0, -20.0), Vector(0.0, 40.0, 0.0),
		Vector(), Color(0.75, 0.75, 0.75)), /* Right */
		Rectangle(Vector(10.0, 0.0, 80.0), Vector(0.0, 0.0, 20.0), Vector(0.0, 40.0, 0.0),
		Vector(), Color(0.75, 0.75, 0.75)), /* Left */
		Rectangle(Vector(10.0, 0.0, 100.0), Vector(20.0, 0.0, 0.0), Vector(0.0, 40.0, 0.0),
		Vector(), Color(0.75, 0.75, 0.75)), /* Front */
		Rectangle(Vector(30.0, 0.0, 80.0), Vector(-20.0, 0.0, 0.0), Vector(0.0, -40.0, 0.0),
		Vector(), Color(0.75, 0.75, 0.75)), /* Back */
		Rectangle(Vector(10.0, 40.0, 100.0), Vector(20.0, 0.0, 0.0), Vector(0.0, 0.0, -20.0),
		Vector(), Color(0.75, 0.75, 0.75)), /* Top */
	};

	return std::vector<Rectangle>(rects, rects + sizeof(rects) / sizeof(rects[0]));
}

#endif