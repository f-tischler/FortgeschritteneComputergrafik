#ifndef Math_H_included
#define Math_H_included

#include <cmath>
#include <cassert>
#include <sstream>
#include <iomanip>

/*------------------------------------------------------------------
| Struct for standard vector operations in 3D
| (used for points, vectors, and colors)
------------------------------------------------------------------*/
class Vector
{
public:
	double x, y, z;           /* Position XYZ or color RGB */

	Vector(const Vector &b) : x(b.x), y(b.y), z(b.z) {}
	Vector(double x_ = 0, double y_ = 0, double z_ = 0) : x(x_), y(y_), z(z_) {}

	Vector& operator+=(const Vector &b)
	{
		x += b.x;
		y += b.y;
		z += b.z;
		return *this;
	}

	Vector operator+(const Vector &b) const
	{
		return Vector(x + b.x, y + b.y, z + b.z);
	}

	Vector operator-(const Vector &b) const
	{
		return Vector(x - b.x, y - b.y, z - b.z);
	}

	Vector operator/(double c) const
	{
		return Vector(x / c, y / c, z / c);
	}

	Vector operator*(double c) const
	{
		return Vector(x * c, y * c, z * c);
	}

	friend Vector operator*(double c, const Vector &b)
	{
		return b * c;
	}

	Vector MultComponents(const Vector &b) const
	{
		return Vector(x * b.x, y * b.y, z * b.z);
	}

	const double LengthSquared() const
	{
		return x*x + y*y + z*z;
	}

	const double Length() const
	{
		return sqrt(LengthSquared());
	}

	const Vector Normalized() const
	{
		return Vector(x, y, z) / sqrt(x*x + y*y + z*z);
	}

	const double Dot(const Vector &b) const
	{
		return x * b.x + y * b.y + z * b.z;
	}

	const Vector Cross(const Vector &b) const
	{
		return Vector((y * b.z) - (z * b.y),
			(z * b.x) - (x * b.z),
			(x * b.y) - (y * b.x));
	}

	std::string toString() const
	{
		std::stringstream stream;
		stream << std::setprecision(3);
		stream << "(" << x << ", " << y << ", " << z << ")";
		return stream.str();
	}

	static bool AreEqual(const Vector& lhs, const Vector& rhs, const Vector& epsilon)
	{
		assert(epsilon.x > 0 && "epsilon values must be greater than zero");
		assert(epsilon.y > 0 && "epsilon values must be greater than zero");
		assert(epsilon.z > 0 && "epsilon values must be greater than zero");

		return  (std::abs(lhs.x - rhs.x) < epsilon.x) &&
				(std::abs(lhs.y - rhs.y) < epsilon.y) &&
				(std::abs(lhs.z - rhs.z) < epsilon.z);
	}
};

typedef Vector Color;

/******************************************************************
* Helper functions for smooth bicubic (Catmull-Rom) interpolation
* using 4x4 color patches;
* First interpolate in y, followed by interpolation of results in x
*******************************************************************/
inline Color cubicInterpolate(Color p[4], double x)
{
	return p[1] + 0.5 * x * (p[2] - p[0] + x * (2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] +
		x * (3.0*(p[1] - p[2]) + p[3] - p[0])));
}

inline Color bicubicInterpolate(Color p[4][4], double x, double y)
{
	Color arr[4];

	arr[0] = cubicInterpolate(p[0], y);
	arr[1] = cubicInterpolate(p[1], y);
	arr[2] = cubicInterpolate(p[2], y);
	arr[3] = cubicInterpolate(p[3], y);

	return cubicInterpolate(arr, x);
}


#endif