#ifndef Camera_H
#define Camera_H

#include "Ray.hpp"

class Camera
{
public:
	explicit Camera(const Vector& pos, const Vector& lookAt)
		: _ray(pos, (lookAt - pos).Normalized()), 
		  _lookAt(lookAt), 
		  _focalDistance((_lookAt - _ray.org).Length())
	{
		
	}

	void Move(const Vector& offset)
	{
		SetLookAt(_lookAt + offset);
		SetPosition(_ray.org + offset);
	}

	void ChangeFocus(const double offset)
	{
		SetLookAt(GetLookAt() + _ray.dir * offset);
	}

	void SetPosition(const Vector& newPos)
	{
		_ray = Ray(newPos, (_lookAt - newPos).Normalized());
		_focalDistance = (_lookAt - _ray.org).Length();
	}

	void SetLookAt(const Vector& lookAt)
	{
		_lookAt = lookAt;
		_ray = Ray(_ray.org, (_lookAt - _ray.org).Normalized());
		_focalDistance = (_lookAt - _ray.org).Length();
	}

	const Vector& GetPosition() const { return _ray.org; }
	const Vector& GetLookAt() const { return  _lookAt; }
	const Ray& GetRay() const { return _ray; }

	double GetFocalDistance() const { return _focalDistance; }

private:
	Ray _ray;
	Vector _lookAt;
	double _focalDistance;
};

#endif // Camera_H

