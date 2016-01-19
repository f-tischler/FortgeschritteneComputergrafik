#ifndef Camera_H
#define Camera_H

#include "Ray.hpp"
#include "glm/common.hpp"
#include "glm/glm.hpp"

class Camera
{
public:
	explicit Camera(const Vector& pos, const Vector& lookAt, const float fovinRadians)
		: _ray(pos, glm::normalize(lookAt - pos)),
		  _lookAt(lookAt),
		  _focalDistance(glm::length(_lookAt - _ray.GetOrigin())),
		  _fov(fovinRadians)
	{
		
	}
	const auto& GetPosition() const { return _ray.GetOrigin(); }
	const auto& GetLookAt() const { return  _lookAt; }
	const auto& GetRay() const { return _ray; }

	void Move(const Vector& offset)
	{
		SetLookAt(_lookAt + offset);
		SetPosition(_ray.GetOrigin() + offset);
	}

	void SetLookAt(const Vector& lookAt)
	{
		_lookAt = lookAt;
		_ray = Ray(_ray.GetOrigin(), glm::normalize(_lookAt - _ray.GetOrigin()));
		_focalDistance = glm::length(_lookAt - _ray.GetOrigin());
	}

	void ChangeFocus(const float offset)
	{
		SetLookAt(GetLookAt() + _ray.GetDirection() * offset);
	}

	void SetPosition(const Vector& newPos)
	{
		_ray = Ray(newPos, glm::normalize(_lookAt - newPos));
		_focalDistance = glm::length(_lookAt - _ray.GetOrigin());
	}



	auto GetFocalDistance() const { return _focalDistance; }
	auto GetFov() const { return _fov; }
	
private:
	Ray _ray;
	Vector _lookAt;
	float _focalDistance;
	float _fov;
};

#endif // Camera_H

