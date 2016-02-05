#ifndef Ray_H_included
#define Ray_H_included

/*------------------------------------------------------------------
| Structure for rays (e.g. viewing ray, ray tracing)
------------------------------------------------------------------*/

using Vector = glm::vec3;

class Ray
{
public:
	Ray() : _org(0, 0, 0), _dir(0, 0, 0) {}
	Ray(const Vector& org, const Vector& dir) : _org(org), _dir(dir) {}

	const auto& GetOrigin() const { return _org; }
	const auto& GetDirection() const { return _dir; }

private:
	Vector _org, _dir;
};


#endif