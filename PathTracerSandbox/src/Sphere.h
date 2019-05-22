
#pragma once

#include "Ray.h"

struct Sphere
{
	double radius;
	double position[3];
	double emission[3];
	double color[3];
	ReflectionType::Enum reflectionType;
	bool isLight;

	inline Vec3 GetPosition(void) const
	{
		return Vec3(position[0], position[1], position[2]);
	}

	inline Vec3 GetEmission(void) const
	{
		return Vec3(emission[0], emission[1], emission[2]);
	}

	inline Vec3 GetColor(void) const
	{
		return Vec3(color[0], color[1], color[2]);
	}

	// returns distance, 0 if nothing was hit
	inline double intersect(const Ray& r) const
	{
		// Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
		const Vec3 op = GetPosition() - r.origin;

		const double eps = 1e-5;
#if USE_OPTIMIZED_VEC 
		const double b = dot_product(op, r.direction);
		const double det = b * b - dot_product(op, op) + radius * radius;
#else
		const double b = op.dot(r.direction);
		const double det = b * b - op.dot(op) + radius * radius;
#endif
		if (det < 0)
		{
			return 0;
		}

		{
			const double t  = b - sqrt(det);
			if (t > eps)
			{
				return t;
			}
		}
		{
			const double t  = b + sqrt(det);
			if (t > eps)
			{
				return t;
			}
		}

		return 0;
	}
};
