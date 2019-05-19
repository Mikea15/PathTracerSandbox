
#pragma once

#include <math.h>
#include "XorShift.h"

#include "vectorclass/vectorclass.h"
#include "vectorclass/vector3d.h"

// using SIMD makes it slower :/ as expected.
#define USE_OPTIMIZED_VEC 1

struct ReflectionType
{
	enum Enum
	{
		DIFFUSE,
		SPECULAR,
		REFRACTIVE
	};
};

#if USE_OPTIMIZED_VEC 

using Vec3 = Vec3d;

#else 

struct Vec3
{
	double x, y, z;

	inline Vec3()
		: x(0.0), y(0.0), z(0.0)
	{
	}

	inline Vec3(double x_, double y_, double z_)
	{
		x = x_;
		y = y_;
		z = z_;
	}

	inline Vec3 operator+(const Vec3& other) const
	{
		return Vec3(x + other.x, y + other.y, z + other.z);
	}

	inline Vec3 operator-(const Vec3& other) const
	{
		return Vec3(x - other.x, y - other.y, z - other.z);
	}

	inline Vec3 operator*(double scalar) const
	{
		return Vec3(x * scalar, y * scalar, z * scalar);
	}

	inline Vec3 mult(const Vec3& other) const
	{
		return Vec3(x * other.x, y * other.y, z * other.z);
	}

	inline double length() const
	{
		return sqrt(x * x + y * y + z * z);
	}

	inline double lengthSq() const
	{
		return x * x + y * y + z * z;
	}

	inline Vec3 normalize(void) const
	{
		const double oneByLength = 1 / length();

		return Vec3(x*oneByLength, y*oneByLength, z*oneByLength);
	}

	inline double dot(const Vec3& b) const
	{
		return x * b.x + y * b.y + z * b.z;
	}

	inline Vec3 cross(const Vec3& other) const
	{
		return Vec3(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x);
	}
};

#endif

static Vec3 RandomInUnitSphere()
{
	Vec3 p(0, 0, 0);
#if USE_OPTIMIZED_VEC 
	do {
		p = Vec3(random::XorShift(0.0f, 1.0f), random::XorShift(0.0f, 1.0f), random::XorShift(0.0f, 1.0f)) * 2.0f - Vec3(1, 1, 1);
	} while (vector_length(p) >= 1.0f);
#else
	do {
		p = Vec3(random::XorShift(0.0f, 1.0f), random::XorShift(0.0f, 1.0f), random::XorShift(0.0f, 1.0f)) * 2.0f - Vec3(1, 1, 1);
	} while (p.lengthSq() >= 1.0f);
#endif
	return p;
}

static Vec3 RandomInUnitHemisphere(const Vec3& n)
{
	Vec3 p = RandomInUnitSphere();
#if USE_OPTIMIZED_VEC 
	if ( dot_product(p, n) <= 0.0)
		p = p - p;
#else
	if (p.dot(n) <= 0.0)
		p = p - p;
#endif
	return p;
}
