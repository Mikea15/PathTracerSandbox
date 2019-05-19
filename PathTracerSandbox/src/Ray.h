
#pragma once

#include "Vec3.h"

struct Ray
{
	Ray(Vec3 o, Vec3 d)
		: origin(o)
		, direction(d)
	{ }

	Vec3 origin;
	Vec3 direction;
};
