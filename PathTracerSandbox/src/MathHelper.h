
#pragma once

#include <cmath>
#include <algorithm>

static const double PI = 3.14159265358979323846264338327950288;
static const double ONE_OVER_PI = 0.318309886183790671538;

namespace MathUtils
{
	inline double Clamp(double x)
	{
		return x < 0 ? 0 : x > 1 ? 1 : x;
	}

	inline float Clamp(float value, float min, float max)
	{
		if (value > max) { value = max; }
		if (value < min) { value = min; }
		return value;
	}

	inline int Gamma(double x)
	{
		const double clamped = Clamp(x);

		const double scale = 255;

		return static_cast<int>(0.5 + scale * 1.138 * sqrt(clamped) - scale * 0.138 * clamped);
	}

	inline unsigned int Linear2sRGB(float color, float gamma)
	{
		const float gcolor = std::pow(color, 1.0f / gamma);
		return static_cast<unsigned int>(Clamp(255.0f * gcolor, 0.0f, 255.0f));
	}

	inline float Reflectance0(float n1, float n2)
	{
		const float sqrt_r0 = (n1 - n2) / (n1 + n2);
		return sqrt_r0 * sqrt_r0;
	}

	inline float Schlick(float cosine, float refIndex)
	{
		float r0 = (1 - refIndex) / (1 + refIndex);
		r0 = r0 * r0;
		return r0 + (1 - r0) * std::pow((1 - cosine), 5);
	}

	inline float SchlickReflectance(float n1, float n2, float c)
	{
		const float r0 = Reflectance0(n1, n2);
		return r0 + (1.0f - r0) * c * c * c * c * c;
	}
}
