#include "Raytracer.h"

#include "XorShift.h"
#include "MathHelper.h"
#include "Vec3.h"
#include "Scene.h"
#include "Sphere.h"

Raytracer::Raytracer(int width, int height, int subSamples, int samplesPerPixel, int raysPerFrame, int aoSamples, float aoRayLength)
	: m_imgWidth(width)
	, m_imgHeight(height)
	, m_pixelCount(m_imgWidth* m_imgHeight)
	, m_invWidth(1.0f / m_imgWidth)
	, m_invHeight(1.0f / m_imgHeight)
	, m_subSamples(subSamples)
	, m_invSubSamples(1.0f / m_subSamples)
	, m_samplesPerPixel(samplesPerPixel)
	, m_invSamplesPerPixel(1.0f / m_samplesPerPixel)
	, m_raysPerFrame(raysPerFrame)
	, m_renderMode(RenderMode::SubPixelSamplingAccumulation)
	, m_renderType(RenderType::Diffuse)
	, m_enableLightSampling(false)
	, m_rayCount(0)
	, m_aoSamples(aoSamples)
	, m_invAoSamples(1.0f / aoSamples)
	, m_aoRayLength(aoRayLength)
{
	m_cols.resize(m_pixelCount, Vec3(0, 0, 0));
	m_sampleAccumulation.resize(m_pixelCount, 0);

	m_colorData.resize(m_pixelCount, Vec3());
	m_rgbData.resize(m_pixelCount * 4u, 0);
}

Raytracer::~Raytracer()
{
}

void Raytracer::Trace(const Ray& cameraRay)
{
	m_rayCount = 0;
	switch (m_renderMode)
	{
	case RenderMode::SubPixelSamplingAccumulation:
		TraceSampleAccumulation(cameraRay);
		break;
	case RenderMode::MultiSampling:
		// TraceMultiSampling(world, lights, cam);
		break;
	case RenderMode::SubPixelSampling:
		// TraceSubPixel(world, lights, cam);
		break;
	default:
		break;
	}
}

void Raytracer::Clear()
{
	for (unsigned int p = 0; p < m_pixelCount; ++p)
	{
		m_rgbData[p * 4u + 2u] = 0;
		m_rgbData[p * 4u + 1u] = 0;
		m_rgbData[p * 4u + 0u] = 0;

		m_colorData[p] = Vec3(0, 0, 0);
		m_cols[p] = Vec3(0, 0, 0);
		m_sampleAccumulation[p] = 0;
	}
}

float Raytracer::GetAvgSamplesPerPixel() const
{
	float agvSpp = 0.0f;
	for (unsigned int p = 0; p < m_pixelCount; ++p)
	{
		agvSpp += m_sampleAccumulation[p];
	}
	agvSpp /= static_cast<float>(m_pixelCount);
	return agvSpp;
}

void Raytracer::SetRenderMode(RenderMode mode)
{
	if (m_renderMode != mode)
	{
		Clear();
		m_renderMode = mode;
	}
}

void Raytracer::SetRenderType(RenderType type)
{
	if (m_renderType != type)
	{
		Clear();
		m_renderType = type;
	}
}

void Raytracer::ToggleLightSampling()
{
	m_enableLightSampling = !m_enableLightSampling;
}

#pragma optimize("", off)
void Raytracer::TraceSampleAccumulation(const Ray& cameraRay)
{
	m_rayCount = 0;
	// camera x/y
	const float fov = 0.5135f;
	const Vec3 cx = Vec3(m_imgWidth * fov / m_imgHeight, 0, 0);
#if USE_OPTIMIZED_VEC 
	const Vec3 cy = normalize_vector(cross_product(cx, cameraRay.direction)) * fov;
#else
	const Vec3 cy = (cx.cross(cameraRay.direction)).normalize() * fov;
#endif

// #pragma omp parallel for schedule(dynamic, 1) 
	for (int rays = 0; rays < m_raysPerFrame; ++rays)
	{
		const unsigned int x = random::XorShift(0, m_imgWidth);
		const unsigned int y = random::XorShift(0, m_imgHeight);

		const unsigned int index = (m_imgHeight - y - 1) * m_imgWidth + x;
		++m_sampleAccumulation[index];

		for (unsigned int sy = 0; sy < 2u; sy++)
		{
			for (unsigned int sx = 0; sx < 2u; sx++)
			{
				const double r1 = 2 * random::XorShift(0.0f, 1.0f);
				const double dx = r1 < 1
					? sqrt(r1) - 1
					: 1 - sqrt(2 - r1);

				const double r2 = 2 * random::XorShift(0.0f, 1.0f);
				const double dy = r2 < 1
					? sqrt(r2) - 1
					: 1 - sqrt(2 - r2);

				Vec3 d =
					cx * (((sx + 0.5 + dx) * 0.5 + x) / m_imgWidth - 0.5) +
					cy * (((sy + 0.5 + dy) * 0.5 + y) / m_imgHeight - 0.5) +
					cameraRay.direction;

#if USE_OPTIMIZED_VEC 
				Ray ray = Ray(cameraRay.origin + d * 130, normalize_vector(d));
#else
				Ray ray = Ray(cameraRay.origin + d * 130, d.normalize());
#endif
				Vec3 r;
				switch (m_renderType)
				{
				case RenderType::AmbientOcclusion: r = TraceAO(ray); break;
				case RenderType::ShadowRays: r = TraceShadowRay(ray); break;
				case RenderType::Diffuse:
				default: r = TraceRay(ray, 0); break;
				}

#if USE_OPTIMIZED_VEC 
				m_colorData[index] = m_colorData[index] + Vec3(MathUtils::Clamp(r.get_x()) * 0.25, MathUtils::Clamp(r.get_y()) * 0.25, MathUtils::Clamp(r.get_z()) * 0.25);
#else
				m_colorData[index] = m_colorData[index] + Vec3(MathUtils::Clamp(r.x) * 0.25, MathUtils::Clamp(r.y) * 0.25, MathUtils::Clamp(r.z) * 0.25);
#endif
			}
		}

		const float oneBySampleCount = 1.0f / m_sampleAccumulation[index];

		// convert BGR to RGB
#if USE_OPTIMIZED_VEC 
		m_rgbData[index * 4u + 2u] = MathUtils::Linear2sRGB(m_colorData[index].get_x() * oneBySampleCount, 2.2f);
		m_rgbData[index * 4u + 1u] = MathUtils::Linear2sRGB(m_colorData[index].get_y() * oneBySampleCount, 2.2f);
		m_rgbData[index * 4u + 0u] = MathUtils::Linear2sRGB(m_colorData[index].get_z() * oneBySampleCount, 2.2f);
#else
		m_rgbData[index * 4u + 2u] = MathUtils::Linear2sRGB(m_colorData[index].x * oneBySampleCount, 2.2f);
		m_rgbData[index * 4u + 1u] = MathUtils::Linear2sRGB(m_colorData[index].y * oneBySampleCount, 2.2f);
		m_rgbData[index * 4u + 0u] = MathUtils::Linear2sRGB(m_colorData[index].z * oneBySampleCount, 2.2f);
#endif
		}
	}


Vec3 Raytracer::TraceAO(const Ray& r)
{
	double t = 0; // distance to intersection
	int id = 0; // id of intersected object
	if (!scene.Intersect(r, t, id))
	{
		// if miss, return black
		return Vec3(1, 1, 1);
	}

	// the hit object
	const Sphere& obj = scene.g_spheres[id];
	const Vec3 x = r.origin + r.direction * t;

#if USE_OPTIMIZED_VEC 
	const Vec3 n = normalize_vector(x - obj.GetPosition());
	const Vec3 nl = dot_product(n, r.direction) < 0 ? n : n * -1;
#else
	const Vec3 n = (x - obj.GetPosition()).normalize();
	const Vec3 nl = n.dot(r.direction) < 0 ? n : n * -1;
#endif

	float accum = 0.0f;
	for (int i = 0; i < m_aoSamples; ++i)
	{
		Vec3 rHemisphereVec = RandomInUnitHemisphere(n);

		double intersectTime = 0;
		int secId = 0;
		Ray secondRay(x, rHemisphereVec);
		if (scene.Intersect(secondRay, intersectTime, secId) && intersectTime <= m_aoRayLength)
		{
#if USE_OPTIMIZED_VEC 
			float nDotL = dot_product(n, rHemisphereVec);
#else
			float nDotL = rHemisphereVec.dot(n);
#endif
			accum += nDotL;
		}
	}

	float col = accum * m_invAoSamples;

	return Vec3(1.0f, 1.0f, 1.0f) - Vec3(col, col, col);
}

Vec3 Raytracer::TraceShadowRay(const Ray& r)
{
	double t = 0; // distance to intersection
	int id = 0; // id of intersected object
	if (!scene.Intersect(r, t, id))
	{
		// if miss, return black
		return Vec3(0, 0, 0);
	}

	// the hit object
	const Sphere& obj = scene.g_spheres[id];
	const Vec3 x = r.origin + r.direction * t;

#if USE_OPTIMIZED_VEC 
	const Vec3 n = normalize_vector(x - obj.GetPosition());
	const Vec3 nl = dot_product(n, r.direction) < 0 ? n : n * -1;
#else
	const Vec3 n = (x - obj.GetPosition()).normalize();
	const Vec3 nl = n.dot(r.direction) < 0 ? n : n * -1;
#endif#

	Vec3 f = obj.GetColor();

	// ideal DIFFUSE reflection
	if (obj.reflectionType == ReflectionType::DIFFUSE)
	{
		double r1 = 2 * PI * random::XorShift(0.0f, 1.0f);
		double r2 = random::XorShift(0.0f, 1.0f);
		double r2s = sqrt(r2);

		Vec3 w = nl;
#if USE_OPTIMIZED_VEC 
		Vec3 u = normalize_vector(cross_product(fabs(w.get_x()) > .1 ? Vec3(0, 1, 0) : Vec3(1, 0, 0), w));
		Vec3 v = cross_product(w, u);
		Vec3 d = normalize_vector(u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2));
#else
		Vec3 u = ((fabs(w.x) > .1 ? Vec3(0, 1, 0) : Vec3(1, 0, 0)).cross(w)).normalize();
		Vec3 v = w.cross(u);
		Vec3 d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).normalize();
#endif
		// loop over any lights
		Vec3 e(0, 0, 0);
		for (int i = 0; i < scene.GetSphereCount(); i++)
		{
			const Sphere& s = scene.g_spheres[i];

			// skip non-lights
			if (!s.isLight)
			{
				continue;
			}

			Vec3 sw = s.GetPosition() - x;

#if USE_OPTIMIZED_VEC 
			Vec3 su = normalize_vector(cross_product(fabs(sw.get_x()) > .1 ? Vec3(0, 1, 0) : Vec3(1, 0, 0), sw));
			Vec3 sv = cross_product(sw, su);
			double cos_a_max = sqrt(1 - s.radius * s.radius / dot_product((x - s.GetPosition()), (x - s.GetPosition())));
#else
			Vec3 su = ((fabs(sw.x) > .1 ? Vec3(0, 1, 0) : Vec3(1, 0, 0)).cross(sw)).normalize();
			Vec3 sv = sw.cross(su);
			double cos_a_max = sqrt(1 - s.radius * s.radius / (x - s.GetPosition()).dot(x - s.GetPosition()));
#endif

			double eps1 = random::XorShift(0.0f, 1.0f), eps2 = random::XorShift(0.0f, 1.0f);
			double cos_a = 1 - eps1 + eps1 * cos_a_max;
			double sin_a = sqrt(1 - cos_a * cos_a);
			double phi = 2 * PI * eps2;
			Vec3 l = su * cos(phi) * sin_a + sv * sin(phi) * sin_a + sw * cos_a;

#if USE_OPTIMIZED_VEC 
			l = normalize_vector(l);
#else
			l = l.normalize();
#endif

			// shadow ray
			if (scene.Intersect(Ray(x, l), t, id) && id == i)
			{
				double omega = 2 * PI * (1 - cos_a_max);

				// 1/pi for BRDF
				e = e + s.GetEmission() /* * l.dot(nl) * omega *  ONE_OVER_PI */;
			}
			else
			{
				e = Vec3(0, 0, 0);
			}
			}

		return e;
		}
	}

Vec3 Raytracer::TraceRay(const Ray& r, int depth, int E)
{
	++m_rayCount;

	double t = 0; // distance to intersection
	int id = 0; // id of intersected object
	if (!scene.Intersect(r, t, id))
	{
		// if miss, return black
		return Vec3(0, 0, 0);
	}

	// the hit object
	const Sphere& obj = scene.g_spheres[id];
	const Vec3 x = r.origin + r.direction * t;

#if USE_OPTIMIZED_VEC 
	const Vec3 n = normalize_vector(x - obj.GetPosition());
	const Vec3 nl = dot_product(n, r.direction) < 0 ? n : n * -1;
#else
	const Vec3 n = (x - obj.GetPosition()).normalize();
	const Vec3 nl = n.dot(r.direction) < 0 ? n : n * -1;
#endif#

	Vec3 f = obj.GetColor();

	// max refl
#if USE_OPTIMIZED_VEC 
	const double p = f.get_x() > f.get_y() && f.get_x() > f.get_z() ? f.get_x() : f.get_y() > f.get_z() ? f.get_y() : f.get_z();
#else
	const double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z;
#endif

	if (++depth > 5 || !p)
	{
		if (random::XorShift(0.0f, 1.0f) < p)
		{
			f = f * (1 / p);
		}
		else
		{
			return obj.GetEmission() * E;
		}
	}

	if (depth > 5)
	{
		return obj.GetEmission() * E;
	}

	// ideal DIFFUSE reflection
	if (obj.reflectionType == ReflectionType::DIFFUSE)
	{
		double r1 = 2 * PI * random::XorShift(0.0f, 1.0f);
		double r2 = random::XorShift(0.0f, 1.0f);
		double r2s = sqrt(r2);

		Vec3 w = nl;

#if USE_OPTIMIZED_VEC
		Vec3 u = normalize_vector(cross_product(fabs(w.get_x()) > .1 ? Vec3(0, 1, 0) : Vec3(1, 0, 0), w));
		Vec3 v = cross_product(w, u);
		Vec3 d = normalize_vector(u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2));
#else
		Vec3 u = ((fabs(w.x) > .1 ? Vec3(0, 1, 0) : Vec3(1, 0, 0)).cross(w)).normalize();
		Vec3 v = w.cross(u);
		Vec3 d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).normalize();
#endif

		// loop over any lights
		Vec3 e(0, 0, 0);
		for (int i = 0; i < scene.GetSphereCount(); i++)
		{
			const Sphere& s = scene.g_spheres[i];

			// skip non-lights
			if (!s.isLight)
			{
				continue;
			}

			Vec3 sw = s.GetPosition() - x;

#if USE_OPTIMIZED_VEC 
			Vec3 su = normalize_vector(cross_product(fabs(sw.get_x()) > .1 ? Vec3(0, 1, 0) : Vec3(1, 0, 0), sw));
			Vec3 sv = cross_product(sw, su);
			double cos_a_max = sqrt(1 - s.radius * s.radius / dot_product((x - s.GetPosition()), (x - s.GetPosition())));
#else
			Vec3 su = ((fabs(sw.x) > .1 ? Vec3(0, 1, 0) : Vec3(1, 0, 0)).cross(sw)).normalize();
			Vec3 sv = sw.cross(su);
			double cos_a_max = sqrt(1 - s.radius * s.radius / (x - s.GetPosition()).dot(x - s.GetPosition()));
#endif

			double eps1 = random::XorShift(0.0f, 1.0f), eps2 = random::XorShift(0.0f, 1.0f);
			double cos_a = 1 - eps1 + eps1 * cos_a_max;
			double sin_a = sqrt(1 - cos_a * cos_a);
			double phi = 2 * PI * eps2;
			Vec3 l = su * cos(phi) * sin_a + sv * sin(phi) * sin_a + sw * cos_a;
#if USE_OPTIMIZED_VEC 
			l = normalize_vector(l);
#else
			l = l.normalize();
#endif

			// shadow ray
			++m_rayCount;
			if (scene.Intersect(Ray(x, l), t, id) && id == i)
			{
				double omega = 2 * PI * (1 - cos_a_max);

				// 1/pi for BRDF
#if USE_OPTIMIZED_VEC 
				e = e + f * (s.GetEmission() * dot_product(l, nl) * omega) * ONE_OVER_PI;
#else
				e = e + f.mult(s.GetEmission() * l.dot(nl) * omega) * ONE_OVER_PI;
#endif
			}
		}
#if USE_OPTIMIZED_VEC 
		return obj.GetEmission() * E + e + f * (TraceRay(Ray(x, d), depth, 0));
#else
		return obj.GetEmission() * E + e + f.mult(TraceRay(Ray(x, d), depth, 0));
#endif
	}
	// Ideal SPECULAR reflection
	else if (obj.reflectionType == ReflectionType::SPECULAR)
	{
#if USE_OPTIMIZED_VEC 
		return obj.GetEmission() + f * (TraceRay(Ray(x, r.direction - n * 2 * dot_product(n, r.direction)), depth));
#else
		return obj.GetEmission() + f.mult(TraceRay(Ray(x, r.direction - n * 2 * n.dot(r.direction)), depth));
#endif
	}

	// Ideal dielectric REFRACTION
#if USE_OPTIMIZED_VEC 
	Ray reflRay(x, r.direction - n * 2 * dot_product(n, r.direction));
	// Ray from outside going in?
	bool into = dot_product(n, nl) > 0;
	double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = dot_product(r.direction, nl), cos2t;
	// Total internal reflection
	if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)
	{
		return obj.GetEmission() + f * (TraceRay(reflRay, depth));
	}
#else
	Ray reflRay(x, r.direction - n * 2 * n.dot(r.direction));
	// Ray from outside going in?
	bool into = n.dot(nl) > 0;
	double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.direction.dot(nl), cos2t;
	// Total internal reflection
	if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)
	{
		return obj.GetEmission() + f.mult(TraceRay(reflRay, depth));
	}
#endif

#if USE_OPTIMIZED_VEC 
	Vec3 tdir = normalize_vector(r.direction * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t))));
	double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : dot_product(tdir, n));
#else
	Vec3 tdir = (r.direction * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normalize();
	double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : tdir.dot(n));
#endif
	double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);

	// Russian roulette
#if USE_OPTIMIZED_VEC 
	return obj.GetEmission() + f * (depth > 2 ? (random::XorShift(0.0f, 1.0f) < P ?
		TraceRay(reflRay, depth) * RP
		: TraceRay(Ray(x, tdir), depth) * TP)
		: TraceRay(reflRay, depth) * Re + TraceRay(Ray(x, tdir), depth) * Tr);

#else
	return obj.GetEmission() + f.mult(depth > 2 ? (random::XorShift(0.0f, 1.0f) < P ?
		TraceRay(reflRay, depth) * RP
		: TraceRay(Ray(x, tdir), depth) * TP)
		: TraceRay(reflRay, depth) * Re + TraceRay(Ray(x, tdir), depth) * Tr);
#endif
	}
