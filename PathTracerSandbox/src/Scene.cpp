
#include "Scene.h"

Scene::Scene()
{
	g_spheres.push_back(Sphere({ 1e5,	{ 1e5 + 1, 40.8, 81.6 },	{ 0, 0, 0 },		{ 0.75, 0.25, 0.25 },		ReflectionType::DIFFUSE,	false }));	// Left
	g_spheres.push_back(Sphere({ 1e5,	{ -1e5 + 99, 40.8, 81.6 },	{ 0, 0, 0 },		{ 0.25, 0.25, 0.75 },		ReflectionType::DIFFUSE,	false }));	// Right
	g_spheres.push_back(Sphere({ 1e5,	{ 50, 40.8, 1e5 },			{ 0, 0, 0 },		{ 0.75, 0.75, 0.75 },		ReflectionType::DIFFUSE,	false }));	// Back
	g_spheres.push_back(Sphere({ 1e5,	{ 50, 40.8, -1e5 + 170 },	{ 0, 0, 0 },		{ 0.00, 0.50, 0.00 },		ReflectionType::DIFFUSE,	false }));	// Front
	g_spheres.push_back(Sphere({ 1e5,	{ 50, 1e5, 81.6 },			{ 0, 0, 0 },		{ 0.75, 0.75, 0.75 },		ReflectionType::DIFFUSE,	false }));	// Bottom
	g_spheres.push_back(Sphere({ 1e5,	{ 50, -1e5 + 81.6, 81.6 },	{ 0, 0, 0 },		{ 0.75, 0.75, 0.75 },		ReflectionType::DIFFUSE,	false }));	// Top
	g_spheres.push_back(Sphere({ 16.5,	{ 27, 16.5, 47 },			{ 0, 0, 0 },		{ 0.999, 0.999, 0.999 },	ReflectionType::DIFFUSE,	false }));	// Mirror
	g_spheres.push_back(Sphere({ 16.5,	{ 73, 16.5, 78 },			{ 0, 0, 0 },		{ 0.999, 0.999, 0.999 },	ReflectionType::REFRACTIVE,	false }));	// Glass
	g_spheres.push_back(Sphere({ 1.5,	{ 50, 81.6 - 16.5, 81.6 },	{ 400, 400, 400 },	{ 0, 0, 0 },				ReflectionType::DIFFUSE,	true }));	// Light
}

bool Scene::Intersect(const Ray& r, double& t, int& id)
{
	double inf = t = 1e20;
	for (unsigned int i = 0u; i < 9; ++i)
	{
		double d = g_spheres[i].intersect(r);
		if (d && d < t)
		{
			t = d;
			id = i;
		}
	}

	return t < inf;
}

