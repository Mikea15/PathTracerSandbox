
#pragma once

#include "Sphere.h"
#include <vector>

struct Ray;

class Scene
{
public:
	Scene();

	bool Intersect(const Ray& r, double& t, int& id);

	unsigned int GetSphereCount() const { return static_cast<unsigned int>(g_spheres.size()); }

	std::vector<Sphere> g_spheres;
};