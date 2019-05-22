
#pragma once

#include <vector>

#include "Sphere.h"
#include "BVH.h"


struct Ray;

class Scene
{
public:
	Scene();

	bool Intersect(const Ray& r, double& t, int& id);
	bool IntersectBVH(const Ray& r, double& t, int& id);

	unsigned int GetSphereCount() const { return static_cast<unsigned int>(g_spheres.size()); }

	std::vector<Sphere> g_spheres;
	BVH m_bvh;
};