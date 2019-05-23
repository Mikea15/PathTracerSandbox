#pragma once

#include <memory>
#include <vector>

#include "Sphere.h"

struct AABB
{
	Vec3 min;
	Vec3 max;

	bool Intersect(const Ray& r, float tmin, float tmax) const
	{
		for (int a = 0; a < 3; ++a) 
		{
			float invD = 1.0f / r.direction[a];

			float t0 = (min[a] - r.origin[a]) * invD;
			float t1 = (max[a] - r.origin[a]) * invD;
			
			if (invD < 0.0f) 
			{
				float temp = t1;
				t1 = t0;
				t0 = temp;
			}

			tmin = t0 > tmin ? t0 : tmin;
			tmax = t1 < tmax ? t1 : tmax;

			if (tmax <= tmin)
				return false;
		}
		return true;
	}
};

class BVHNode
{
public:
	BVHNode* left;
	BVHNode* right;

	int sphereIndex;
	AABB boundingBox;
};

class BVH
{
public:
	BVH() {}

	void Add(const Sphere& sphere, int index) 
	{
		BVHNode node;
		Vec3 rad(sphere.radius, sphere.radius, sphere.radius);
		AABB boundingBox;
		boundingBox.min = sphere.GetPosition() - rad;
		boundingBox.max = sphere.GetPosition() + rad;

		node.boundingBox = boundingBox;
		node.sphereIndex = index;

		leafs.push_back(node);
	}

	bool Intersect(const std::vector<Sphere>& world, const Ray& r, double& t, int& id) 
	{
		double inf = t = 1e20;
		int bvhId = id;
		
		for (const BVHNode& n : leafs)
		{
			if (n.boundingBox.Intersect(r, 1e-5, 1e10))
			{
				bvhId = n.sphereIndex;
				double d = world[bvhId].intersect(r);
				if (d && d < t)
				{
					t = d;
					id = bvhId;
				}
			}
		}

		return t < inf;
	}

private:
	std::unique_ptr<BVHNode> root;

	std::vector<BVHNode> leafs;
};