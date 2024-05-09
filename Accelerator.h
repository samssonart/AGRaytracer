#pragma once
#include <vector>
#include <algorithm>
#include <string>

#include "glm\glm.hpp"
#include "BBox.h"

class Scene;
class Intersection;
class Ray;

// node for storing state information during traversal.
struct BVHTraversal {
	uint i; // node
	float mint; // Minimum hit time for this node.
	BVHTraversal() { }
	BVHTraversal(int _i, float _mint) : i(_i), mint(_mint) { }
};

class Accelerator {

public:

	std::string name;	// name of the accelerator
	Scene* scene;		// pointer to scene on which acc. structure is built
	float build_time;	// build time of the structure
#if _DEBUG
	uint isecs = 0;
#endif
	virtual Intersection getIntersection(const Ray& ray, bool fast_occlusion) = 0;
	virtual uint getNodeCount() = 0;
	virtual void printStats() = 0;
};