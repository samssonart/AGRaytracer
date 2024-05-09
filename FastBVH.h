#pragma once
#include "Accelerator.h"

struct BVHFlatNode {
	BBox bbox;
	uint start_offset, n_prims, right_offset;
};

class FastBVH : public Accelerator {
private:
	uint node_count = 0, leaf_count = 0, prims_in_leaf = 4;
	BVHFlatNode* flat_tree;

	void build();

public:
	uint FastBVH::getNodeCount();
	FastBVH(Scene* scene);
	~FastBVH();
	Intersection FastBVH::getIntersection(const Ray& ray, bool fast_occlusion);
	void FastBVH::printStats();
	BVHFlatNode* getTree() { return flat_tree; }

};