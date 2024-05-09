#pragma once
#include "Accelerator.h"
#include "Intersection.h"

#include <iostream>
#include <functional>
#include <algorithm>
#include <limits>

using std::bind2nd;
using std::ptr_fun;

struct BVHAccelTreeNode {
	BBox bbox;
	unsigned int primitive;
	BVHAccelTreeNode *leftChild;
	BVHAccelTreeNode *rightSibling;
};

struct BVHAccelArrayNode {
	BBox bbox;
	unsigned int primitive;
	unsigned int skipIndex;
};




class LuxBVH : public Accelerator {

public:

	BVHAccelArrayNode* nodes;
	std::vector<Primitive*> primitives;
	unsigned int totalNodes = 0;
	unsigned int treeType = 16; // 4/8 is best on CPU, 16 best on GPU 
	int costSamples = 0, isectCost = 100, traversalCost = 10;
	float emptyBonus = 0.5f;

	LuxBVH(Scene* scene) {
		this->name = "Luxrays bvhaccel.cpp";
		this->primitives = scene->objects;

		Tmpl8::timer build_timer;
		build_timer.reset();

		std::vector<BVHAccelTreeNode*> bvList;
		for (unsigned int i = 0; i < scene->primitive_count; ++i) {
			BVHAccelTreeNode *ptr = new BVHAccelTreeNode();
			ptr->bbox = this->primitives[i]->getBBox();

			ptr->primitive = i;
			ptr->leftChild = NULL;
			ptr->rightSibling = NULL;
			bvList.push_back(ptr);
		}

		totalNodes = 0;
		BVHAccelTreeNode *rootNode = BuildHierarchy(bvList, 0, bvList.size(), 2);

		nodes = new BVHAccelArrayNode[totalNodes];
		BuildArray(rootNode, 0);
		FreeHierarchy(rootNode);

		build_time = build_timer.elapsed();
	}

	~LuxBVH() {
		delete nodes;
	}

	void FreeHierarchy(BVHAccelTreeNode *node) {
		if (node) {
			FreeHierarchy(node->leftChild);
			FreeHierarchy(node->rightSibling);

			delete node;
		}
	}

	BVHAccelTreeNode* LuxBVH::BuildHierarchy(std::vector<BVHAccelTreeNode *> &list, unsigned int begin, unsigned int end, unsigned int axis);

	void FindBestSplit(std::vector<BVHAccelTreeNode *> &list, unsigned int begin, unsigned int end, float *splitValue, unsigned int *bestAxis) {
		if (end - begin == 2) {
			// Trivial case with two elements
			*splitValue = (list[begin]->bbox.max[0] + list[begin]->bbox.min[0] +
				list[end - 1]->bbox.max[0] + list[end - 1]->bbox.min[0]) / 2;
			*bestAxis = 0;
		}
		else {
			// Calculate BBs mean center (times 2)
			vec3 mean2(0, 0, 0), var(0, 0, 0);
			for (unsigned int i = begin; i < end; i++)
				mean2 += list[i]->bbox.max + list[i]->bbox.min;
			mean2 /= static_cast<float>(end - begin);

			// Calculate variance
			for (unsigned int i = begin; i < end; i++) {
				vec3 v = list[i]->bbox.max + list[i]->bbox.min - mean2;
				v.x *= v.x;
				v.y *= v.y;
				v.z *= v.z;
				var += v;
			}
			// Select axis with more variance
			if (var.x > var.y && var.x > var.z)
				*bestAxis = 0;
			else if (var.y > var.z)
				*bestAxis = 1;
			else
				*bestAxis = 2;

			if (costSamples > 1) {
				BBox nodeBounds;
				for (unsigned int i = begin; i < end; i++)
					nodeBounds.expandToInclude(list[i]->bbox);

				vec3 d = nodeBounds.max - nodeBounds.min;
				const float invTotalSA = 1.f / nodeBounds.surfaceArea();

				// Sample cost for split at some points
				float increment = 2 * d[*bestAxis] / (costSamples + 1);
				float bestCost = INFINITY;
				for (float splitVal = 2 * nodeBounds.min[*bestAxis] + increment; splitVal < 2 * nodeBounds.max[*bestAxis]; splitVal += increment) {
					int nBelow = 0, nAbove = 0;
					BBox bbBelow, bbAbove;
					for (unsigned int j = begin; j < end; j++) {
						if ((list[j]->bbox.max[*bestAxis] + list[j]->bbox.min[*bestAxis]) < splitVal) {
							nBelow++;
							bbBelow.expandToInclude(list[j]->bbox);
						}
						else {
							nAbove++;
							bbAbove.expandToInclude(list[j]->bbox);
						}
					}
					const float pBelow = bbBelow.surfaceArea() * invTotalSA;
					const float pAbove = bbAbove.surfaceArea() * invTotalSA;
					float eb = (nAbove == 0 || nBelow == 0) ? emptyBonus : 0.f;
					float cost = traversalCost + isectCost * (1.f - eb) * (pBelow * nBelow + pAbove * nAbove);
					// Update best split if this is lowest cost so far
					if (cost < bestCost) {
						bestCost = cost;
						*splitValue = splitVal;
					}
				}
			}
			else {
				// Split in half around the mean center
				*splitValue = mean2[*bestAxis];
			}
		}
	}

	unsigned int BuildArray(BVHAccelTreeNode *node, unsigned int offset) {
		// Build array by recursively traversing the tree depth-first
		while (node) {
			BVHAccelArrayNode *p = &nodes[offset];

			p->bbox = node->bbox;
			p->primitive = node->primitive;
			offset = BuildArray(node->leftChild, offset + 1);
			p->skipIndex = offset;

			node = node->rightSibling;
		}

		return offset;
	}

	Intersection getIntersection(const Ray& ray, bool fast_occlusion) {

		if (!nodes) return Intersection();

		unsigned int currentNode = 0;
		unsigned int stopNode = nodes[0].skipIndex;
		Intersection result;

		while (currentNode < stopNode) {
			if (nodes[currentNode].bbox.intersect(ray)) {
				if (nodes[currentNode].primitive != 0xffffffffu) {
					Intersection current;
					Primitive* obj = primitives[nodes[currentNode].primitive];
					current = obj->intersect(ray);
					if (current.didIntersect) {

						if (fast_occlusion)
							return current;

						if (current.distance < result.distance) {
							result = current;
						}

					}
				}

				currentNode++;
			}
			else
				currentNode = nodes[currentNode].skipIndex;
		}

		return result;
	}

	unsigned int getNodeCount() {
		return totalNodes;
	}

	void printStats() {
		std::cout << "[Core] " << this->name << " with " << totalNodes << " nodes for " << (int)primitives.size()
			<< " primitives was build in " << this->build_time << "ms. Taking up "
			<< float(totalNodes * sizeof(BVHAccelArrayNode)) / (1024.f*1024.f) << "MB.\n";
	}
};