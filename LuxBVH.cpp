#include "precomp.h"

bool bvh_ltf_x(BVHAccelTreeNode *n, float v) {
	return n->bbox.max.x + n->bbox.min.x < v;
}

bool bvh_ltf_y(BVHAccelTreeNode *n, float v) {
	return n->bbox.max.y + n->bbox.min.y < v;
}

bool bvh_ltf_z(BVHAccelTreeNode *n, float v) {
	return n->bbox.max.z + n->bbox.min.z < v;
}

bool(*const bvh_ltf[3])(BVHAccelTreeNode *n, float v) = { bvh_ltf_x, bvh_ltf_y, bvh_ltf_z };


BVHAccelTreeNode* LuxBVH::BuildHierarchy(std::vector<BVHAccelTreeNode *> &list, unsigned int begin, unsigned int end, unsigned int axis) {
	unsigned int splitAxis = axis;
	float splitValue;

	totalNodes += 1;
	if (end - begin == 1) // Only a single item in list so return it
		return list[begin];

	BVHAccelTreeNode *parent = new BVHAccelTreeNode();
	parent->primitive = 0xffffffffu;
	parent->leftChild = NULL;
	parent->rightSibling = NULL;

	std::vector<unsigned int> splits;
	splits.reserve(treeType + 1);
	splits.push_back(begin);
	splits.push_back(end);
	for (unsigned int i = 2; i <= treeType; i *= 2) { // Calculate splits, according to tree type and do partition
		for (unsigned int j = 0, offset = 0; j + offset < i && splits.size() > j + 1; j += 2) {
			if (splits[j + 1] - splits[j] < 2) {
				j--;
				offset++;
				continue; // Less than two elements: no need to split
			}

			FindBestSplit(list, splits[j], splits[j + 1], &splitValue, &splitAxis);

			std::vector<BVHAccelTreeNode *>::iterator it =
				partition(list.begin() + splits[j], list.begin() + splits[j + 1], bind2nd(ptr_fun(bvh_ltf[splitAxis]), splitValue));
			unsigned int middle = std::distance(list.begin(), it);
			middle = std::max(splits[j] + 1, std::min(splits[j + 1] - 1, middle)); // Make sure coincidental BBs are still split
			splits.insert(splits.begin() + j + 1, middle);
		}
	}

	BVHAccelTreeNode *child, *lastChild;
	// Left Child
	child = BuildHierarchy(list, splits[0], splits[1], splitAxis);
	parent->leftChild = child;
	parent->bbox = child->bbox;
	lastChild = child;

	// Add remaining children
	for (unsigned int i = 1; i < splits.size() - 1; i++) {
		child = BuildHierarchy(list, splits[i], splits[i + 1], splitAxis);
		lastChild->rightSibling = child;
		parent->bbox.expandToInclude(child->bbox);
		lastChild = child;
	}

	return parent;
}