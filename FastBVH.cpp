#include "precomp.h"
#include "FastBVH.h"

//struct BVHFlatNode {
//	BBox bbox;
//	uint start_offset, n_prims, right_offset;
//};

struct BVHBuildEntry {
	// if non-zero then this is the index of the parent
	uint parent;
	// the range of objects in the object list covered by this node
	uint start, end;
};

void FastBVH::printStats() {
	std::cout << "[Core] " << this->name << " with " << node_count << " nodes for " << (int)scene->objects.size()
		<< " primitives was build in " << this->build_time << "ms. Taking up "
		<< float(node_count * sizeof(BVHFlatNode)) / (1024.f*1024.f) << "MB.\n";
}

FastBVH::FastBVH(Scene* scene) {
	this->scene = scene;
	this->name = "Fast-BVH";

	build();
}

void FastBVH::build() {
	Tmpl8::timer build_timer;
	build_timer.reset();

	BVHBuildEntry todo[128];
	uint stackptr = 0;
	const uint Parent = 0xfffffffc;
	const uint Untouched = 0xffffffff;
	const uint TouchedTwice = 0xfffffffd;

	// Push the root
	todo[stackptr].start = 0;
	todo[stackptr].end = (uint)scene->objects.size();
	todo[stackptr].parent = Parent;
	stackptr++;

	BVHFlatNode node;
	vector<BVHFlatNode> buildnodes;
	buildnodes.reserve(scene->objects.size() * 2);

	while (stackptr > 0) {
		// pop the next item off of the stack
		BVHBuildEntry &bnode(todo[--stackptr]);
		uint start = bnode.start;
		uint end = bnode.end;
		uint nPrims = end - start;

		node_count++;
		node.start_offset = start;
		node.n_prims = nPrims;
		node.right_offset = Untouched;

		// calculate the bounding box for this node
		BBox bb(scene->objects[start]->getBBox());
		BBox bc(scene->objects[start]->getCentroid());
		for (uint p = start + 1; p < end; ++p) {
			bb.expandToInclude(scene->objects[p]->getBBox());
			bc.expandToInclude(scene->objects[p]->getCentroid());
		}
		node.bbox = bb;

		// if the number of primitives at this point is less than the leaf
		// size, then this will become a leaf. (Signified by rightOffset == 0)
		if (nPrims <= prims_in_leaf) {
			node.right_offset = 0;
			leaf_count++;
		}

		buildnodes.push_back(node);

		// cild touches parent...
		// special case: don't do this for the root.
		if (bnode.parent != Parent) {
			buildnodes[bnode.parent].right_offset--;

			// when this is the second touch, this is the right child.
			// the right child sets up the offset for the flat tree.
			if (buildnodes[bnode.parent].right_offset == TouchedTwice) {
				buildnodes[bnode.parent].right_offset = node_count - 1 - bnode.parent;
			}
		}

		// if this is a leaf, no need to subdivide.
		if (node.right_offset == 0)
			continue;

		uint split_dim;
		float split_coord;

		// set the split dimensions
		split_dim = bc.maxDimension();

		// split on the center of the longest axis
		split_coord = .5f * (bc.min[split_dim] + bc.max[split_dim]);

		// partition the list of objects on this split
		uint mid = start;
		for (uint i = start; i<end; ++i) {
			if (scene->objects[i]->getCentroid()[split_dim] < split_coord) {
				swap(scene->objects[i], scene->objects[mid]);
				++mid;
			}
		}

		// if we get a bad split, just choose the center...
		if (mid == start || mid == end) {
			mid = start + (end - start) / 2;
		}


		// push right child
		todo[stackptr].start = mid;
		todo[stackptr].end = end;
		todo[stackptr].parent = node_count - 1;
		stackptr++;

		// push left child
		todo[stackptr].start = start;
		todo[stackptr].end = mid;
		todo[stackptr].parent = node_count - 1;
		stackptr++;
	}

	// copy the temp node data to a flat array
	flat_tree = new BVHFlatNode[node_count];
	for (uint n = 0; n<node_count; ++n)
		flat_tree[n] = buildnodes[n];

	this->build_time = build_timer.elapsed();
}

FastBVH::~FastBVH() {
	delete[] flat_tree;
}

uint FastBVH::getNodeCount() {
	return node_count;
}

Intersection FastBVH::getIntersection(const Ray& ray, bool fast_occlusion) {
	Intersection result;

	float bbhits[4];
	std::int32_t closer, other;

	// working set
	BVHTraversal todo[64];
	std::int32_t stackptr = 0;

	// pushing the root node
	todo[stackptr].i = 0;
	todo[stackptr].mint = -INFINITY;

	while (stackptr >= 0) {
		// pop the next working node
		int ni = todo[stackptr].i;
		float near_t = todo[stackptr].mint;
		stackptr--;
		const BVHFlatNode &node(flat_tree[ni]);

		// if it's further than the closest found intersection, continue
		if (near_t > result.distance)
			continue;

		// node is leaf, intersect it
		if (node.right_offset == 0) {
			for (uint o = 0; o<node.n_prims; o++) {
				Intersection current;

				Primitive* obj = scene->objects[node.start_offset + o];
				current = obj->intersect(ray);
#if _DEBUG
				isecs++;
#endif
				if (current.didIntersect) {

					// any intersection is good enough?
					if (fast_occlusion) {
						return current;
					}

					// keep the closest intersection only
					if (current.distance < result.distance) {
						result = current;
					}
				}
			}

		}
		else { // not a leaf

			bool hitc0 = flat_tree[ni + 1].bbox.intersect(ray, bbhits, bbhits + 1);
			bool hitc1 = flat_tree[ni + node.right_offset].bbox.intersect(ray, bbhits + 2, bbhits + 3);

			// both nodes
			if (hitc0 && hitc1) {

				// assume left is closer
				closer = ni + 1;
				other = ni + node.right_offset;

				// if the right is closer, swap it with left
				if (bbhits[2] < bbhits[0]) {
					swap(bbhits[0], bbhits[2]);
					swap(bbhits[1], bbhits[3]);
					swap(closer, other);
				}

				// push the farther first
				todo[++stackptr] = BVHTraversal(other, bbhits[2]);
				todo[++stackptr] = BVHTraversal(closer, bbhits[0]);
			}
			else if (hitc0) {
				todo[++stackptr] = BVHTraversal(ni + 1, bbhits[0]);
			}
			else if (hitc1) {
				todo[++stackptr] = BVHTraversal(ni + node.right_offset, bbhits[2]);
			}

		}
	}

	return result;
}
