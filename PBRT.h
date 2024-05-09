#pragma once
#include "Accelerator.h"
#include "ppl.h"

// BVH Builder ported and adapted from Physically Based Rendering v2 book by Matt Pharr & Greg Humphreys

// Originally it has three split methods:
// SPLIT_MIDDLE, which first computes the midpoint of the primitives’ centroids along the longest axis
// SPLIT_EQUAL_COUNTS, which partition primitives into equally - sized subsets
// SPLIT_SAH, which uses SAH, choosing a partition of the primitives along the chosen axis that gives a minimal cost estimate

// I've implemented one more split method for PBRT builder - SPLIT_AAC, using paper
// "Efficient BVH Construction via Approximate Agglomerative Clustering" by Gu et al.
// Paper: http://www.cs.cmu.edu/~ygu1/paper/HPG13/HPG13.pdf
// Slides: http://highperformancegraphics.org/wp-content/uploads/Gu-AAC.pdf

#define MORTON_CODE_END 0
#define MORTON_CODE_START 29

enum SplitMethod { SPLIT_MIDDLE, SPLIT_EQUAL_COUNTS, SPLIT_SAH, SPLIT_AAC_HQ, SPLIT_AAC_FAST };

struct BVHPrimitiveInfo {
	BVHPrimitiveInfo() { }
	BVHPrimitiveInfo(int pn, const BBox &b)
		: primitiveNumber(pn), bounds(b) {
		centroid = .5f * b.min + .5f * b.max;
	}
	int primitiveNumber;
	vec3 centroid;
	BBox bounds;
};

struct BVHBuildNode {
	// BVHBuildNode Public Methods
	BVHBuildNode() { children[0] = children[1] = NULL; }
	void InitLeaf(uint first, uint n, const BBox &b) {
		firstPrimOffset = first;
		nPrimitives = n;
		bounds = b;
	}
	void InitInterior(uint axis, BVHBuildNode *c0, BVHBuildNode *c1) {
		children[0] = c0;
		children[1] = c1;
		bounds.expandToInclude(c0->bounds);
		bounds.expandToInclude(c1->bounds);
		splitAxis = axis;
		nPrimitives = 0;
	}
	BBox bounds;
	BVHBuildNode *children[2];
	uint splitAxis, firstPrimOffset, nPrimitives;
};

struct CompareToMid {
	CompareToMid(int d, float m) { dim = d; mid = m; }
	int dim;
	float mid;
	bool operator()(const BVHPrimitiveInfo &a) const {
		return a.centroid[dim] < mid;
	}
};

struct ComparePoints {
	ComparePoints(int d) { dim = d; }
	int dim;
	bool operator()(const BVHPrimitiveInfo &a,
		const BVHPrimitiveInfo &b) const {
		return a.centroid[dim] < b.centroid[dim];
	}
};

struct CompareToBucket {
	CompareToBucket(int split, int num, int d, const BBox &b)
		: centroidBounds(b)
	{
		splitBucket = split; nBuckets = num; dim = d;
	}
	bool operator()(const BVHPrimitiveInfo &p) const {
		int b = nBuckets * ((p.centroid[dim] - centroidBounds.min[dim]) /
			(centroidBounds.max[dim] - centroidBounds.min[dim]));
		if (b == nBuckets) b = nBuckets - 1;
		assert(b >= 0 && b < nBuckets);
		return b <= splitBucket;
	}

	int splitBucket, nBuckets, dim;
	const BBox &centroidBounds;
};

struct LinearBVHNode {
	BBox bounds;
	union {
		unsigned int primitivesOffset;    // leaf
		unsigned int secondChildOffset;   // interior
	};

	unsigned short nPrimitives;  // 0 -> interior node
	unsigned short axis;         // interior node: xyz
};


class PBRT : public Accelerator {

public:

	uint maxPrimsInNode = 4;
	LinearBVHNode *nodes;
	uint totalNodes = 0;
	vector<Primitive*> primitives;

	uint AAC_DELTA;
	float AAC_EPSILON;
	float AAC_ALPHA;
	
	SplitMethod splitMethod;

	void build() {

		Tmpl8::timer build_timer;
		build_timer.reset();

		if (primitives.size() == 0) {
			nodes = NULL;
			return;
		}

		vector<BVHPrimitiveInfo> buildData;
		buildData.reserve(primitives.size());

		for (uint i = 0; i < primitives.size(); ++i) {
			BBox bbox = primitives[i]->getBBox();
			buildData.push_back(BVHPrimitiveInfo(i, bbox));
		}

		vector<Primitive*> orderedPrims;
		orderedPrims.reserve(primitives.size());
				
		if (splitMethod == SPLIT_AAC_FAST || splitMethod == SPLIT_AAC_HQ) {
			buildAAC(buildData, 0, (uint)primitives.size(), orderedPrims);
		}
		else {
			BVHBuildNode *root = recursiveBuild(buildData, 0,
				(uint)primitives.size(), orderedPrims);

			primitives.swap(orderedPrims);

			nodes = new LinearBVHNode[totalNodes];

			uint offset = 0;
			flattenBVHTree(root, &offset);
			assert(offset == totalNodes);
		}

		build_time = build_timer.elapsed();
	}

	bool IntersectP(const BBox &bounds, const Ray &ray,
		const vec3 &invDir, const uint dirIsNeg[3]) {
		// Check for ray intersection against $x$ and $y$ slabs
		float tmin = (bounds[dirIsNeg[0]].x - ray.origin.x) * invDir.x;
		float tmax = (bounds[1 - dirIsNeg[0]].x - ray.origin.x) * invDir.x;
		float tymin = (bounds[dirIsNeg[1]].y - ray.origin.y) * invDir.y;
		float tymax = (bounds[1 - dirIsNeg[1]].y - ray.origin.y) * invDir.y;
		if ((tmin > tymax) || (tymin > tmax))
			return false;
		if (tymin > tmin) tmin = tymin;
		if (tymax < tmax) tmax = tymax;

		// Check for ray intersection against $z$ slab
		float tzmin = (bounds[dirIsNeg[2]].z - ray.origin.z) * invDir.z;
		float tzmax = (bounds[1 - dirIsNeg[2]].z - ray.origin.z) * invDir.z;
		if ((tmin > tzmax) || (tzmin > tmax))
			return false;
		if (tzmin > tmin)
			tmin = tzmin;
		if (tzmax < tmax)
			tmax = tzmax;
		return (tmin < INFINITY) && (tmax > 0); // return (tmin < ray.maxt) && (tmax > ray.mint);
	}

	BVHBuildNode* recursiveBuild(vector<BVHPrimitiveInfo> &buildData, uint start,
		uint end, vector<Primitive*> &orderedPrims) {
		assert(start != end);

		totalNodes++;

		BVHBuildNode* node = new BVHBuildNode();

		// Compute bounds of all primitives in BVH node
		BBox bbox;
		for (uint i = start; i < end; ++i)
			bbox.expandToInclude(buildData[i].bounds);

		uint nPrimitives = end - start;

		if (nPrimitives == 1) {
			// Create leaf _BVHBuildNode_
			uint firstPrimOffset = (uint)orderedPrims.size();
			for (uint i = start; i < end; ++i) {
				uint primNum = buildData[i].primitiveNumber;
				orderedPrims.push_back(primitives[primNum]);
			}
			node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
		}
		else {
			// Compute bound of primitive centroids, choose split dimension _dim_
			BBox centroidBounds;
			for (uint i = start; i < end; ++i)
				centroidBounds.expandToInclude(buildData[i].centroid);
			int dim = centroidBounds.maxDimension();

			// Partition primitives into two sets and build children
			uint mid = (start + end) / 2;
			if (centroidBounds.max[dim] == centroidBounds.min[dim]) {
				// If nPrimitives is no greater than maxPrimsInNode,
				// then all the nodes can be stored in a compact bvh node.
				if (nPrimitives <= maxPrimsInNode) {
					// Create leaf _BVHBuildNode_
					uint firstPrimOffset = (uint)orderedPrims.size();
					for (uint i = start; i < end; ++i) {
						uint primNum = buildData[i].primitiveNumber;
						orderedPrims.push_back(primitives[primNum]);
					}
					node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
					return node;
				}
				else {
					// else if nPrimitives is greater than maxPrimsInNode, we
					// need to split it further to guarantee each node contains
					// no more than maxPrimsInNode primitives.
					node->InitInterior(dim,
						recursiveBuild(buildData, start, mid, orderedPrims),
						recursiveBuild(buildData, mid, end, orderedPrims));
					return node;
				}
			}

			// Partition primitives based on _splitMethod_
			switch (splitMethod) {
			case SPLIT_MIDDLE: {
				// Partition primitives through node's midpoint
				float pmid = .5f * (centroidBounds.min[dim] + centroidBounds.max[dim]);
				BVHPrimitiveInfo *midPtr = std::partition(&buildData[start],
					&buildData[end - 1] + 1,
					CompareToMid(dim, pmid));
				mid = midPtr - &buildData[0];
				if (mid != start && mid != end)
					// for lots of prims with large overlapping bounding boxes, this
					// may fail to partition; in that case don't break and fall through
					// to SPLIT_EQUAL_COUNTS
					break;
			}
			case SPLIT_EQUAL_COUNTS: {
				// Partition primitives into equally-sized subsets
				mid = (start + end) / 2;
				std::nth_element(&buildData[start], &buildData[mid],
					&buildData[end - 1] + 1, ComparePoints(dim));
				break;
			}
			case SPLIT_SAH: default: {
				// Partition primitives using approximate SAH
				if (nPrimitives <= 4) {
					// Partition primitives into equally-sized subsets
					mid = (start + end) / 2;
					std::nth_element(&buildData[start], &buildData[mid],
						&buildData[end - 1] + 1, ComparePoints(dim));
				}
				else {
					// Allocate _BucketInfo_ for SAH partition buckets
					const int nBuckets = 12;
					struct BucketInfo {
						BucketInfo() { count = 0; }
						int count;
						BBox bounds;
					};
					BucketInfo buckets[nBuckets];

					// Initialize _BucketInfo_ for SAH partition buckets
					for (uint i = start; i < end; ++i) {
						int b = nBuckets *
							((buildData[i].centroid[dim] - centroidBounds.min[dim]) /
							(centroidBounds.max[dim] - centroidBounds.min[dim]));
						if (b == nBuckets) b = nBuckets - 1;
						assert(b >= 0 && b < nBuckets);
						buckets[b].count++;
						buckets[b].bounds.expandToInclude(buildData[i].bounds);
					}

					// Compute costs for splitting after each bucket
					float cost[nBuckets - 1];
					for (int i = 0; i < nBuckets - 1; ++i) {
						BBox b0, b1;
						int count0 = 0, count1 = 0;
						for (int j = 0; j <= i; ++j) {
							b0.expandToInclude(buckets[j].bounds);
							count0 += buckets[j].count;
						}
						for (int j = i + 1; j < nBuckets; ++j) {
							b1.expandToInclude(buckets[j].bounds);
							count1 += buckets[j].count;
						}
						cost[i] = .125f + (count0*b0.surfaceArea() + count1*b1.surfaceArea()) /
							bbox.surfaceArea();
					}

					// Find bucket to split at that minimizes SAH metric
					float minCost = cost[0];
					uint minCostSplit = 0;
					for (int i = 1; i < nBuckets - 1; ++i) {
						if (cost[i] < minCost) {
							minCost = cost[i];
							minCostSplit = i;
						}
					}

					// Either create leaf or split primitives at selected SAH bucket
					if (nPrimitives > maxPrimsInNode ||
						minCost < nPrimitives) {
						BVHPrimitiveInfo *pmid = std::partition(&buildData[start],
							&buildData[end - 1] + 1,
							CompareToBucket(minCostSplit, nBuckets, dim, centroidBounds));
						mid = pmid - &buildData[0];
					}

					else {
						// Create leaf _BVHBuildNode_
						uint firstPrimOffset = (uint)orderedPrims.size();
						for (uint i = start; i < end; ++i) {
							uint primNum = buildData[i].primitiveNumber;
							orderedPrims.push_back(primitives[primNum]);
						}
						node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
						return node;
					}
				}
				break;
			}
			}
			node->InitInterior(dim,
				recursiveBuild(buildData, start, mid, orderedPrims),
				recursiveBuild(buildData, mid, end, orderedPrims));
		}
		return node;
	}

	uint flattenBVHTree(BVHBuildNode *node, uint *offset) {
		LinearBVHNode *linearNode = &nodes[*offset];
		linearNode->bounds = node->bounds;
		uint myOffset = (*offset)++;
		if (node->nPrimitives > 0) {
			assert(!node->children[0] && !node->children[1]);
			linearNode->primitivesOffset = node->firstPrimOffset;
			linearNode->nPrimitives = node->nPrimitives;
		}
		else {
			// Creater interior flattened BVH node
			linearNode->axis = node->splitAxis;
			linearNode->nPrimitives = 0;
			flattenBVHTree(node->children[0], offset);
			linearNode->secondChildOffset = flattenBVHTree(node->children[1],
				offset);
		}

		delete node;
		return myOffset;
	}

	PBRT(Scene* scene, SplitMethod method, string name = "PBRT v2 BVH") {

		this->scene = scene;
		this->name = name;
		this->splitMethod = method;

		if (this->splitMethod == SPLIT_AAC_FAST) {
			AAC_DELTA = 4;
			AAC_EPSILON = 0.2f;
			AAC_ALPHA = 0.5f - AAC_EPSILON;
		}
		else if (this->splitMethod == SPLIT_AAC_HQ) {
			AAC_DELTA = 20;
			AAC_EPSILON = 0.1f;
			AAC_ALPHA = 0.5f - AAC_EPSILON;
		}

		this->primitives = scene->objects;

		build();
	}

	~PBRT() {
		delete[] nodes;
	}

	Intersection getIntersection(const Ray& ray, bool fast_occlusion) {

		if (!nodes) return Intersection();

		Intersection result;

		vec3 invDir(1.f / ray.direction.x, 1.f / ray.direction.y, 1.f / ray.direction.z);
		uint dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
		// Follow ray through BVH nodes to find primitive intersections
		uint todoOffset = 0, nodeNum = 0;
		uint todo[64];
		while (true) {
			const LinearBVHNode *node = &nodes[nodeNum];
			// Check ray against BVH node
			//float a, b;
			if (
				//node->bounds.intersect(ray, &a, &b) //28fps
				//node->bounds.intersect(ray) //23fps
				IntersectP(node->bounds, ray, invDir, dirIsNeg) //33fps
				) {
				if (node->nPrimitives > 0) {
					// Intersect ray with primitives in leaf BVH node
					for (uint i = 0; i < node->nPrimitives; ++i) {
/*
						printf("i : %u \n", i);
						printf("BBmin : %f %f %f \n", node->bounds.min.x, node->bounds.min.y, node->bounds.min.z);
						printf("BBmax : %f %f %f \n", node->bounds.max.x, node->bounds.max.y, node->bounds.max.z);
						printf("Offset : %u %u %u \n\n", node->nPrimitives, node->primitivesOffset, node->primitivesOffset + i);*/

						Intersection current;
						Primitive* obj = primitives[node->primitivesOffset + i];
						current = obj->intersect(ray);
#if _DEBUG
						isecs++;
#endif

						if (current.didIntersect) {

							if (fast_occlusion)
								return current;

							if (current.distance < result.distance) {
								result = current;
							}

						}
					}
					if (todoOffset == 0) break;
					nodeNum = todo[--todoOffset];
				}
				else {
					// Put far BVH node on _todo_ stack, advance to near node
					if (dirIsNeg[node->axis]) {
						todo[todoOffset++] = nodeNum + 1;
						nodeNum = node->secondChildOffset;
					}
					else {
						todo[todoOffset++] = node->secondChildOffset;
						nodeNum = nodeNum + 1;
					}
				}
			}
			else {
				if (todoOffset == 0) break;
				nodeNum = todo[--todoOffset];
			}
		}

		return result;
	}

	uint getNodeCount() {

		return totalNodes;
	}

	void printStats() {
		std::cout << "[Core] " << this->name << " with " << totalNodes << " nodes for " << (int)primitives.size()
			<< " primitives was build in " << this->build_time << "ms. Taking up "
			<< float(totalNodes * sizeof(LinearBVHNode)) / (1024.f*1024.f) << "MB.\n";
	}

	// Expands a 10-bit integer into 30 bits
	// by inserting 2 zeros after each bit.
	std::uint32_t expandBits(std::uint32_t v) {
		v = (v * 0x00010001u) & 0xFF0000FFu;
		v = (v * 0x00000101u) & 0x0F00F00Fu;
		v = (v * 0x00000011u) & 0xC30C30C3u;
		v = (v * 0x00000005u) & 0x49249249u;
		return v;
	}

	// Calculates a 30-bit Morton code for the given 3D point located within the unit cube [0,1]
	std::uint32_t morton3D(float x, float y, float z) {
		x = std::min(std::max(x * 1024.0f, 0.0f), 1023.0f);
		y = std::min(std::max(y * 1024.0f, 0.0f), 1023.0f);
		z = std::min(std::max(z * 1024.0f, 0.0f), 1023.0f);

		std::uint32_t xx = expandBits((std::uint32_t)x);
		std::uint32_t yy = expandBits((std::uint32_t)y);
		std::uint32_t zz = expandBits((std::uint32_t)z);
		return ((zz << 2) | (yy << 1) | xx);
	}

	float AAC_C() {
		return (0.5f * powf(AAC_DELTA, 0.5f + AAC_EPSILON));
		
	}

	std::uint32_t AAC_F(std::uint32_t x) {
		return (std::uint32_t)(ceil(AAC_C() * powf(x, AAC_ALPHA)));
	}


	std::uint32_t makePartition(std::pair<std::uint32_t, std::uint32_t> *sortedMC,
		std::uint32_t start, std::uint32_t end, int partitionbit) {

		std::uint32_t curFind = (1 << partitionbit);

		if (((sortedMC[start].first & curFind) == (sortedMC[end - 1].first & curFind))
			|| (partitionbit < MORTON_CODE_END)) {
			if (partitionbit < MORTON_CODE_END)
				cout << "makePartition spliting in half. SUPPOSED TO BE RARE";
			return start + (end - start) / 2;
		}

		std::uint32_t lower = start;
		std::uint32_t upper = end;

		//std::uint32_t curFind = (1 << partitionbit);
		while (lower < upper) {
			std::uint32_t mid = lower + (upper - lower) / 2;
			if ((sortedMC[mid].first & curFind) == 0) {
				lower = mid + 1;
			}
			else {
				upper = mid;
			}
		}

		return lower;
	}

	std::uint32_t findBestMatch(std::vector<BVHBuildNode*> &clusters, std::uint32_t i) {
		float closestDist = INFINITY;
		std::uint32_t idx = i;
		for (std::uint32_t j = 0; j < clusters.size(); ++j) {
			if (i == j) continue;

			BBox combined;
			combined.expandToInclude(clusters[i]->bounds);
			combined.expandToInclude(clusters[j]->bounds);

			float d = combined.surfaceArea();
			if (d < closestDist) {
				closestDist = d;
				idx = j;
			}
		}

		return idx;
	}

	std::vector<BVHBuildNode*> combineCluster(std::vector<BVHBuildNode*> &clusters,
		std::uint32_t n, int dim) {

		std::vector<std::uint32_t> closest(clusters.size(), 0);

		for (std::uint32_t i = 0; i < clusters.size(); ++i) {
			closest[i] = findBestMatch(clusters, i);
		}

		while (clusters.size() > n) {
			float bestDist = INFINITY;
			std::uint32_t leftIdx = 0;
			std::uint32_t rightIdx = 0;

			for (std::uint32_t i = 0; i < clusters.size(); ++i) {


				BBox combined;
				combined.expandToInclude(clusters[i]->bounds);
				combined.expandToInclude(clusters[closest[i]]->bounds);
				
				float d = combined.surfaceArea();
				if (d < bestDist) {
					bestDist = d;
					leftIdx = i;
					rightIdx = closest[i];
				}
			}

			totalNodes++;
			BVHBuildNode* node = new BVHBuildNode();

			node->InitInterior(dim,
				clusters[leftIdx],
				clusters[rightIdx]);
			clusters[leftIdx] = node;
			clusters[rightIdx] = clusters.back();
			closest[rightIdx] = closest.back();
			clusters.pop_back();
			closest.pop_back();

			closest[leftIdx] = findBestMatch(clusters, leftIdx);

			for (std::uint32_t i = 0; i < clusters.size(); ++i) {
				if (closest[i] == leftIdx || closest[i] == rightIdx)
					closest[i] = findBestMatch(clusters, i);
				else if (closest[i] == closest.size()) {
					closest[i] = rightIdx;
				}
			}
		}


		return clusters;
	}

	void recursiveBuildAAC(vector<BVHPrimitiveInfo> &buildData,
		std::pair<std::uint32_t, std::uint32_t> *mortonCodes,
		std::uint32_t start, std::uint32_t end, int partitionBit,
		std::vector<BVHBuildNode*> *clusterData) {

		std::vector<BVHBuildNode*> clusters;
		if (end - start == 0) {
			return;
		}

		int dim = partitionBit % 3; // need 0-2 for dimension

		if (end - start < AAC_DELTA) {
			std::vector<BVHBuildNode*> clusters;
			totalNodes += (end - start);

			for (std::uint32_t i = start; i < end; ++i) {
				// Create leaf _BVHBuildNode_
				BVHBuildNode* node = new BVHBuildNode();
				std::uint32_t primIdx = mortonCodes[i].second;

				node->InitLeaf(primIdx, 1, buildData[primIdx].bounds); // deal with firstPrimOffset later with DFS
				clusters.push_back(node);
			}

			*clusterData = combineCluster(clusters, AAC_F(AAC_DELTA), dim);
			return;
		}


		std::uint32_t splitIdx = makePartition(mortonCodes, start, end, partitionBit);

		int newPartionBit = partitionBit - 1;
		std::vector<BVHBuildNode*> leftC;
		std::vector<BVHBuildNode*> rightC;

		std::uint32_t rightTotalnodes = 0;

		recursiveBuildAAC(buildData, mortonCodes, start, splitIdx, newPartionBit, &leftC);
		recursiveBuildAAC(buildData, mortonCodes, splitIdx, end, newPartionBit, &rightC);

		leftC.insert(leftC.end(), rightC.begin(), rightC.end());
		*clusterData = combineCluster(leftC, AAC_F(end - start), dim);
	}

	std::uint32_t bvhDfs(BVHBuildNode* node, vector<BVHPrimitiveInfo> &buildData,
		vector<Primitive*> &orderedPrims, std::uint32_t *offset) {

		LinearBVHNode *linearNode = &nodes[*offset];
		linearNode->bounds = node->bounds;
		std::uint32_t myOffset = (*offset)++;

		if (node->nPrimitives > 0) {
			assert(!node->children[0] && !node->children[1]);
			std::uint32_t firstPrimOffset = orderedPrims.size();
			std::uint32_t primNum = buildData[node->firstPrimOffset].primitiveNumber;
			orderedPrims.push_back(primitives[primNum]);
			node->firstPrimOffset = firstPrimOffset;
			linearNode->primitivesOffset = node->firstPrimOffset;
			linearNode->nPrimitives = node->nPrimitives;
		}
		else {
			linearNode->axis = node->splitAxis;
			linearNode->nPrimitives = 0;
			bvhDfs(node->children[0], buildData, orderedPrims, offset);
			linearNode->secondChildOffset = bvhDfs(node->children[1], buildData, orderedPrims, offset);
		}

		delete node;
		return myOffset;
	}

	void buildAAC(vector<BVHPrimitiveInfo> &buildData, std::uint32_t start,
		std::uint32_t end, vector<Primitive*> &orderedPrims) {

		std::vector<vec3> bboxCenters;
		std::pair<std::uint32_t, std::uint32_t> *mortonCodes = new std::pair<std::uint32_t, std::uint32_t>[end - start];

		BBox centerBBox = BBox();
		for (std::uint32_t i = start; i < end; ++i) {
			centerBBox.expandToInclude(buildData[i].centroid);
		}

		for (std::uint32_t i = start; i < end; ++i) {
			vec3 c = buildData[i].centroid;
			float newX = (c.x - centerBBox.min.x) / (centerBBox.max.x - centerBBox.min.x);
			float newY = (c.y - centerBBox.min.y) / (centerBBox.max.y - centerBBox.min.y);
			float newZ = (c.z - centerBBox.min.z) / (centerBBox.max.z - centerBBox.min.z);
			std::uint32_t mc = morton3D(newX, newY, newZ);

			mortonCodes[i - start] = std::make_pair(mc, i);
		}

		//concurrency::parallel_sort(mortonCodes, mortonCodes + (end - start)); // TODO is taking longer than std::sort somehow
		std::sort(mortonCodes, mortonCodes+(end-start));

		std::vector<BVHBuildNode*> clusters;
		recursiveBuildAAC(buildData, mortonCodes, start, end, MORTON_CODE_START, &clusters);
		BVHBuildNode* root = combineCluster(clusters, 1, 2)[0];

		delete[] mortonCodes;

		// Compute representation of depth-first traversal of BVH tree
		nodes = new LinearBVHNode[totalNodes];
		for (std::uint32_t i = 0; i < totalNodes; ++i)
			new (&nodes[i]) LinearBVHNode;
		std::uint32_t offset = 0;

		bvhDfs(root, buildData, orderedPrims, &offset);
		primitives.swap(orderedPrims);
	}
};