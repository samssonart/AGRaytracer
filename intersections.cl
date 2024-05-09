bool intersectTriangle(const Triangle tr, const Ray ray, Intersection* intersection, const int tr_index) {

	float3 e1 = tr.vertices[1] - tr.vertices[0];
	float3 e2 = tr.vertices[2] - tr.vertices[0];
	float3 s1 = cross(ray.direction, e2);

	const float divisor = dot(s1, e1);
	if (divisor == 0.f) return false;

	// native_recip is significantly faster than just division
	const float invDivisor = 1.0f / divisor; //todo native_recip

	// First barycentric coord
	const float3 d = ray.origin - tr.vertices[0];
	const float b1 = dot(d, s1) * invDivisor;
	if (b1 < 0.f || b1 > 1.f) return false;

	// Second barycentric coord
	const float3 s2 = cross(d, e1);
	const float b2 = dot(ray.direction, s2) * invDivisor;
	if (b2 < 0.f || (b1 + b2) > 1.f) return false;

	// Distance to intersection point
	const float t = dot(e2, s2) * invDivisor;
	if (t < 0.0f || t > intersection->distance) return false;

	intersection->did_intersect = true;
	intersection->distance = t;
	intersection->position = ray.origin + (ray.direction * t);
	intersection->normal = tr.vert_normals[0] + b1 * (tr.vert_normals[1] - tr.vert_normals[0]) + b2 * (tr.vert_normals[2] - tr.vert_normals[0]);
	intersection->bary = (float3)(b1, b2, 0.0f);
	intersection->mat_id = tr.mat_id;
	
	return true;
}

bool intersectBBox(const float3 b_min, const float3 b_max, const Ray ray) {

	float t1 = (b_min.x - ray.origin.x) * ray.inverse_direction.x;
	float t2 = (b_max.x - ray.origin.x) * ray.inverse_direction.x;
	float t3 = (b_min.y - ray.origin.y) * ray.inverse_direction.y;
	float t4 = (b_max.y - ray.origin.y) * ray.inverse_direction.y;
	float t5 = (b_min.z - ray.origin.z) * ray.inverse_direction.z;
	float t6 = (b_max.z - ray.origin.z) * ray.inverse_direction.z;
	
	float tmin = max(max(min(t1, t2), min(t3, t4)), min(t5, t6));
	float tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6));
	
	// if tmax < 0, ray (line) is intersecting aabb, but whole aabb is behing us
	if (tmax < 0) {
		//t = tmax;
		return false;
	}
	
	// if tmin > tmax, ray doesn't intersect aabb
	if (tmin > tmax) {
		//t = tmax;
		return false;
	}
	
	//t = tmin;
	return true;
}

bool intersectBVH(__global const BVHNode* nodes, __global const Triangle* triangles,
	const Ray ray, Intersection* intersection) {

	bool wasFound = false;
	uint currentNode = 0;
	const uint stopNode = nodes[0].skipIndex;

	while (currentNode < stopNode) {
		const __global BVHNode *node = &nodes[currentNode];

		if (intersectBBox(node->bounds[0], node->bounds[1], ray)) {

			const unsigned int triIndex = node->primitive;

			if (triIndex != NULL_INDEX)
				intersectTriangle(triangles[triIndex], ray, intersection, triIndex);

			if (intersection->did_intersect) {
				if (ray.fast_occlusion) {
					return true;
				}

				wasFound = true;
			}

			currentNode++;
		}
		else {
			currentNode = node->skipIndex;
		}
	}

	return wasFound;

}

Intersection getClosestIntersection(__global const BVHNode* bvh, __global const Triangle* triangles, const Ray ray) {

	Intersection intersection = { false, INFINITY };

	intersectBVH(bvh, triangles, ray, &intersection);

	// naive intersection, for debugging purposes, add general struct to the method to use
	//for (int i=0; i < general->triangle_count; i++) {
	//	intersectTriangle(triangles[i], ray, &intersection, i);
	//}

	return intersection;
}