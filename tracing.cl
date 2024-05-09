float3 cast_ray(
	__global const General* general,
	__global const PointLight* lights,
	__global const Triangle* triangles,
	__global const BVHNode* bvh,
	__global const Material* materials,
	const Ray ray) {

	Intersection intersection = getClosestIntersection(bvh, triangles, ray);

	if (intersection.did_intersect) {
		return performLighting(general, bvh, triangles, lights, materials, ray, intersection);
	}
	else {
		return general->background_color;
	}
}
