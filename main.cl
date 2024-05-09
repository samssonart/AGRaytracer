__kernel void raytrace(
	__global const General* general,
	__global const Camera* camera,
	__global const PointLight* lights,
	__global const Triangle* triangles,
	__global const BVHNode* bvh,
	__global int* buffer,
	__global const Material* materials) {

	// kernel task id
	int idx = get_global_id(0);
	int idy = get_global_id(1);
	int id = idx + idy * camera->width;

	// general check
	if (id >= camera->width * camera->height) return;

	// primary ray generation
	float normalized_x = ((float)idx / camera->width) - 0.5f;
	float normalized_y = ((float)idy / camera->height) - 0.5f;

	float3 image_point = camera->direction +
		normalized_x * camera->right * camera->aspect_ratio +
		normalized_y * camera->up;
	image_point = normalize(image_point);

	Ray ray = { camera->position, image_point, 1.f / image_point, false, general->max_reflections };

	// casting primary rays, does lighting etc
	float3 result = cast_ray(general, lights, triangles, bvh, materials, ray);

	// convert to 0-255 RGB, result float3 is clamped 0-1
	int r = (int)(result.x * 255.0f);
	int g = (int)(result.y * 255.0f);
	int b = (int)(result.z * 255.0f);

	// write the screen buffer
	buffer[id] = (r << 16) + (g << 8) + b;
}

__kernel void generate_primary_rays(
	__global const General* general,
	__global const Camera* camera,
	__global Ray* rays) {

	// kernel task id
	int idx = get_global_id(0);
	int idy = get_global_id(1);
	int id = idx + idy * camera->width;

	// general bounds check
	if (id >= camera->width * camera->height) return;

	// primary ray generation
	float normalized_x = ((float)idx / camera->width) - 0.5f;
	float normalized_y = ((float)idy / camera->height) - 0.5f;

	float3 image_point = camera->direction +
		normalized_x * camera->right * camera->aspect_ratio +
		normalized_y * camera->up;
	image_point = normalize(image_point);

	Ray ray = { camera->position, image_point, 1.f / image_point, false, general->max_reflections, id };

	// write back
	rays[id] = ray;
}

__kernel void compute_intersections(
	__global const Ray* rays,
	__global const Triangle* triangles,
	__global const BVHNode* bvh,
	__global const Camera* camera,
	__global Intersection* intersections) {

	// kernel task id
	int idx = get_global_id(0);
	int idy = get_global_id(1);
	int id = idx + idy * camera->width;

	// general bounds check
	if (id >= camera->width * camera->height) return;
	intersections[id] = getClosestIntersection(bvh, triangles, rays[id]);

}

__kernel void compute_shading(
	__global const Intersection* intersections,
	__global const Material* materials,
	__global const Ray* rays,
	__global const PointLight* lights,
	__global const General* general,
	__global const Camera* camera,
	__global const BVHNode* bvh,
	__global const Triangle* triangles,
	__global int* buffer) {

	// kernel task id
	int idx = get_global_id(0);
	int idy = get_global_id(1);
	int id = idx + idy * camera->width;

	// general bounds check
	if (id >= camera->width * camera->height) return;

	// casting primary rays, does lighting etc
	float3 result;
	if (intersections[id].did_intersect) {
		result = performLighting(general, bvh, triangles, lights, materials, rays[id], intersections[id]);
	}
	else {
		result = general->background_color;
	}

	// convert to 0-255 RGB, result float3 is clamped 0-1
	int r = (int)(result.x * 255.0f);
	int g = (int)(result.y * 255.0f);
	int b = (int)(result.z * 255.0f);

	// write the screen buffer
	buffer[id] = (r << 16) + (g << 8) + b;
}
