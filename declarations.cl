float3 cast_ray(
	__global const General* general,
	__global const PointLight* lights,
	__global const Triangle* triangles,
	__global const BVHNode* bvh,
	__global const Material* materials,
	const Ray ray);

bool intersectTriangle(const Triangle tr, const Ray ray, Intersection* intersection, const int tr_index);

bool intersectBBox(const float3 b_min, const float3 b_max, const Ray ray);

bool intersectBVH(__global const BVHNode* nodes, __global const Triangle* triangles,
	const Ray ray, Intersection* intersection);

Intersection getClosestIntersection(__global const BVHNode* bvh, __global const Triangle* triangles, const Ray ray);

float3 performLighting(
	__global const General* general,
	__global const BVHNode* bvh,
	__global const Triangle* triangles,
	__global const PointLight* lights,
	__global const Material* materials,
	const Ray ray,
	const Intersection intersection);

float3 get_reflectance(float3 normal, float3 incident, float n1, float n2);

float3 reflect_vector(float3 direction, float3 normal);

