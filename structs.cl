#define NULL_INDEX 0xffffffffu

typedef struct General {
	int triangle_count;
	int lights_count;

	bool cast_shadows;
	float ambient_light;
	int super_sampling;
	int max_reflections;
	float3 background_color;
	float epsilon;
} General;

typedef struct Camera {
	float3 position;
	float3 direction;
	float3 right;
	float3 up;
	int width, height;
	float aspect_ratio;
} Camera;

typedef struct PointLight {
	float3 position;
	float3 color;
	float intensity;
} PointLight;

typedef struct Material {
	float3 diffuse;
	float refractive_index;
	float shininess;
	float reflectivity;
} Material;

typedef struct Triangle {
	float3 vertices[3];
	float3 vert_normals[3];
	float3 normal;
	int mat_id;
} Triangle;

typedef struct {
	float3 bounds[2];
	unsigned int primitive;
	unsigned int skipIndex;
} BVHNode;

typedef struct Ray {
	float3 origin;
	float3 direction;
	float3 inverse_direction;
	bool fast_occlusion;
	int reflections_remaining;
} Ray;

typedef struct Intersection {
	bool did_intersect;
	float distance;
	float3 position;
	float3 normal;
	float3 bary;
	int mat_id;
} Intersection;