#pragma once
#include "glm/glm.hpp"
#include "Material.h"
using namespace glm;

class Ray {
public:
	vec3 origin; // Ray Origin
	vec3 direction; // Ray Direction
	size_t matId;
	int reflections_remaining;

	Ray(const vec3& o, const vec3& d, int reflections_remaining, size_t matId)
		: origin(o), direction(normalize(d)), matId(matId), reflections_remaining(reflections_remaining) {
	}

	Ray() : origin(vec3(0.0f)), direction(vec3(0.0f)), matId(0), reflections_remaining(-1) {}
};