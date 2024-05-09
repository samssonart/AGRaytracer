#pragma once
#include "glm/glm.hpp"

class Light {
public:
	vec3 location;
	float intensity;

	Light() {}

	Light(vec3 loc, float inten) {
		location = loc;
		intensity = inten;
	}
};