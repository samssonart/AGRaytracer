#pragma once
#include "Material.h"

class Air : public Material {
public:
	Air() {}

	vec3 getColor(vec3 point) {
		return vec3(0.0f);
	}

	vec3 getColor() {
		return vec3(0.0f);
	}

	float getRefractiveIndex() {
		return 1.0f; // air refraction index
	}
};