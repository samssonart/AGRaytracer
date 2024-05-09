#pragma once
#include "Material.h"

class Glass : public Material {
private:
	float refractive_index;
	float shininess;

public:
	Glass(float refractiveIndex, float shininess) {
		this->refractive_index = refractiveIndex;
		this->shininess = shininess;
	}

	vec3 getColor(vec3 point) {
		return vec3(0.0f);
	}

	vec3 getColor() {
		return vec3(0.0f);
	}

	float getRefractiveIndex() {
		return refractive_index;
	}

	float getShininess() {
		return shininess;
	}
};