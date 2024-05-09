#pragma once
#include "Material.h"

class Flat : public Material {
private:
	vec3 color;

public:
	Flat(vec3 color) {
		this->color = color;
	}

	vec3 getColor(vec3 point) {
		return this->color;
	}

	vec3 getColor() {
		return this->color;
	}
};