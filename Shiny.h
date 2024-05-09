#pragma once
#include "Material.h"

class Shiny : public Material {
private:
	vec3 color;
	float shininess;
	float reflectivity;

public:
	Shiny(vec3 color, float shininess, float reflectivity) {
		this->color = color;
		this->shininess = shininess;
		this->reflectivity = reflectivity;
	}

	vec3 getColor(vec3 point) {
		return this->color;
	}

	vec3 getColor() {
		return this->color;
	}

	float getShininess() {
		return shininess;
	}

	float getReflectivity() {
		return reflectivity;
	}
};