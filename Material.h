#pragma once
#include "glm/glm.hpp"

#define NOT_SHINY -1
#define NOT_REFLECTIVE -1
#define NOT_REFRACTIVE -1
#define AIR_REFRACTIVE_INDEX 1

class Material {

public:

	virtual glm::vec3 getColor(glm::vec3 point) {
		return glm::vec3(0.0f);
	}

	virtual glm::vec3 getColor() {
		return glm::vec3(0.0f);
	}

	virtual float getShininess() {
		return NOT_SHINY;
	}

	virtual float getReflectivity() {
		return NOT_REFLECTIVE;
	}

	virtual float getRefractiveIndex() {
		return NOT_REFRACTIVE;
	}

};