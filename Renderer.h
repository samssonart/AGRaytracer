#pragma once
#include "glm\glm.hpp"

typedef unsigned long Pixel;

class Renderer {

protected:

	// Converts glm::vec3 to Tmpl8::Pixel (unsigned long)
	Tmpl8::Pixel toPixel(glm::vec3 color) const {
		return ((int)(color.x * 255) << 16) + ((int)(color.y * 255) << 8) + (int)(color.z * 255);
	}

public:

	std::string name, device, platform;

	std::vector<Accelerator*> accelerators;
	std::size_t current_accelerator_id;
	std::string accelerators_ui_string;

	virtual void render(Tmpl8::Surface* screen) = 0;

};