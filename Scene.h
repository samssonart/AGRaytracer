#pragma once
#include <vector>
#include "Accelerator.h"
#include "Primitive.h"
#include "Renderer.h"
#include "Mesh.h"
#include "Camera.h"
#include "Light.h"
#include "Air.h"

class Scene {

public:
	std::vector<Primitive*> objects;
	std::vector<Mesh*> meshes;
	std::vector<Light*> lights;
	std::vector<Renderer*> renderers;
	std::vector<Material*> materials;
	Camera camera;

	int primitive_count = 0;

	// run time changable params
	size_t current_renderer_id = 0;
	bool cast_shadows = true;
	float ambient_light_coeff = 0.0f;
	int supersampling = 1;
	int max_reflections = 3;
	float background_color[3] = { 0.0f, 0.75f, 1.0f };
	float epsilon = 0.01f;

	// fps
	float calculated_fps = 0.0f;

	Scene() {
		materials.push_back(new Air()); // starting medium for correct refractions
	}

};