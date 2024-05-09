#pragma once
#include "Accelerator.h"

class NoAcceleration : public Accelerator {

public:

	NoAcceleration(Scene* scene) {

		this->name = "Naive";
		
		this->scene = scene;
		this->build_time = 0.0f;
	}
		
	Intersection getIntersection(const Ray& ray, bool occlusion) {
		float current_distance = INFINITY;
		Intersection result;

		for (size_t i = 0; i < scene->objects.size(); i++) {
			Intersection current_intersection = scene->objects[i]->intersect(ray);

			if (current_intersection.didIntersect) {
				if (current_intersection.distance < current_distance) {
					current_distance = current_intersection.distance;
					result = current_intersection;
				}
			}
		}

		return result;
	}

	uint getNodeCount() {
		return (uint)scene->primitive_count;
	}

	void printStats() {}




};