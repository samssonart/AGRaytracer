#pragma once
#include "glm/glm.hpp"
#include "Material.h"
class Primitive;

class Intersection {
public:
	float distance; // intersection point distance from ray origin
	Primitive* object = NULL; // intersected object
	vec3 point; // intersection point
	bool didIntersect; // intersection confirmation
	vec3 normal; // intersection normal
	Ray ray; // ray that caused intersection
	size_t startMaterial;
	size_t endMaterial;

	// Constructor for hits
	Intersection(float distance, Primitive* object, vec3 point, vec3 normal, Ray ray, size_t startMaterial, size_t endMaterial)
		: distance(distance), object(object), point(point), didIntersect(true),
		normal(normal), ray(ray), startMaterial(startMaterial), endMaterial(endMaterial) { }

	// Empty constructor, for no hits
	Intersection()
		: distance(INFINITY), object(NULL), point(vec3(0.0f)), didIntersect(false) { }
};