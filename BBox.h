#pragma once
#include "Ray.h"
#include "glm\glm.hpp"
#include <algorithm>

class BBox {
public:
	vec3 min, max;
	BBox();
	BBox(const vec3& min, const vec3& max);
	BBox(const vec3& p);
	BBox(const vec3& a, const vec3& b, const vec3& c);

	bool intersect(const Ray& ray, float *tnear, float *tfar) const;
	bool intersect(const Ray& ray) const;
	void expandToInclude(const vec3& p);
	void expandToInclude(const BBox& b);
	short maxDimension() const;
	float surfaceArea() const;
	vec3 BBox::getCenter() const;
	float distance(const BBox& b);
	const vec3 &operator[](int i) const;
	vec3 &operator[](int i);
};