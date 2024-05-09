#pragma once
#include "Primitive.h"

class Sphere : public Primitive {
public:
	vec3 center;
	float r, r2;

	Sphere(const vec3& center, const float radius, size_t matId)
		: center(center), r(radius), r2(radius*radius) {
		this->matId = matId;
		this->type = SPHERE;
	}

	// intersect the given ray with this sphere
	Intersection intersect(const Ray& ray) {
		float t0, t1;

		vec3 L = ray.origin - center;
		float a = dot(ray.direction, ray.direction);
		float b = 2 * dot(ray.direction, L);
		float c = dot(L, L) - r2;

		if (!solveQuadratic(a, b, c, t0, t1))
			return Intersection(); // no hit

		if (t0 > t1) std::swap(t0, t1);

		if (t0 < 0.0f) {
			t0 = t1; // using t1, because t0 is negative
			if (t0 < 0.0f || std::isnan(t0)) return Intersection(); // both t0 and t1 are negative
		}

		vec3 intersection = ray.origin + (ray.direction * t0);
		vec3 normal = getNormal(intersection);

		// flip normal if this is a refractive ray
		if (dot(normal, ray.direction) > 0) {
			normal = -normal;
		}

		return Intersection(t0, this, intersection, normal, ray, ray.matId, matId);
	}

	// returns normal at the given point, assume point is on surface
	vec3 getNormal(const vec3& point) const {
		return normalize(point - center);
	}

	// returns bounding box for this sphere
	BBox getBBox() const {
		vec3 min = center - vec3(r);
		vec3 max = center + vec3(r);
		return BBox(min, max);
	}

	// returns centroid (center) of this sphere
	vec3 getCentroid() const {
		return center;
	}

	int getPrimitiveCount() const {
		return 1;
	}
};
