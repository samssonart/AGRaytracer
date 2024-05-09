#pragma once
#include "BBox.h"
#include "Intersection.h"
#include "Ray.h"
#include "Material.h"

// enum that has types of available primitives
enum PrimitiveType { SPHERE, TRIANGLE, MESH };

class Primitive {
protected:
	// generic quadratic equation solver
	bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1) const
	{
		float discr = b * b - 4 * a * c;
		if (discr < 0) return false;
		else if (discr == 0.0f) {
			x0 = x1 = -0.5f * b / a;
		}
		else {
			float q = (b > 0.0f) ?
				-0.5f * (b + std::sqrt(discr)) :
				-0.5f * (b - std::sqrt(discr));
			x0 = q / a;
			x1 = c / q;
		}

		return true;
	}

public:
	PrimitiveType type; // enum with type
	size_t matId;
	Primitive* parent = NULL;

	// intersect given ray with this primitive object
	virtual Intersection intersect(const Ray& ray) = 0;

	// return primitive object's normal at the given point
	virtual vec3 getNormal(const vec3& point) const = 0;

	// return bounding box for this primitive object
	virtual BBox getBBox() const = 0;

	// return centroid of this primitive object
	virtual vec3 getCentroid() const = 0;

	// return primitive count
	virtual int getPrimitiveCount() const = 0;
};
