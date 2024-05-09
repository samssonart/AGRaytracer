#pragma once
#include "Primitive.h"

class Triangle : public Primitive {
public:
	vec3 a, b, c;
	vec3 a_norm, b_norm, c_norm, face_normal;

	Triangle(vec3 a, vec3 b, vec3 c, size_t matId)
		: a(a), b(b), c(c) {
		this->type = TRIANGLE;
		this->matId = matId;
		this->face_normal = normalize(cross((b - a), (c - a)));
		a_norm = b_norm = c_norm = face_normal;
	}

	Triangle(vec3 a, vec3 b, vec3 c, vec3 a_norm, vec3 b_norm, vec3 c_norm, size_t matId)
		: a(a), b(b), c(c), a_norm(a_norm), b_norm(b_norm), c_norm(c_norm) {
		this->type = TRIANGLE;
		this->matId = matId;
		this->face_normal = normalize(cross((b - a), (c - a)));
	}

	Intersection intersect(const Ray& ray) {
		float epsilon = 0.000001f;
		vec3 e1, e2; // edges
		vec3 P, Q, T;
		float det, inv_det, u, v;
		float t;

		// find vectors for two edges sharing V1
		e1 = b - a;
		e2 = c - a;

		// begin calculating determinant - also used to calculate u parameter
		P = cross(ray.direction, e2);
		// if determinant is near zero, ray lies in plane of triangle
		det = dot(e1, P);
		// NOT CULLING
		if (det > -epsilon && det < epsilon)
			return Intersection();

		inv_det = 1.f / det;

		// calculate distance from V1 to ray origin
		T = ray.origin - a;

		// calculate u parameter and test bound
		u = dot(T, P) * inv_det;
		// the intersection lies outside of the triangle
		if (u < 0.f || u > 1.f) 
			return Intersection();

		// prepare to test v parameter
		Q = cross(T, e1); 


		// calculate v parameter and test bound
		v = dot(ray.direction, Q) *  inv_det;
		// the intersection lies outside of the triangle
		if (v < 0.f || u + v  > 1.f) 
			return Intersection();

		t = dot(e2, Q) * inv_det;

		if (t > epsilon) { // ray intersection
			vec3 point = ray.origin + (ray.direction * t);

			//return Intersection(t, this, point, getNormal(point, u, v), ray, ray.material, material);
			vec3 interpolated_normal = getNormal(point, u, v);

			// corrects refractions, flips normals for NOT front facing 
			// TODO: do this only for refracted rays?
			if (dot(-ray.direction, interpolated_normal) <= 0.0f) {
				interpolated_normal = -interpolated_normal;
			}

			return Intersection(t, this, point, interpolated_normal, ray, ray.matId, matId);
		}

		// no hit, no win
		return Intersection();
	}

	vec3 getNormal(const vec3& point, float u, float v) {
		vec3 interpolated = a_norm + u * (b_norm - a_norm) + v * (c_norm - a_norm);
		return normalize(interpolated);
	}

	vec3 getNormal(const vec3& point) const {
		return face_normal;
	}

	BBox getBBox() const {
		return BBox(a, b, c);
	}

	vec3 getCentroid() const {
		return (a + b + c) / 3.0f;
	}

	int getPrimitiveCount() const {
		return 1;
	}

	float Triangle::getArea() const {
		vec3 e1 = b - a;
		vec3 e2 = c - a;
		vec3 n = cross(e1, e2);
		float a = 0.5f * std::fabs(length(n));
		return a;
	}
};
