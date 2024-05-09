#include "precomp.h"
#include "BBox.h"
#define QUALITY_OVER_PERFORMANCE 1 // 1 uses non-SSE version, 0 uses SSE version with much lower precision

BBox::BBox()
	: min(vec3(INFINITY)), max(vec3(-INFINITY)) {}

BBox::BBox(const glm::vec3& min, const glm::vec3& max)
	: min(min), max(max) {
}

BBox::BBox(const glm::vec3& p)
	: min(p), max(p) {
}

BBox::BBox(const glm::vec3& a, const glm::vec3& b, const glm::vec3& c) {
	expandToInclude(a);
	expandToInclude(b);
	expandToInclude(c);
}

const vec3 &BBox::operator[](int i) const {
	assert(i == 0 || i == 1);
	if (i == 0) {
		return min;
	}
	else {
		return max;
	}
}

vec3 &BBox::operator[](int i) {
	assert(i == 0 || i == 1);
	if (i == 0) {
		return min;
	}
	else {
		return max;
	}
}

void BBox::expandToInclude(const glm::vec3& p) {
	min = glm::min(min, p);
	max = glm::max(max, p);
}

void BBox::expandToInclude(const BBox& b) {
	min = glm::min(min, b.min);
	max = glm::max(max, b.max);
}

short BBox::maxDimension() const {
	short result = 0;
	vec3 extent = this->max - this->min;

	if (extent.y > extent.x) result = 1;
	if (extent.z > extent.y) result = 2;
	return result;
}

float BBox::surfaceArea() const {
	vec3 extent = this->max - this->min;
	return 2.0f * (extent.x * extent.y + extent.z * extent.y + extent.x * extent.z);
}

vec3 BBox::getCenter() const {
	return (min + max) * 0.5f;
}

float BBox::distance(const BBox& other) {
	float result = 0;
	for (length_t i = 0; i < 3; ++i) {
		const float amin = min[i];
		const float amax = max[i];
		const float bmin = other.min[i];
		const float bmax = other.max[i];

		if (amin > bmax) {
			float delta = bmax - amin;
			result += delta * delta;
		}
		else if (bmin > amax) {
			float delta = amax - bmin;
			result += delta * delta;
		}
	}
	return std::sqrt(result); //todo why not glm::sqrt? precision?
}

// https://tavianator.com/fast-branchless-raybounding-box-intersections-part-2-nans/
bool BBox::intersect(const Ray& ray) const {
	vec3 inv_dir = vec3(1.0f) / ray.direction;
	double tx0 = (min.x - ray.origin.x) * inv_dir.x;
	double tx1 = (max.x - ray.origin.x) * inv_dir.x;

	double tmin = std::min(tx0, tx1);
	double tmax = std::max(tx0, tx1);

	double ty0 = (min.y - ray.origin.y) * inv_dir.y;
	double ty1 = (max.y - ray.origin.y) * inv_dir.y;

	tmin = std::max(tmin, std::min(ty0, ty1));
	tmax = std::min(tmax, std::max(ty0, ty1));

	double tz0 = (min.z - ray.origin.z) * inv_dir.z;
	double tz1 = (max.z - ray.origin.z) * inv_dir.z;

	tmin = std::max(tmin, std::min(tz0, tz1));
	tmax = std::min(tmax, std::max(tz0, tz1));

	if (tmax < std::max(tmin, 0.0))
		return false;

	return true;
}

#if QUALITY_OVER_PERFORMANCE

// reference fps: 12.7fps
bool BBox::intersect(const Ray& ray, float *tnear, float *tfar) const {
	float t0 = -INFINITY;
	float t1 = INFINITY;

	// Loop over the three axes and compute the hit time for the
	// two axis-aligned bounding box planes in each, decreasing the
	// parametric range of the ray until start>end time, which means
	// the ray missed the box, or until we finish which means there
	// is an intersection.
	for (int i = 0; i < 3; i++) {
		float invDir = 1.0f / ray.direction[i]; // suprisingly faster to do calculation here, than precompute ray inverse direction (probably because of cache)
		float tNear = (min[i] - ray.origin[i]) * invDir;
		float tFar = (max[i] - ray.origin[i]) * invDir;

		if (tNear > tFar) swap(tNear, tFar);

		if (tNear > t0) t0 = tNear;
		if (tFar < t1) t1 = tFar;
		if (t0 > t1) return false;
	}

	*tnear = t0;
	*tfar = t1;
	return true;
}

#else

// reference fps: 17.5fps, but has some weird moments due to glm::min/max and SSE precision
bool BBox::intersect(const Ray& ray, float *tnear, float *tfar) const {
	// Smits ray-box intersection test using slabs
	// http://www.cs.utah.edu/~awilliam/box/box.pdf
	// *** SSE ***
	__m128 min_ps = _mm_set_ps(0.0f, min[2], min[1], min[0]);
	__m128 max_ps = _mm_set_ps(0.0f, max[2], max[1], max[0]);
	__m128 position_ps = _mm_set_ps(0.0f, ray.origin[2], ray.origin[1], ray.origin[0]);
	__m128 div_ps = _mm_set_ps(0.0f, 1.0f / ray.direction[2], 1.0f / ray.direction[1], 1.0f / ray.direction[0]);

	// Compute intersections
	__m128 t1_ps = _mm_mul_ps(_mm_sub_ps(min_ps, position_ps), div_ps);
	__m128 t2_ps = _mm_mul_ps(_mm_sub_ps(max_ps, position_ps), div_ps);

	__m128 tmin_ps = _mm_min_ps(t1_ps, t2_ps);
	__m128 tmax_ps = _mm_max_ps(t1_ps, t2_ps);
	float *tmin_a = (float *)&tmin_ps;
	float *tmax_a = (float *)&tmax_ps;

	if (tmin_a[0] > tmax_a[1] || tmin_a[1] > tmax_a[0])
		return false;

	float tmin = glm::max(tmin_a[0], tmin_a[1]);
	float tmax = glm::min(tmax_a[0], tmax_a[1]);

	if (tmin > tmax_a[2] || tmin_a[2] > tmax)
		return false;

	tmin = glm::max(tmin, tmin_a[2]);
	tmax = glm::min(tmax, tmax_a[2]);

	*tnear = tmin;
	*tfar = tmax;

	return tmax > 0.0f;
}

#endif
