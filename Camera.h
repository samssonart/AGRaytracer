#pragma once
#include "glm/glm.hpp"
#include <glm/gtx/rotate_vector.hpp>
using namespace glm;

static const float look_speed = 6.0f;
static const float move_speed = 2.0f;

class Camera {
public:
	vec3 pos;
	vec3 dir;
	vec3 up;
	vec3 right;
	float aspect_ratio, width, height;

	Camera() {}

	Camera(vec3 p, vec3 d) {
		pos = p;
		dir = d;
		up = vec3(0.0f, -1.0f, 0.0f);
		aspect_ratio = (float)SCRWIDTH / (float)SCRHEIGHT;
		width = (float)SCRWIDTH;
		height = (float)SCRHEIGHT;

		update();
	}

	void update() {
		right = cross(dir, up);
		up = cross(right, dir);
	}

	// movement methods
	void moveForward(float distance) {
		this->pos += this->dir * distance * move_speed;
	}
	void moveBackward(float distance) {
		moveForward(-distance);
	}
	void moveRight(float distance) {
		vec3 right = cross(this->dir, this->up);
		this->pos += right * distance * move_speed;
	}
	void moveLeft(float distance) {
		moveRight(-distance);
	}
	void moveUp(float distance) {
		this->pos += this->up * (-distance) * move_speed;
	}
	void moveDown(float distance) {
		moveUp(-distance);
	}

	// looking methods
	void lookDown(float degrees) {
		vec3 right = cross(this->dir, this->up);
		this->dir = rotate(this->dir, glm::radians(degrees * look_speed), right);
	}
	void lookUp(float degrees) {
		lookDown(-degrees);
	}
	void lookLeft(float degrees) {
		this->dir = rotate(this->dir, glm::radians(degrees * look_speed), this->up);
	}
	void lookRight(float degrees) {
		lookLeft(-degrees);
	}

	void reset() {
		pos = vec3(0, 30, -168);
		dir = vec3(0, 0, 1);
		up = vec3(0.0f, -1.0f, 0.0f);
	}


};
