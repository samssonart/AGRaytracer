#pragma once
#include "Renderer.h"
#include "Scene.h"

class RayCPU : public Renderer {

	// Members
	Scene* scene = NULL;
	
	// Ray tracing functions
	glm::vec3 castRay(Ray& ray) const {
		Intersection intersection = getClosestIntersection(ray, false);

		if (intersection.didIntersect) {
			return performLighting(intersection);
		}
		else {
			return glm::vec3(scene->background_color[0], scene->background_color[1], scene->background_color[2]);
			// no intersections, out of scene
		}
	}

	// Returns closest intersection, handles accelerator selection
	// If isForShadow is true, returns true on first found intersection
	Intersection getClosestIntersection(Ray& ray, bool isForShadow) const {
		return this->accelerators[this->current_accelerator_id]->getIntersection(ray, false);
	}

	// Returns true if given ray doesn't reach light and false otherwise
	bool isInShadow(Ray& ray, float light_distance) const {
		Intersection intersection = getClosestIntersection(ray, true);

		return intersection.didIntersect && intersection.distance < light_distance;
	}

	// Performs shading of the intersected point
	// Combining result color from ambient (can be 0), diffuse and specular, and reflected and refracted color
	glm::vec3 performLighting(Intersection& intersection) const {

		glm::vec3 color = scene->materials[intersection.object->matId]->getColor(intersection.point);

		glm::vec3 ambientColor = color * scene->ambient_light_coeff;

		glm::vec3 diffuseAndSpecularColor = getDiffuseAndSpecularLighting(intersection);
		glm::vec3 reflectedColor = getReflectiveRefractiveLighting(intersection);

		// Phong equation = ambient + diffuse + specular
		return glm::clamp(ambientColor + diffuseAndSpecularColor + reflectedColor, 0.0f, 1.0f);
	}

	// Returns color of diffuse and specular components at given intersection point
	glm::vec3 getDiffuseAndSpecularLighting(Intersection& intersection) const {

		vec3 diffuse_color(0.0f);
		vec3 specular_color(0.0f);
		bool shadowed = false;

		for (Light* light : scene->lights) {
			vec3 light_offset = light->location - intersection.point;
			float light_distance = length(light_offset);

			vec3 light_direction = normalize(light_offset);
			float dot_product = dot(light_direction, intersection.normal);

			if (dot_product >= 0.0f) {
				if (scene->cast_shadows) {

					Ray shadowRay(intersection.point + (light_direction * 2.0f * scene->epsilon), light_direction, 1, intersection.ray.matId);

					shadowed = isInShadow(shadowRay, light_distance);
					if (shadowed) {
						continue;
					}
				}

				diffuse_color = diffuse_color + ((scene->materials[intersection.endMaterial]->getColor(intersection.point) * dot_product) * light->intensity);
				specular_color = specular_color + getSpecularLighting(intersection, light);
			}
		}

		return diffuse_color + specular_color;
	}

	// Returns direction of reflection ray
	glm::vec3 reflectVector(glm::vec3 vector, glm::vec3 normal) const {
		return glm::normalize(2.0f * normal * glm::dot(normal, vector) - vector);
	}

	// Returns direction for refraction ray
	glm::vec3 refractVector(const vec3& normal, const vec3& incident, float n1, float n2) const {
		// http://steve.hollasch.net/cgindex/render/refraction.txt
		float eta = n1 / n2;
		float c1 = -dot(normal, incident);
		float cs2 = 1 - eta * eta * (1 - c1 * c1);
		if (cs2 < 0) {
			// total internal reflection 
			return vec3(0.0f);
		}
		glm::vec3 refractedVector = eta * incident + (eta * c1 - sqrtf(cs2)) * normal;
		return refractedVector;
	}

	// Returns color of specular component at given intersection point from given light
	vec3 getSpecularLighting(Intersection& intersection, Light* light) const {
		vec3 specular_color(0.0f);

		float shininess = scene->materials[intersection.endMaterial]->getShininess();

		if (shininess == NOT_SHINY) {
			return specular_color;
		}

		vec3 view = normalize(intersection.ray.origin - intersection.point);
		vec3 light_offset = light->location - intersection.point;
		vec3 reflected = reflectVector(normalize(light_offset), intersection.normal);

		float dot_product = dot(reflected, view);

		if (dot_product <= 0) {
			return specular_color;
		}

		float specular_amount = pow(dot_product, shininess) * light->intensity;
		specular_color = vec3(specular_amount);

		return specular_color;
	}

	// Returns reflectance either by Schlick's approximation or original Fresnel equaions
	float getReflectance(const vec3& normal, const vec3& incident, float n1, float n2) const {

		//// Fresnel equations
		//float n = n1 / n2;
		//float cosI = -dot(normal, incident);
		//float sinT2 = n * n * (1.0f - cosI * cosI);
		//if (sinT2 > 1.0f) {
		//	// Total Internal Reflection.
		//	return 1.0f;
		//}
		//float cosT = sqrtf(1.0f - sinT2);
		//// https://upload.wikimedia.org/math/a/b/8/ab8211bd07c24f9f7ccf9d0503e582c7.png
		//float Rs = (n1 * cosI - n2 * cosT) / (n1 * cosI + n2 * cosT);
		//Rs *= Rs;
		//// https://upload.wikimedia.org/math/f/2/c/f2c7ba9030a97a3e9f43bd2973c0eea0.png
		//float Rp = (n2 * cosI - n1 * cosT) / (n2 * cosI + n1 * cosT);
		//Rp *= Rp;
		//// if the incident light is unpolarised (containing an equal mix of s- and p-polarisations)
		//// the reflectance is sum of both divided by two
		//return (Rs + Rp) / 2.0f;

		// Computing reflection coefficent with Schlick's approximation 
		// https://en.wikipedia.org/wiki/Schlick%27s_approximation
		float r0 = ((n1 - n2) / (n1 + n2));
		r0 *= r0;
		float r = r0 + (1 - r0) * powf(1 - dot(normalize(-normal), incident), 5.0f);
		return clamp(r, 0.0f, 1.0f);
	}

	// Returns color of reflection and refraction at given intersection point
	vec3 getReflectiveRefractiveLighting(Intersection& intersection) const {
		float reflectivity = scene->materials[intersection.endMaterial]->getReflectivity();
		float startRefractiveIndex = scene->materials[intersection.startMaterial]->getRefractiveIndex();
		float endRefractiveIndex = scene->materials[intersection.endMaterial]->getRefractiveIndex();
		int reflectionsRemaining = intersection.ray.reflections_remaining;

		if ((reflectivity == NOT_REFLECTIVE && endRefractiveIndex == NOT_REFRACTIVE) || reflectionsRemaining <= 0) {
			return vec3(0.0f);
		}

		float reflectivePercentage = reflectivity;
		float refractivePercentage = 0.0f;

		if (endRefractiveIndex != NOT_REFRACTIVE) {
			reflectivePercentage = getReflectance(intersection.normal, intersection.ray.direction, startRefractiveIndex, endRefractiveIndex);
			refractivePercentage = 1 - reflectivePercentage;
		}

		if (refractivePercentage <= 0 && reflectivePercentage <= 0) {
			return vec3(0.0f);
		}

		vec3 reflectiveColor(0.0f);
		vec3 refractiveColor(0.0f);

		if (reflectivePercentage > 0) {
			vec3 reflected = reflectVector(normalize(intersection.ray.origin - intersection.point), intersection.normal);
			Ray reflectedRay(intersection.point + (reflected * scene->epsilon), reflected, reflectionsRemaining - 1, intersection.ray.matId);

			reflectiveColor = castRay(reflectedRay) * reflectivePercentage;
		}

		if (refractivePercentage > 0) {
			vec3 refracted = refractVector(intersection.normal, intersection.ray.direction, startRefractiveIndex, endRefractiveIndex);
			Ray refractedRay = Ray(intersection.point + (refracted * scene->epsilon), refracted, 1, intersection.endMaterial);

			refractiveColor = castRay(refractedRay) * refractivePercentage;
		}

		return reflectiveColor + refractiveColor;
	}

public:

	~RayCPU() {
		// Freeing up acceleration structures
		for (Accelerator* acc : this->accelerators) {
			delete acc;
		}
		this->accelerators.clear();
	}

	// Populates scene pointer 
	RayCPU(Scene* scene) {
		this->scene = scene;
		this->name = "CPU Ray Tracing";
		cout << "[Renderer] Initializing " << this->name << " renderer.\n";

		// Building renderers acceleration structures
		this->accelerators.push_back(new NoAcceleration(scene)); // brute force intersection
		this->accelerators.push_back(new FastBVH(scene)); // Fast-BVH by Brandon Pelfrey, https://github.com/brandonpelfrey/Fast-BVH

		// PBRTv2 by Matt Pharr & Greg Humphreys, http://pbrt.org/
		this->accelerators.push_back(new PBRT(scene, SPLIT_EQUAL_COUNTS, "PBRT2 Grid")); // EQUAL_COUNTS, grid
		this->accelerators.push_back(new PBRT(scene, SPLIT_MIDDLE, "PBRT2 Middle")); // SPLIT_MIDDLE, middle point of centroids along the longest axis
		this->accelerators.push_back(new PBRT(scene, SPLIT_SAH, "PBRT2 SAH BIN")); // SPLIT_SAH, SAH with binning
		
		// Approximate Agglomerative Clustering, on base of PBRTv2 BVH Builder
		this->accelerators.push_back(new PBRT(scene, SPLIT_AAC_HQ, "PBRT2 AAC HQ")); // SPLIT_AAC_HQ
		this->accelerators.push_back(new PBRT(scene, SPLIT_AAC_FAST, "PBRT2 AAC Fast")); // SPLIT_AAC_FAST

		// Luxrays BVH builder
		this->accelerators.push_back(new LuxBVH(scene));

		// Setting default to PBRT HQ
		this->current_accelerator_id = 5;

		// Building UI string of accelerators
		for (Accelerator* a : this->accelerators) {
			accelerators_ui_string += a->name;
			accelerators_ui_string.push_back('\0');
		}
		accelerators_ui_string.push_back('\0');

		// Reporting of structures build times
		for (Accelerator* a : this->accelerators) {
			a->printStats();
		}

		cout << "[Renderer] Finished initialization of " << this->name << " renderer.\n\n";
	}

	// Returns color of the pixel with X/Y coords, handles supersampling
	glm::vec3 trace_ray(int x, int y) const {

		if (scene == NULL) {
			std::cout << "[RayCPU] Initialize the renderer first!";
			exit(EXIT_FAILURE);
		}

		glm::vec3 result_color(0.0f);

		float x_origin = x - 1.0f;
		float y_origin = y - 1.0f;

		float super_width = 1.0f / scene->supersampling;
		float super_coeff = 1.0f / (scene->supersampling * scene->supersampling);

		for (int i = 1; i < scene->supersampling + 1; i++) {
			for (int j = 1; j < scene->supersampling + 1; j++) {
				float normalized_i = ((x_origin + (super_width * i)) / scene->camera.width) - 0.5f;
				float normalized_j = ((y_origin + (super_width * j)) / scene->camera.height) - 0.5f;

				vec3 image_point = scene->camera.dir +
					normalized_i * scene->camera.right * scene->camera.aspect_ratio +
					normalized_j * scene->camera.up;

				Ray viewRay(scene->camera.pos, image_point, scene->max_reflections, 0); // 0 is index of Air, starting material
				result_color += castRay(viewRay) * super_coeff;
			}
		}

		return result_color;

	}

	// Traces all primary rays and plots them on screen
	void render(Tmpl8::Surface* screen) {
	
		#pragma omp parallel for
		for (int x = 0; x < screen->GetWidth(); x++) {
			for (int y = 0; y < screen->GetHeight(); y++) {
				screen->Plot(x, y, toPixel(trace_ray(x, y)));
			}
		}
	}

}; 