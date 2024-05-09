float3 reflect_vector(float3 direction, float3 normal) {
	return normalize(2.0f * normal * dot(normal, direction) - direction);
}

float3 get_reflectance(float3 normal, float3 incident, float n1, float n2) {
	float r0 = ((n1 - n2) / (n1 + n2));
	r0 *= r0;
	float r = r0 + (1 - r0) * pow(1 - dot(normalize(-normal), incident), 5.0f);
	return clamp(r, 0.0f, 1.0f);
}

float3 performLighting(
	__global const General* general,
	__global const BVHNode* bvh,
	__global const Triangle* triangles,
	__global const PointLight* lights,
	__global const Material* materials,
	const Ray ray,
	const Intersection intersection) {

	// material diffuse color
	float3 color = materials[intersection.mat_id].diffuse;
	
	// ambient component, depending on runtime changeable ambient light coefficent
	float3 ambient_color = color * general->ambient_light;

	// diffuse and specular
	float3 diffuse_color = (float3)(0.0f, 0.0f, 0.0f);
	float3 specular_color = (float3)(0.0f, 0.0f, 0.0f);

	for (uint i = 0; i < general->lights_count; i++) {

		PointLight light = lights[i];
		float3 light_offset = light.position - intersection.position;
		float light_distance = length(light_offset);
		float3 light_direction = normalize(light_offset);

		float dot_product = dot(light_direction, normalize(intersection.normal));

		if (dot_product >= 0.0f) {
			if (general->cast_shadows) {

				Ray shadow_ray = {
					intersection.position + (light_direction * general->epsilon),
					light_direction,
					1.f / light_direction,
					true, // fast occlusion flag, stops searching intersections after first hit
					1
				};

				Intersection shadow_intersection = getClosestIntersection(bvh, triangles, shadow_ray);

				if (shadow_intersection.did_intersect) {
					continue;
				}
			}

			diffuse_color += (materials[intersection.mat_id].diffuse * dot_product) * light.intensity;
			
			// specular part
			float shininess = materials[intersection.mat_id].shininess;

			if (shininess == -1.0f) {
				continue;
			}

			float3 view = normalize(ray.origin - intersection.position);
			float3 reflected = reflect_vector(light_direction, intersection.normal);

			float reflected_dot = dot(reflected, view);

			if (reflected_dot <= 0) {
				continue;
			}

			float specular_amount = pow(reflected_dot, shininess) * light.intensity;
			specular_color += (float3)(specular_amount, specular_amount, specular_amount);

		}
	}

	// reflection and refraction
	float3 reflective_color = (float3)(0.0f, 0.0f, 0.0f);
	//float3 refractive_color = (float3)(0.0f, 0.0f, 0.0f);

	float reflectivity = materials[intersection.mat_id].reflectivity;

	if (reflectivity > 0.0f && ray.reflections_remaining > 0) {
		float3 reflected = reflect_vector(normalize(ray.origin - intersection.position), intersection.normal);
		Ray reflected_ray = { intersection.position + (reflected * general->epsilon), reflected, 1.f / reflected, false, ray.reflections_remaining - 1 };
		
		reflective_color = cast_ray(general, lights, triangles, bvh, materials, reflected_ray) * reflectivity;
	}

	return clamp(ambient_color + diffuse_color
		+ specular_color + reflective_color /*+ refractive_color*/,
		0.0f, 1.0f);
}

