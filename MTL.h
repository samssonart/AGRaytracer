#pragma once
#include "Material.h"

class MTL : public Material {
private:
	glm::vec3 diffuse = glm::vec3(0.0f);
	float refractive_index = NOT_REFRACTIVE;
	float shininess = NOT_SHINY;
	float reflectivity = NOT_REFLECTIVE;

public:
	MTL(float ambient[3], float diffuse[3], float specular[3], float transmittance[3], float emission[3],
		float shininess, float ior, float dissolve, int illum) {
		
		// MTL specifications, useful links
		// http://paulbourke.net/dataformats/mtl/
		// http://www.fileformat.info/format/material/

		// TODO: improve material system, implement all types of materials instead of simplification
		// in a perfect world:
		//this->ambient = glm::vec3(ambient[0], ambient[1], ambient[2]);
		//this->diffuse = glm::vec3(diffuse[0], diffuse[1], diffuse[2]);
		//this->specular = glm::vec3(specular[0], specular[1], specular[2]);
		//this->transmittance = glm::vec3(transmittance[0], transmittance[1], transmittance[2]);
		//this->emission = glm::vec3(emission[0], emission[1], emission[2]);
		//this->shininess = shininess;
		//this->ior = ior;
		//this->dissolve = dissolve;
		//this->illum = illum;



		switch (illum) {
		default:
		case 0: // This is a constant color illumination model.
		case 1: // This is a diffuse illumination model using Lambertian shading.
			this->diffuse = glm::vec3(diffuse[0], diffuse[1], diffuse[2]);
			break;


		case 2: // This is a diffuse and specular illumination model using Lambertian 
				// shading and Blinn's interpretation of Phong's specular illumination
				// model(BLIN77).
			this->diffuse = glm::vec3(diffuse[0], diffuse[1], diffuse[2]);
			this->shininess = 100.0f;
			break;


		case 3: // This is a diffuse and specular illumination model with reflection 
				// using Lambertian shading, Blinn's interpretation of Phong's specular
				// illumination model(BLIN77), and a reflection term similar to that in
				// Whitted's illumination model (WHIT80).
		case 8: // This illumination model is similar to illumination model 3 without 
				// ray tracing.
			this->diffuse = glm::vec3(diffuse[0], diffuse[1], diffuse[2]);
			this->shininess = 100.0f;
			this->reflectivity = 0.5f;
			break;


		case 4: //The diffuse and specular illumination model used to simulate glass 
				// is the same as illumination model 3.
		case 5: // This is a diffuse and specular shading models similar to 
				// illumination model 3, except that reflection due to Fresnel effects is 
				// introduced into the equation.
		case 6: // This is a diffuse and specular illumination model similar to that 
				// used by Whitted (WHIT80) that allows rays to refract through a surface.
		case 7: // This illumination model is similar to illumination model 6, except 
				// that reflection and transmission due to Fresnel effects has been
				// introduced to the equation.
		case 9: // This illumination model is similar to illumination model 4without 
				// ray tracing.
		case 10:// This illumination model is used to cast shadows onto an invisible 
				// surface.
			this->diffuse = glm::vec3(0.0f);
			this->shininess = 50.0f;
			this->refractive_index = 2.0f;
			break;
		}
	}

	glm::vec3 getColor(glm::vec3 point) {
		return this->diffuse;
	}

	glm::vec3 getColor() {
		return this->diffuse;
	}

	float getShininess() {
		return shininess;
	}

	float getReflectivity() {
		return reflectivity;
	}

	float getRefractiveIndex() {
		return refractive_index;
	}

};