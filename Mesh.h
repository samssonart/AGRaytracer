#pragma once
#include <vector>
#include "Primitive.h"
#include "Triangle.h"
#include "Flat.h"
#include "Shiny.h"

class Mesh : public Primitive {
public:
	BBox bbox;
	vec3 center;
	string name;
	int primitive_count;
	int first_index;

	Mesh(std::vector<Primitive*>* triangles,
		std::vector<float> positions, std::vector<float> normals, std::vector<float> texcoords,
		std::vector<unsigned int> indices, std::vector<unsigned char> num_vertices, 
		std::vector<int> material_ids, std::string name) {

		// type of primitive
		this->type = MESH;

		// name, passed from .obj file
		this->name = name;

		// setting index of the first mesh triangle in geometry buffer
		first_index = (int)(*triangles).size();

		// holding vectors
		std::vector<vec3> vertices;
		std::vector<vec3> vert_normals;

		// populating vertices and vertex normals if available
		for (size_t v = 0; v < positions.size() / 3; v++) {
			vertices.push_back(vec3(positions[3 * v + 0], positions[3 * v + 1], positions[3 * v + 2]));
			if (normals.size() > 0) {
				vert_normals.push_back(vec3(normals[3 * v + 0], normals[3 * v + 1], normals[3 * v + 2]));
			}
		}

		for (size_t f = 0; f < num_vertices.size(); f++) {
			// vertices for a face
			vec3 a = vertices[indices[3 * f + 0]];
			vec3 b = vertices[indices[3 * f + 1]];
			vec3 c = vertices[indices[3 * f + 2]];

			// normals
			vec3 a_norm, b_norm, c_norm;
			if (normals.size() > 0) {
				// vertex normals found in .obj
				a_norm = vert_normals[indices[3 * f + 0]];
				b_norm = vert_normals[indices[3 * f + 1]];
				c_norm = vert_normals[indices[3 * f + 2]];
			}

			// Materials
			// TODO add smoothed to constructor
			int material = material_ids[f] + 1; // +1 because zeroth material is always air, starting material
			if (material_ids[f] < 0) material = 1;

			// Assembling the triangle
			Triangle* tri = new Triangle(a, b, c, a_norm, b_norm, c_norm, material);
			tri->parent = this;

			(*triangles).push_back(tri);
		}

		// if vertex normals are not in file, generate them here
		if (normals.size() == 0) {
			cout << "[OBJ] Mesh " << name << " didn't have vertex normals in .obj file. Generating them.\n";
			vert_normals.resize(vertices.size());

			for (size_t f = 0; f < num_vertices.size(); f++) {
				// normalized and with appropriate weights, depending on area of triangles
				vert_normals[indices[3 * f + 0]] += normalize(((Triangle*)(*triangles)[first_index + f])->face_normal * ((Triangle*)(*triangles)[first_index + f])->getArea());
				vert_normals[indices[3 * f + 1]] += normalize(((Triangle*)(*triangles)[first_index + f])->face_normal * ((Triangle*)(*triangles)[first_index + f])->getArea());
				vert_normals[indices[3 * f + 2]] += normalize(((Triangle*)(*triangles)[first_index + f])->face_normal * ((Triangle*)(*triangles)[first_index + f])->getArea());
			}

			for (size_t f = 0; f < num_vertices.size(); f++) {
				vec3 a_norm = vert_normals[indices[3 * f + 0]];
				vec3 b_norm = vert_normals[indices[3 * f + 1]];
				vec3 c_norm = vert_normals[indices[3 * f + 2]];

				((Triangle*)(*triangles)[first_index + f])->a_norm = a_norm;
				((Triangle*)(*triangles)[first_index + f])->b_norm = b_norm;
				((Triangle*)(*triangles)[first_index + f])->c_norm = c_norm;
			}
		}

		// creating bbox
		for (Primitive* tr : (*triangles)) {
			bbox.expandToInclude(tr->getBBox());
		}

		// computing center
		center = (bbox.max - bbox.min) / 2.0f;

		// primitive count
		this->primitive_count = (int)num_vertices.size();
	}

	~Mesh() {
	}

	Intersection intersect(const Ray& ray) {
		// If everything is correct, it should never be invoked directly
		cout << "[Core] Something is trying to intersect(ray) with Mesh class directly, check it out.\n";
		return Intersection();
	}

	vec3 getNormal(const vec3& point) const {
		// If everything is correct, it should never be invoked directly
		cout << "[Core] Something is trying to getNormal(vec3) with Mesh class directly, check it out.\n";
		return vec3(0.0f);
	}

	BBox getBBox() const {
		return this->bbox;
	}

	vec3 getCentroid() const {
		return this->center;
	}

	int getPrimitiveCount() const {
		return this->primitive_count;
	}


};