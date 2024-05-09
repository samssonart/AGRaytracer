#pragma once
#include <CL\cl.hpp>
#include "Renderer.h"
#include "Scene.h"

class Primitive;
class Triangle;
class BVH;
class Light;

typedef struct CL_General {
	cl_int triangle_count;
	cl_int lights_count;

	cl_bool cast_shadows;
	cl_float ambient_light;
	cl_int supersampling;
	cl_int max_reflections;
	cl_float3 background_color;
	cl_float epsilon;
} CL_General;

typedef struct CL_Ray {
	cl_float3 origin;
	cl_float3 direction;
	cl_float3 inverse_direction;
	cl_bool fast_occlusion;
	cl_int reflections_remaining;
} CL_Ray;

typedef struct CL_Intersection {
	cl_bool did_intersect;
	cl_float distance;
	cl_float3 position;
	cl_float3 normal;
	cl_float3 bary;
	cl_int mat_id;
} CL_Intersection;

typedef struct CL_Camera {
	cl_float3 position;
	cl_float3 direction;
	cl_float3 right;
	cl_float3 up;
	cl_int width, height;
	cl_float aspect_ratio;
} CL_Camera;

typedef struct CL_PointLight {
	cl_float3 position;
	cl_float3 color;
	cl_float intensity;
} CL_PointLight;

typedef struct CL_Material {
	cl_float3 diffuse;
	cl_float refractive_index;
	cl_float shininess;
	cl_float reflectivity;
} CL_Material;

typedef struct CL_Triangle {
	cl_float3 vertices[3];
	cl_float3 vert_normals[3];
	cl_float3 normal;
	cl_int mat_id;
} CL_Triangle;

typedef struct CL_BVHNode {
	cl_float3 bounds[2];
	union {
		cl_uint primitives_offset;
		cl_uint second_child_offset;
	};
	cl_ushort n_primitives;
	cl_ushort axis;
} CL_BVHNode;

typedef struct CL_BVHAccelArrayNode {
	cl_float3 bounds[2];
	cl_uint primitive;
	cl_uint skipIndex;
} CL_BVHAccelArrayNode;


class RayGPU : public Renderer {

protected:

	// Pointer to scene
	Scene* scene = NULL;
	
	// Error handling
	cl_int err;

	// Device related
	vector<cl::Platform> platforms;
	cl::Platform default_platform;
	cl::Device default_device;
	bool hasDedicatedGPU = false; // TODO?
	string used_device_name, used_platform_name;

	// Buffers
	cl::Buffer general_buffer;
	cl::Buffer camera_buffer;
	cl::Buffer lights_buffer;
	cl::Buffer triangles_buffer;
	cl::Buffer bvh_buffer;
	cl::Buffer screen_buffer;
	cl::Buffer materials_buffer;

	cl::Buffer rays_buffer;
	cl::Buffer intersections_buffer;

	// Context, Program and Queue
	cl::Context context;
	cl::Program program;
	cl::Program::Sources sources;
	cl::CommandQueue queue;

	// Kernels
	cl::Kernel raytracing_kernel;

	cl::Kernel primary_ray_generator_kernel;
	cl::Kernel intersection_kernel;
	cl::Kernel shading_kernel;
	
	// Loaders
	void loadProgram(vector<string> source_files);
	void loadGeneral();
	void loadMaterials();
	void loadCamera();
	void loadLights();
	void loadTriangles();
	void loadBVH();

	// Helpers
	// Converts vec3 to cl_float3
	cl_float3 vec3tofloat3(vec3 vector) {

		cl_float3 value;

		value.s[0] = vector.x;
		value.s[1] = vector.y;
		value.s[2] = vector.z;
		value.s[3] = 0.0f;

		return value;
	}

	// Returns OpenCL error code translation 
	string getOpenCLErrorMessage(cl_int err) {
		switch (err)
		{
		case 0: return "CL_SUCCESS";
		case -1: return "CL_DEVICE_NOT_FOUND";
		case -2: return "CL_DEVICE_NOT_AVAILABLE";
		case -3: return "CL_COMPILER_NOT_AVAILABLE";
		case -4: return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
		case -5: return "CL_OUT_OF_RESOURCES";
		case -6: return "CL_OUT_OF_HOST_MEMORY";
		case -7: return "CL_PROFILING_INFO_NOT_AVAILABLE";
		case -8: return "CL_MEM_COPY_OVERLAP";
		case -9: return "CL_IMAGE_FORMAT_MISMATCH";
		case -10: return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
		case -11: return "CL_BUILD_PROGRAM_FAILURE";
		case -12: return "CL_MAP_FAILURE";
		case -13: return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
		case -14: return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
		case -15: return "CL_COMPILE_PROGRAM_FAILURE";
		case -16: return "CL_LINKER_NOT_AVAILABLE";
		case -17: return "CL_LINK_PROGRAM_FAILURE";
		case -18: return "CL_DEVICE_PARTITION_FAILED";
		case -19: return "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";
		case -30: return "CL_INVALID_VALUE";
		case -31: return "CL_INVALID_DEVICE_TYPE";
		case -32: return "CL_INVALID_PLATFORM";
		case -33: return "CL_INVALID_DEVICE";
		case -34: return "CL_INVALID_CONTEXT";
		case -35: return "CL_INVALID_QUEUE_PROPERTIES";
		case -36: return "CL_INVALID_COMMAND_QUEUE";
		case -37: return "CL_INVALID_HOST_PTR";
		case -38: return "CL_INVALID_MEM_OBJECT";
		case -39: return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
		case -40: return "CL_INVALID_IMAGE_SIZE";
		case -41: return "CL_INVALID_SAMPLER";
		case -42: return "CL_INVALID_BINARY";
		case -43: return "CL_INVALID_BUILD_OPTIONS";
		case -44: return "CL_INVALID_PROGRAM";
		case -45: return "CL_INVALID_PROGRAM_EXECUTABLE";
		case -46: return "CL_INVALID_KERNEL_NAME";
		case -47: return "CL_INVALID_KERNEL_DEFINITION";
		case -48: return "CL_INVALID_KERNEL";
		case -49: return "CL_INVALID_ARG_INDEX";
		case -50: return "CL_INVALID_ARG_VALUE";
		case -51: return "CL_INVALID_ARG_SIZE";
		case -52: return "CL_INVALID_KERNEL_ARGS";
		case -53: return "CL_INVALID_WORK_DIMENSION";
		case -54: return "CL_INVALID_WORK_GROUP_SIZE";
		case -55: return "CL_INVALID_WORK_ITEM_SIZE";
		case -56: return "CL_INVALID_GLOBAL_OFFSET";
		case -57: return "CL_INVALID_EVENT_WAIT_LIST";
		case -58: return "CL_INVALID_EVENT";
		case -59: return "CL_INVALID_OPERATION";
		case -60: return "CL_INVALID_GL_OBJECT";
		case -61: return "CL_INVALID_BUFFER_SIZE";
		case -62: return "CL_INVALID_MIP_LEVEL";
		case -63: return "CL_INVALID_GLOBAL_WORK_SIZE";
		case -64: return "CL_INVALID_PROPERTY";
		case -65: return "CL_INVALID_IMAGE_DESCRIPTOR";
		case -66: return "CL_INVALID_COMPILER_OPTIONS";
		case -67: return "CL_INVALID_LINKER_OPTIONS";
		case -68: return "CL_INVALID_DEVICE_PARTITION_COUNT";
		case -1000: return "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR";
		case -1001: return "CL_PLATFORM_NOT_FOUND_KHR";
		case -1002: return "CL_INVALID_D3D10_DEVICE_KHR";
		case -1003: return "CL_INVALID_D3D10_RESOURCE_KHR";
		case -1004: return "CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR";
		case -1005: return "CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR";
		default: return "Unknown error";
		}
	}

	// If something is wrong throws error code with translation and exits application
	void checkError(cl_int err, string message) {
		if (err != CL_SUCCESS) {
			cout << "[OpenCL] Error " << err << ": " << message << " - " << getOpenCLErrorMessage(err) << "\n";
			system("pause");
			exit(EXIT_FAILURE);
		}
	}

public:

	RayGPU(Scene* scene);
	~RayGPU();
	void render(Tmpl8::Surface* screen);

};