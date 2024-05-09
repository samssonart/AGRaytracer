#include "precomp.h"
#define WAVEFRONT 0

/// Loaders
// Loads kernels specified in vector argument and builds them
void RayGPU::loadProgram(vector<string> source_files) {
	cout << "[OpenCL] Building OpenCL program.\n";

	// pushing all CL source files into one string stream
	stringstream source_files_ss;
	for (const string& source_file : source_files) {
		ifstream file(source_file);

		if (!file.good()) {
			cout << "[OpenCL] Error reading OpenCL file: " << source_file << "\n";
			system("pause");
			exit(EXIT_FAILURE);
		}

		source_files_ss << file.rdbuf();
		file.close();
	}

	string combined_source = source_files_ss.str();
	sources.push_back({ combined_source.c_str(), combined_source.length() });


	program = cl::Program(context, sources);
	err = program.build({ default_device }, "-cl-nv-verbose");
	string build_log = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device);

	if (!build_log.empty())
		cout << "[OpenCL] Build log:\n" << build_log  << "\n";

	checkError(err, "Can't build kernels!");
}

// Fills materials_buffer with CL_Material objects
void RayGPU::loadMaterials() {
	size_t material_count = scene->materials.size();
	CL_Material* materials = new CL_Material[material_count];

	for (size_t i = 0; i < material_count; i++) {
		materials[i].diffuse = vec3tofloat3(scene->materials[i]->getColor());
		materials[i].reflectivity = scene->materials[i]->getReflectivity();
		materials[i].refractive_index = scene->materials[i]->getRefractiveIndex();
		materials[i].shininess = scene->materials[i]->getShininess();
	}

	checkError(queue.enqueueWriteBuffer(this->materials_buffer, CL_TRUE, 0, sizeof(CL_Material) * material_count, materials),
		"Can't write materials buffer");

	delete materials;
}

// Creates new CL_General object containing all run-time modifiable variables and writes to GPU general_buffer
void RayGPU::loadGeneral() {

	CL_General general;
	general.triangle_count = scene->primitive_count;
	general.lights_count = (int)scene->lights.size();

	general.cast_shadows = scene->cast_shadows;
	general.ambient_light = scene->ambient_light_coeff;
	general.supersampling = scene->supersampling;
	general.max_reflections = scene->max_reflections;
	general.background_color = vec3tofloat3(vec3(scene->background_color[0], scene->background_color[1], scene->background_color[2]));
	general.epsilon = scene->epsilon;

	checkError(queue.enqueueWriteBuffer(this->general_buffer, CL_TRUE, 0, sizeof(CL_General), &general), "Can't write general buffer");
}

// Creates new CL_Camera object containing camera data and write to GPU camera_buffer
void RayGPU::loadCamera() {

	CL_Camera camera = {
		vec3tofloat3(scene->camera.pos),
		vec3tofloat3(scene->camera.dir),
		vec3tofloat3(scene->camera.right),
		vec3tofloat3(scene->camera.up),
		(int)scene->camera.width,
		(int)scene->camera.height,
		scene->camera.aspect_ratio
	};

	checkError(queue.enqueueWriteBuffer(this->camera_buffer, CL_TRUE, 0, sizeof(CL_Camera), &camera), "Can't write camera buffer");
}

// Fills lights_buffer with CL_PointLight objects
void RayGPU::loadLights() {

	size_t lights_count = scene->lights.size();
	CL_PointLight* lights_data = new CL_PointLight[lights_count];

	for (int i = 0; i < lights_count; i++) {
		lights_data[i].position = vec3tofloat3(scene->lights[i]->location);
		lights_data[i].color = vec3tofloat3(vec3(1.0f));
		lights_data[i].intensity = scene->lights[i]->intensity;
	}

	checkError(queue.enqueueWriteBuffer(this->lights_buffer, CL_TRUE, 0, sizeof(CL_PointLight) * lights_count, lights_data),
		"Can't write lights buffer");

	delete lights_data;
}

// Fills triangles_buffer with triangles
void RayGPU::loadTriangles() {

	int triangles_count = (int)scene->objects.size();
	CL_Triangle* triangles_data = new CL_Triangle[triangles_count];

	for (int i = 0; i < triangles_count; i++) {
		Triangle* current = (Triangle*)scene->objects[i];

		triangles_data[i].vertices[0] = vec3tofloat3(current->a);
		triangles_data[i].vertices[1] = vec3tofloat3(current->b);
		triangles_data[i].vertices[2] = vec3tofloat3(current->c);

		triangles_data[i].vert_normals[0] = vec3tofloat3(current->a_norm);
		triangles_data[i].vert_normals[1] = vec3tofloat3(current->b_norm);
		triangles_data[i].vert_normals[2] = vec3tofloat3(current->c_norm);

		triangles_data[i].normal = vec3tofloat3(current->face_normal);

		triangles_data[i].mat_id = current->matId;
	}

	checkError(queue.enqueueWriteBuffer(this->triangles_buffer, CL_TRUE, 0, sizeof(CL_Triangle) * triangles_count, triangles_data),
		"Can't write triangles buffer");

	delete triangles_data;
}

// Fills bvh_buffer with currently selected accelerator flattened bvh nodes
void RayGPU::loadBVH() {

	uint node_count = this->accelerators[this->current_accelerator_id]->getNodeCount();

	//this->bvh_buffer = cl::Buffer(this->context, CL_MEM_READ_ONLY, sizeof(CL_BVHAccelArrayNode) * node_count);
	//checkError(this->raytracing_kernel.setArg(4, this->bvh_buffer), "Can't set BVH buffer");


	// TODO proper general BVH, not just LuxBVH
	BVHAccelArrayNode* bvh_root = ((LuxBVH*)(this->accelerators[this->current_accelerator_id]))->nodes;

	CL_BVHAccelArrayNode* bvh_data = new CL_BVHAccelArrayNode[node_count];

	for (uint i = 0; i < node_count; i++) {
		bvh_data[i].bounds[0] = vec3tofloat3(bvh_root[i].bbox.min);
		bvh_data[i].bounds[1] = vec3tofloat3(bvh_root[i].bbox.max);

		bvh_data[i].primitive = bvh_root[i].primitive;
		bvh_data[i].skipIndex = bvh_root[i].skipIndex;
	}

	queue.enqueueWriteBuffer(this->bvh_buffer, CL_TRUE, 0, sizeof(CL_BVHAccelArrayNode) * node_count, bvh_data);

	delete bvh_data;
}

/// Constructor and destructor
RayGPU::RayGPU(Scene* scene) {
	this->scene = scene;
	this->name = "OpenCL Ray Tracing";
	cout << "[Renderer] Initializing " << this->name << " renderer.\n";

	// Building renderers acceleration structures
	this->accelerators.push_back(new LuxBVH(scene)); // Luxrays BVH builder

	// Setting default acceleration structure
	this->current_accelerator_id = 0;

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

	// Querying platforms
	cl::Platform::get(&platforms);
	if (platforms.size() == 0) {
		std::cout << "[RayGPU] No OpenCL platforms found. No GPGPU acceleration will be possible!\n";
		return;
	}

#if _DEBUG
	for (size_t i = 0; i < platforms.size(); i++) {
		vector<cl::Device> devices;
		platforms[i].getDevices(CL_DEVICE_TYPE_ALL, &devices);
		cout << "\n\nFound platform: " << platforms[i].getInfo<CL_PLATFORM_NAME>() << " by " << platforms[i].getInfo<CL_PLATFORM_VENDOR>() << std::endl;
		cout << "It has devices:\n";
		for (size_t j = 0; j < devices.size(); j++) {
			cout << "Device #" << j << ": " << devices[j].getInfo<CL_DEVICE_NAME>() << endl;
		}
	}
#endif

	// TODO auto selection of device
	vector<cl::Device> devices_;
	default_platform = platforms[0];
	default_platform.getDevices(CL_DEVICE_TYPE_ALL, &devices_);
	default_device = devices_[0];

	//
	//	// If more than one platform - figure out which platform to use
	//	if (all_platforms.size() > 1) {
	//
	//		// checkErroring all platforms if there are any dedicated GPUs
	//		for (int it = 0; it < all_platforms.size(); it++) {
	//
	//			vector<cl::Device> devices;
	//			all_platforms[it].getDevices(CL_DEVICE_TYPE_GPU | CL_DEVICE_TYPE_CPU, &devices);
	//
	//#ifdef _DEBUG
	//			// reporting about found platform
	//			std::cout << "[OpenCL] Found platform " << all_platforms[it].getInfo<CL_PLATFORM_NAME>() << " by "
	//				<< all_platforms[it].getInfo<CL_PLATFORM_VENDOR>() << std::endl;
	//			// reporting about devices on the platform
	//			for (int i = 0; i < devices.size(); i++) {
	//				cout << "[OpenCL] Device #" << i << ": " << devices[i].getInfo<CL_DEVICE_NAME>() << endl;
	//			}
	//#endif
	//
	//			// If the vendor is Nvidia or AMD, use this platform
	//			for (cl::Device device : devices) {
	//
	//				if (device.getInfo<CL_DEVICE_VENDOR_ID>() == 4318 
	//					|| device.getInfo<CL_DEVICE_VENDOR_ID>() == 4098) {
	//
	//					default_platform = all_platforms[it];
	//					default_device = device;
	//					hasDedicatedGPU = true;
	//					break;
	//				}
	//			}
	//		}
	//	}
	//	else {
	//#ifdef _DEBUG
	//		// reporting about found platform
	//		std::cout << "[OpenCL] Found platform " << all_platforms[0].getInfo<CL_PLATFORM_NAME>() << " by "
	//			<< all_platforms[0].getInfo<CL_PLATFORM_VENDOR>() << std::endl;
	//#endif
	//		default_platform = all_platforms[0];
	//		vector<cl::Device> devices;
	//		default_platform.getDevices(CL_DEVICE_TYPE_GPU | CL_DEVICE_TYPE_CPU, &devices);
	//
	//		for (cl::Device device : devices) {
	//
	//			if (device.getInfo<CL_DEVICE_VENDOR_ID>() == 4318
	//				|| device.getInfo<CL_DEVICE_VENDOR_ID>() == 4098) {
	//				
	//				default_device = device;
	//				hasDedicatedGPU = true;
	//				break;
	//			}
	//		}
	//
	//		if (!hasDedicatedGPU) {
	//			default_device = devices[0];
	//		}
	//
	//	}

	// Filling strings of used platform and device
	used_device_name = default_device.getInfo<CL_DEVICE_NAME>();
	used_platform_name = default_platform.getInfo<CL_PLATFORM_NAME>();
	device = used_device_name;
	platform = used_platform_name;

	// Creating context
	context = cl::Context({ default_device });

	// Creating queue
	queue = cl::CommandQueue(context, default_device, 0, &err);
	checkError(err, "Error creating command queue!");

	// Buffers
	//cout << "[OpenCL] Allocating buffers on device.\n";
	this->general_buffer = cl::Buffer(this->context, CL_MEM_READ_ONLY, sizeof(CL_General));
	this->camera_buffer = cl::Buffer(this->context, CL_MEM_READ_ONLY, sizeof(CL_Camera));
	this->materials_buffer = cl::Buffer(this->context, CL_MEM_READ_ONLY, sizeof(CL_Material) * scene->materials.size());
	this->lights_buffer = cl::Buffer(this->context, CL_MEM_READ_ONLY, sizeof(CL_PointLight) * scene->lights.size());
	this->triangles_buffer = cl::Buffer(this->context, CL_MEM_READ_ONLY, sizeof(CL_Triangle) * scene->primitive_count);
	this->bvh_buffer = cl::Buffer(this->context, CL_MEM_READ_ONLY, sizeof(CL_BVHAccelArrayNode) * this->accelerators[this->current_accelerator_id]->getNodeCount());

#if WAVEFRONT
	this->rays_buffer = cl::Buffer(this->context, CL_MEM_READ_WRITE, sizeof(CL_Ray) * (SCRWIDTH * SCRHEIGHT));
	this->intersections_buffer = cl::Buffer(this->context, CL_MEM_READ_WRITE, sizeof(CL_Intersection) * (SCRWIDTH * SCRHEIGHT));
#endif

	this->screen_buffer = cl::Buffer(context, CL_MEM_WRITE_ONLY, sizeof(int) * SCRHEIGHT * SCRWIDTH);
	//cout << "[OpenCL] Screen buffer size: " << (sizeof(int) * SCRHEIGHT * SCRWIDTH / 1024) << "Kb\n";

	// Load .cl source files
	vector<string> source_files = {
		"structs.cl", // this should be first
		"declarations.cl",
		"intersections.cl",
		"tracing.cl",
		"shading.cl",
		"main.cl" // this one should be last
	};
	this->loadProgram(source_files);

#if WAVEFRONT
	// Wavefront approach
	this->primary_ray_generator_kernel = cl::Kernel(program, "generate_primary_rays");
	checkError(this->primary_ray_generator_kernel.setArg(0, this->general_buffer), "Can't set general buffer");
	checkError(this->primary_ray_generator_kernel.setArg(1, this->camera_buffer), "Can't set camera buffer");
	checkError(this->primary_ray_generator_kernel.setArg(2, this->rays_buffer), "Can't set rays buffer");

	this->intersection_kernel = cl::Kernel(program, "compute_intersections");
	checkError(this->intersection_kernel.setArg(0, this->rays_buffer), "Can't set rays buffer");
	checkError(this->intersection_kernel.setArg(1, this->triangles_buffer), "Can't set triangles buffer");
	checkError(this->intersection_kernel.setArg(2, this->bvh_buffer), "Can't set BVH buffer");
	checkError(this->intersection_kernel.setArg(3, this->camera_buffer), "Can't set camera buffer");
	checkError(this->intersection_kernel.setArg(4, this->intersections_buffer), "Can't set intersection buffer");

	this->shading_kernel = cl::Kernel(program, "compute_shading");
	checkError(this->shading_kernel.setArg(0, this->intersections_buffer), "Can't set intersection buffer");
	checkError(this->shading_kernel.setArg(1, this->materials_buffer), "Can't set materials buffer");
	checkError(this->shading_kernel.setArg(2, this->rays_buffer), "Can't set rays buffer");
	checkError(this->shading_kernel.setArg(3, this->lights_buffer), "Can't set lights buffer");
	checkError(this->shading_kernel.setArg(4, this->general_buffer), "Can't set general buffer");
	checkError(this->shading_kernel.setArg(5, this->camera_buffer), "Can't set camera buffer");
	checkError(this->shading_kernel.setArg(6, this->bvh_buffer), "Can't set BVH buffer");
	checkError(this->shading_kernel.setArg(7, this->triangles_buffer), "Can't set triangles buffer");
	checkError(this->shading_kernel.setArg(8, this->screen_buffer), "Can't set screen buffer");

#else 

	this->raytracing_kernel = cl::Kernel(program, "raytrace");
	checkError(this->raytracing_kernel.setArg(0, this->general_buffer), "Can't set general buffer");
	checkError(this->raytracing_kernel.setArg(1, this->camera_buffer), "Can't set camera buffer");
	checkError(this->raytracing_kernel.setArg(2, this->lights_buffer), "Can't set lights buffer");
	checkError(this->raytracing_kernel.setArg(3, this->triangles_buffer), "Can't set triangles buffer");
	checkError(this->raytracing_kernel.setArg(4, this->bvh_buffer), "Can't set BVH buffer");
	checkError(this->raytracing_kernel.setArg(5, this->screen_buffer), "Can't set screen buffer");
	checkError(this->raytracing_kernel.setArg(6, this->materials_buffer), "Can't set materials buffer");

#endif
	
	// Loading triangles, materials and BVH
	this->loadTriangles();
	this->loadMaterials();
	this->loadBVH();

	cout << "[Renderer] Finished initialization of " << this->name << " renderer.\n\n";
}

RayGPU::~RayGPU() {
	// Freeing up acceleration structures
	for (Accelerator* acc : this->accelerators) {
		delete acc;
	}
	this->accelerators.clear();

	// TODO: research if it is true that C++ wrapper for OpenCL handles all deallocation stuff
}

/// Public render() function
// Traces all primary rays and plots them on screen
void RayGPU::render(Tmpl8::Surface* screen) {

	// frame buffer
	int* frame_buffer = (int*) screen->GetBuffer();

	// load lights, camera and general
	this->loadGeneral();
	this->loadCamera();
	this->loadLights(); // TODO do this all only on change?

#if WAVEFRONT
	// pseudo wavefront implementation, to be honest
	checkError(queue.enqueueNDRangeKernel(this->primary_ray_generator_kernel, cl::NullRange,
		cl::NDRange(SCRWIDTH, SCRHEIGHT), cl::NullRange), "Can't enqueue primary ray generator");

	checkError(queue.enqueueNDRangeKernel(this->intersection_kernel, cl::NullRange,
		cl::NDRange(SCRWIDTH, SCRHEIGHT), cl::NullRange), "Can't enqueue intersection kernel");

	// TODO: compact misses

	checkError(queue.enqueueNDRangeKernel(this->shading_kernel, cl::NullRange,
		cl::NDRange(SCRWIDTH, SCRHEIGHT), cl::NullRange), "Can't enqueue shading kernel");
#else
	// mega kernel
	checkError(queue.enqueueNDRangeKernel(this->raytracing_kernel, cl::NullRange,
		cl::NDRange(SCRWIDTH, SCRHEIGHT), cl::NullRange), "Can't enqueue mega kernel");
#endif

	// load back to host frame buffer
	queue.enqueueReadBuffer(this->screen_buffer, CL_TRUE, 0, sizeof(int) * SCRHEIGHT * SCRWIDTH, frame_buffer);
}
