#include "precomp.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
#define DEMO 0

// Timer for benchmarks and FPS
static timer time_passed;
static int frames = 0;
static float starttime = 0.0f;
static int frames_threshold = 10;

/// TODOS, prioritized
// OpenCL: refactor device selection, now it's just 0th device on 0th platform
// OpenCL: refraction
// OpenCL: finish the wavefront implementation - add ray compression step
// fix compiler warnings
// textures
// add minT, move maxT into rays instead of holding them in intersection struct
// more BVHs: check out more LuxRays bvh builders, should be doable to add SBVH
// skybox
// bvh visualization?
// factor out interpolated normals as a run time parameter
// Return click debug to RayCPU and maybe try to introduce it to RayGPU
// stratification (jittering at least)
// adaptive supersampling instead of naive supersampling
// Post-pro: tone mapping - http://filmicgames.com/archives/75
// normal maps
// instead of uint put std _t classes
// multiple lights - check correctness with ground truth, not sure it's correct
// packet traversal
// improve AAC BVH building - seems too underwhelming for so much effort - maybe radixsort instead of std::sort?
// Post-pro: rescale filter?

// Initialize the application
void Game::Init() {

#if DEMO 
	/// Demo 1
	// Camera
	scene.camera = Camera(vec3(0, 5, -50), vec3(0, 0, 1));
	// Mesh, 100K triangles
	loadObj("scenes/dragon/dragon.obj");
	// Additional reflective triangle behind
	scene.materials.push_back(new Shiny(vec3(0.0f), 100.0f, 0.8f));
	scene.objects.push_back(new Triangle(vec3(-10, 0, 30), vec3(30, 0, 0), vec3(10, 20, 15), 2));
	// Light
	scene.lights.push_back(new Light(vec3(100, 25, -150), 1.0f));
	// Only OpenCL renderer for the demo
	scene.primitive_count = (int)scene.objects.size();
	scene.renderers.push_back(new RayGPU(&scene));

#else
	// uncomment one of the lines to try it out
	//static string INPUT_FILE = "scenes/cube/cube.obj"; // 12 triangles
	//static string INPUT_FILE = "scenes/teapot-low/teapot-low.obj"; // 240 triangles
	//static string INPUT_FILE = "scenes/teapot/teapot.obj"; // 16K triangles
	static string INPUT_FILE = "scenes/lena/lenaStatue.obj"; // 36K triangles
	//static string INPUT_FILE = "scenes/dragon/dragon.obj"; // 100K triangles

	scene.camera = Camera(vec3(0, 30, -168), vec3(0, 0, 1));
	scene.lights.push_back(new Light(vec3(100, 25, -150), 1.0f));
	loadObj(INPUT_FILE);

	// Renderers, each has a set of acceleration structure they are building on initialization
	scene.primitive_count = (int)scene.objects.size();
	scene.renderers.push_back(new RayGPU(&scene));
	//scene.renderers.push_back(new RayCPU(&scene));
#endif

	// FPS counter initializer
	// Not using ImGui FPS counter, because it's using secondary controls SDL window, might be inaccurate in some cases
	time_passed.reset();
	starttime = time_passed.elapsed();
}

// Converts tinyobjloader's material to my material class
MTL* convertToMTL(tinyobj::material_t mat) {
	return new MTL(mat.ambient, mat.diffuse, mat.specular, mat.transmittance, mat.emission, mat.shininess, mat.ior,
		mat.dissolve, mat.illum);
}

// Loads .obj from given string and pushes the primitives into scene.objects
void Game::loadObj(string input_file) {
	cout << "[OBJ] Loading " << input_file << ":\n";

	vector<tinyobj::shape_t> shapes;
	vector<tinyobj::material_t> materials;

	string err;

	// computing basepath, needed for materials
	string input_file_basepath = input_file.substr(0, input_file.find_last_of('/') + 1);
	tinyobj::LoadObj(shapes, materials, err, input_file.c_str(), input_file_basepath.c_str(), true);

	if (!err.empty()) {
		cerr << "[OBJ] Error loading the file: " << err << endl;
	}

	for (size_t i = 0; i < materials.size(); i++) {
		scene.materials.push_back(convertToMTL(materials[i]));
	}

	for (size_t i = 0; i < shapes.size(); i++) {
		scene.meshes.push_back(new Mesh(&scene.objects, shapes[i].mesh.positions, shapes[i].mesh.normals, shapes[i].mesh.texcoords,
			shapes[i].mesh.indices, shapes[i].mesh.num_vertices, shapes[i].mesh.material_ids, shapes[i].name));
	}
	
#ifdef _DEBUG
	for (size_t i = 0; i < shapes.size(); i++) {
		printf("shape[%ld].name = %s\n", i, shapes[i].name.c_str());
		printf("Size of shape[%ld].indices: %ld\n", i, shapes[i].mesh.indices.size());
		printf("Size of shape[%ld].material_ids: %ld\n", i, shapes[i].mesh.material_ids.size());
		assert((shapes[i].mesh.indices.size() % 3) == 0);
		for (size_t f = 0; f < shapes[i].mesh.indices.size() / 3; f++) {
			printf("  idx[%ld] = %d, %d, %d. mat_id = %d\n", f, shapes[i].mesh.indices[3 * f + 0], shapes[i].mesh.indices[3 * f + 1], shapes[i].mesh.indices[3 * f + 2], shapes[i].mesh.material_ids[f]);
		}

		printf("shape[%ld].vertices: %ld\n", i, shapes[i].mesh.positions.size());
		assert((shapes[i].mesh.positions.size() % 3) == 0);
		for (size_t v = 0; v < shapes[i].mesh.positions.size() / 3; v++) {
			printf("  v[%ld] = (%f, %f, %f)\n", v,
				shapes[i].mesh.positions[3 * v + 0],
				shapes[i].mesh.positions[3 * v + 1],
				shapes[i].mesh.positions[3 * v + 2]);
		}
	}
#endif

	cout << "[OBJ] Finished loading " << input_file << ".\n\n";

}


// Main application tick function
void Game::Tick(float _DT) {
	screen->Clear(0);
	scene.camera.update();

	scene.renderers[scene.current_renderer_id]->render(screen);

	// fps
	frames++;
	float time_passed_f = time_passed.elapsed();
	if (time_passed_f - starttime > 0.25f && frames > frames_threshold) {
		scene.calculated_fps = (float)frames / (time_passed_f - starttime) * 1000;
		starttime = time_passed_f;
		frames = 0;
		frames_threshold = (int)scene.calculated_fps;
	}
}

// Camera and view controls
void Game::KeyDown(int _Key) {

#ifdef _DEBUG
	cout << "Pressed button with ID: " << _Key << endl;
#endif

	// movement
	if (_Key == 26) { // w
		scene.camera.moveForward(1.0f);
	}
	else if (_Key == 4) { // a
		scene.camera.moveLeft(1.0f);
	}
	else if (_Key == 22) { // s
		scene.camera.moveBackward(1.0f);
	}
	else if (_Key == 7) { // d
		scene.camera.moveRight(1.0f);
	}
	else if (_Key == 8) { // q
		scene.camera.moveUp(1.0f);
	}
	else if (_Key == 20) { // e
		scene.camera.moveDown(1.0f);
	}

	// looking
	else if (_Key == 82) { // arrow up
		scene.camera.lookUp(1.0f);
	}
	else if (_Key == 80) { // arrow left
		scene.camera.lookLeft(1.0f);
	}
	else if (_Key == 81) { // arrow down
		scene.camera.lookDown(1.0f);
	}
	else if (_Key == 79) { // arrow right
		scene.camera.lookRight(1.0f);
	}

	// various
	else if (_Key == 5) { // b - next accelerator
		scene.renderers[scene.current_renderer_id]->current_accelerator_id++;
		if (scene.renderers[scene.current_renderer_id]->current_accelerator_id >= scene.renderers[scene.current_renderer_id]->accelerators.size())
			scene.renderers[scene.current_renderer_id]->current_accelerator_id = 1;
	}
	else if (_Key == 25) { // v - prev accelerator
		scene.renderers[scene.current_renderer_id]->current_accelerator_id--;
		if (scene.renderers[scene.current_renderer_id]->current_accelerator_id <= 0)
			scene.renderers[scene.current_renderer_id]->current_accelerator_id = scene.renderers[scene.current_renderer_id]->accelerators.size() - 1;
	}
	else if (_Key == 10) { // g - next renderer
		scene.current_renderer_id++;
		if (scene.current_renderer_id >= scene.renderers.size())
			scene.current_renderer_id = 0;
	}
	else if (_Key == 9) { // f - prev renderer
		scene.current_renderer_id--;
		if (scene.current_renderer_id >= scene.renderers.size())
			scene.current_renderer_id = scene.renderers.size() - 1;
	}
	else if (_Key == 17) { // n - toggle shadows
		scene.cast_shadows = !scene.cast_shadows;
	}
	else if (_Key == 21) { // r - reset camera
		scene.camera.reset();
	}
	else if (_Key == 15) { // l - move light to camera position
		scene.lights[0]->location = vec3(scene.camera.pos.x, scene.camera.pos.y, scene.camera.pos.z);
	}

}

// Close down application
void Game::Shutdown() {

	// destroying all scene primitives
	for (Primitive* primitive : scene.objects) {
		//delete primitive->material;
		delete primitive;
	}
	scene.objects.clear();

	// destroying all materials, including starting air medium
	for (Material* mat : scene.materials) {
		delete mat;
	}
	scene.materials.clear();

	// destroying renderers
	for (Renderer* ren : scene.renderers) {
		delete ren;
	}
	scene.renderers.clear();

	// destroying all scene lights
	for (Light* light : scene.lights) {
		delete light;
	}
	scene.lights.clear();

}

