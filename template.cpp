// Template, major revision 7
// IGAD/NHTV - Jacco Bikker - 2006-2015

// Note:
// this version of the template uses SDL2 for all frame buffer interaction
// see: https://www.libsdl.org
// for glm (OpenGL mathematics) see http://glm.g-truc.net

#ifdef _MSC_VER
#pragma warning (disable : 4530) // complaint about exception handler
#pragma warning (disable : 4273)
#pragma warning (disable : 4311) // pointer truncation from HANDLE to long
#endif

#include "precomp.h"

namespace Tmpl8 { 

double timer::inv_freq = 1;

void NotifyUser( char* s )
{
	HWND hApp = FindWindow( NULL, "Template" );
	MessageBox( hApp, s, "ERROR", MB_OK );
	exit( 0 );
}
}

using namespace Tmpl8;
using namespace std;

int ACTWIDTH, ACTHEIGHT;
static bool firstframe = true;

Surface* surface = 0;
Game* game = 0;
SDL_Window* window = 0;

#ifdef _MSC_VER
void redirectIO()
{
	static const WORD MAX_CONSOLE_LINES = 500;
	int hConHandle;
	HANDLE lStdHandle;
	CONSOLE_SCREEN_BUFFER_INFO coninfo;
	FILE *fp;
	AllocConsole();
	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE),
	&coninfo);
	coninfo.dwSize.Y = MAX_CONSOLE_LINES;
	SetConsoleScreenBufferSize(GetStdHandle(STD_OUTPUT_HANDLE),
	coninfo.dwSize);
	lStdHandle = GetStdHandle(STD_OUTPUT_HANDLE);
	hConHandle = _open_osfhandle((intptr_t)lStdHandle, _O_TEXT);
	fp = _fdopen( hConHandle, "w" );
	*stdout = *fp;
	setvbuf( stdout, NULL, _IONBF, 0 );
	lStdHandle = GetStdHandle(STD_INPUT_HANDLE);
	hConHandle = _open_osfhandle((intptr_t)lStdHandle, _O_TEXT);
	fp = _fdopen( hConHandle, "r" );
	*stdin = *fp;
	setvbuf( stdin, NULL, _IONBF, 0 );
	lStdHandle = GetStdHandle(STD_ERROR_HANDLE);
	hConHandle = _open_osfhandle((intptr_t)lStdHandle, _O_TEXT);
	fp = _fdopen( hConHandle, "w" );
	*stderr = *fp;
	setvbuf( stderr, NULL, _IONBF, 0 );
	ios::sync_with_stdio();
}
#endif

int main( int argc, char **argv ) {  

#ifdef _MSC_VER
	redirectIO();
#endif

	// redirecting IO back to stdout / stderr
	// don't want to touch redirectIO()
	freopen("CON", "w", stdout);
	freopen("CON", "w", stderr);

	// initializing SDL
	if (SDL_Init(SDL_INIT_VIDEO) != 0)
	{
		printf("[Core] SDL error: %s\n", SDL_GetError());
		return EXIT_FAILURE;
	}

	// pixel buffer for ray tracing operations
	surface = new Surface(SCRWIDTH, SCRHEIGHT);
	surface->Clear(0);
	surface->InitCharset();

	// ray tracing stuff
	game = new Game();
	game->SetTarget(surface);

	// main window
	window = SDL_CreateWindow("Ray tracer", 50, 50, SCRWIDTH, SCRHEIGHT, SDL_WINDOW_SHOWN);
	SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
	SDL_Texture* tracingBuffer = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, SCRWIDTH, SCRHEIGHT);

	// UI window
	int w_x, w_y;
	SDL_GetWindowPosition(window, &w_x, &w_y);
	int UI_WINDOW_WIDTH = 250;
	SDL_Window* ui_window = SDL_CreateWindow("UI", w_x + SCRWIDTH, w_y, UI_WINDOW_WIDTH, SCRHEIGHT, SDL_WINDOW_SHOWN | SDL_WINDOW_OPENGL);
	SDL_Renderer* ui_renderer = SDL_CreateRenderer(ui_window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
	SDL_GLContext ui_glContext = SDL_GL_CreateContext(ui_window);
	ImGui_ImplSdl_Init(ui_window);
	
	timer t; t.reset();

	string accelerator_names, renderers_names;
	bool exit_requested = false;
	while (!exit_requested) {

		// Event handling
		SDL_Event event;
		while (SDL_PollEvent(&event)) {

			// for both windows
			if (event.type == SDL_QUIT) {
				exit_requested = true;
				break;
			}
			else if (event.type == SDL_KEYDOWN) {
				if (event.key.keysym.sym == SDLK_ESCAPE) {
					exit_requested = true;
				}
			}

			// for ui window
			if (event.window.windowID == 2) {
				ImGui_ImplSdl_ProcessEvent(&event); // delegate events to imgui sdl bindings
			}
			// for main window
			else { 
				switch (event.type) {
				case SDL_KEYDOWN:
					game->KeyDown(event.key.keysym.scancode);
					break;
				case SDL_KEYUP:
					game->KeyUp(event.key.keysym.scancode);
					break;
				case SDL_MOUSEMOTION:
					game->MouseMove(event.motion.xrel, event.motion.yrel);
					break;
				case SDL_MOUSEBUTTONUP:
					game->MouseUp(event.button.button);
					break;
				case SDL_MOUSEBUTTONDOWN:
					game->MouseDown(event.button.button, event.button.x, event.button.y);
					break;
				default:
					break;
				}
			}
		}

		// Ticking UI stuff
		ImGui_ImplSdl_NewFrame(ui_window);

		// Ticking ray tracing stuff
		if (firstframe) {
			game->Init();
			firstframe = false;

			// build accelerator names string for ui
			for (Renderer* a : game->scene.renderers) {
				renderers_names += a->name;
				renderers_names.push_back('\0');
			}
			renderers_names.push_back('\0');
		}
		game->Tick(t.elapsed());
		t.reset();

		// Rendering primary window
		void* target = 0;
		int pitch;
		SDL_LockTexture(tracingBuffer, NULL, &target, &pitch);
		if (pitch == (surface->GetWidth() * 4)) {
			memcpy(target, surface->GetBuffer(), SCRWIDTH * SCRHEIGHT * 4);
		} else {
			unsigned char* t = (unsigned char*)target;
			for (int i = 0; i < SCRHEIGHT; i++) {
				memcpy(t, surface->GetBuffer() + i * SCRWIDTH, SCRWIDTH * 4);
				t += pitch;
			}
		}
		SDL_UnlockTexture(tracingBuffer);
		SDL_RenderCopy(renderer, tracingBuffer, NULL, NULL);
		SDL_RenderPresent(renderer);

		// Rendering UI window
		ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f));
		ImGui::SetNextWindowSize(ImVec2((float)UI_WINDOW_WIDTH, SCRHEIGHT));
		ImGui::Begin("Controls", (bool*)0, 
			ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
			ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings);

		ImGui::Text("FPS: %.2f (%.2f ms/frame)", game->scene.calculated_fps, 1000 / game->scene.calculated_fps);
//#if _DEBUG
//		ImGui::Text("Intersections: %u", game->scene.accelerators[game->scene.current_accelerator_id]->isecs);
//#endif
		//ImGui::Text("Loaded .obj: %s", INPUT_FILE);
		ImGui::Text("Primitives: %i (%i meshes)", game->scene.primitive_count, game->scene.meshes.size());
		ImGui::Text("Lights: %i", game->scene.lights.size());
		ImGui::Text("Camera pos: (%.3f, %.3f, %.3f)", game->scene.camera.pos.x, game->scene.camera.pos.y, game->scene.camera.pos.z);
		ImGui::Text("Camera dir: (%.3f, %.3f, %.3f)", game->scene.camera.dir.x, game->scene.camera.dir.y, game->scene.camera.dir.z);
		if (ImGui::Button("Reset camera")) {
			game->scene.camera.reset();
		}

		// TODO Ray count?

		ImGui::Separator();
		ImGui::PushItemWidth(120.0f);
		ImGui::Combo("Renderer", (int*)&game->scene.current_renderer_id, renderers_names.c_str());
		if (!game->scene.renderers[game->scene.current_renderer_id]->device.empty()) {
			ImGui::Text("Device: %s", game->scene.renderers[game->scene.current_renderer_id]->device.c_str());
		}
		if (!game->scene.renderers[game->scene.current_renderer_id]->platform.empty()) {
			ImGui::Text("Platform: %s", game->scene.renderers[game->scene.current_renderer_id]->platform.c_str());
		}
		ImGui::Combo("Acceleration", (int*)&game->scene.renderers[game->scene.current_renderer_id]->current_accelerator_id, game->scene.renderers[game->scene.current_renderer_id]->accelerators_ui_string.c_str());
		ImGui::SliderFloat("Ambient light", &game->scene.ambient_light_coeff, 0.0f, 1.0f);
		ImGui::SliderInt("Supersampling", &game->scene.supersampling, 1, 4);
		ImGui::SliderInt("Max reflections", &game->scene.max_reflections, 0, 3);
		ImGui::SliderFloat3("Background color", &game->scene.background_color[0], 0.0f, 1.0f);
		ImGui::SliderFloat("Epsilon", &game->scene.epsilon, 0.001f, 0.05f);
		ImGui::Checkbox("Shadows", &game->scene.cast_shadows);

		// TODO click debugging goes here
		// ImGui::Separator();

		ImGui::Separator();

		ImGui::Text("Keybinds (main window)");
		ImGui::BulletText("Camera position: W/A/S/D/Q/E");
		ImGui::BulletText("Camera direction: arrow keys");
		ImGui::BulletText("Reset camera: R");
		ImGui::BulletText("Toggle shadows: N");
		ImGui::BulletText("Next renderer: G");
		ImGui::BulletText("Prev renderer: F");
		ImGui::BulletText("Next accelerator: B");
		ImGui::BulletText("Prev accelerator: V");
		ImGui::BulletText("Move 0th light to camera: L");

		ImGui::End();

		glViewport(0, 0, (int)ImGui::GetIO().DisplaySize.x, (int)ImGui::GetIO().DisplaySize.y);
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		ImGui::Render();
		SDL_GL_SwapWindow(ui_window);
	}

	game->Shutdown();

	ImGui_ImplSdl_Shutdown();

	SDL_DestroyWindow(window);
	SDL_DestroyWindow(ui_window);
	SDL_GL_DeleteContext(ui_glContext);
	SDL_Quit();

	return 1;
}