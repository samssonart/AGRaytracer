#define SCRWIDTH		640
#define SCRHEIGHT		640
#define GPU_RENDER		0

#define GLM_FORCE_RADIANS

#include <inttypes.h>
extern "C" 
{ 
#include "glew.h" 
}
#include "gl.h"
#include "io.h"
#include <ios>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include "SDL.h"
#include "SDL_syswm.h"
#include "SDL_opengl.h"
#include "wglext.h"
#include "template.h"
#include "counters.h"
#include "surface.h"

using namespace std;
using namespace Tmpl8;				// to use template classes
using namespace glm;				// to use glm vector stuff

#include <vector>
#include "fcntl.h"
#include "threads.h"
#include <glm/gtc/matrix_transform.hpp>
#include "freeimage.h"
#include <string>
#include <CL\cl.hpp>
#include "FastBVH.h"
#include "Sphere.h"
#include "Triangle.h"
#include "Flat.h"
#include "Air.h"
#include "Shiny.h"
#include "Glass.h"
#include "game.h"
#include "Mesh.h"
#include "MTL.h"
#include "NoAcceleration.h"

#include "PBRT.h"
#include "LuxBVH.h"

#include "imgui.h"
#include "imgui_impl_sdl.h"

#include "Renderer.h"
#include "RayCPU.h"
#include "RayGPU.h"