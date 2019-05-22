// PathTracerSandbox
// Release x64
// 7.1		MRays/s
// 6.89		MRays/s / SIMD
// 10.91	MRays/s / OpenMP
// 15.07	MRays/s / SIMD / OpenMP

#include <Windows.h>

#include <chrono>
#include <iostream>
#include <string>
#include <iomanip> // setprecision
#include <sstream> // stringstream

#include "imgui/imgui.h"
#include "imgui/imgui_impl_sdl.h"
#include "imgui/imgui_impl_opengl3.h"

#include <GL/glew.h>

#include "SDL.h"
#include "Vec3.h"
#include "Ray.h"
#include "Raytracer.h"

static const unsigned int SCREEN_WIDTH = 800u;
static const unsigned int SCREEN_HEIGHT = 600u;
static const unsigned int PIXEL_COUNT = SCREEN_WIDTH * SCREEN_HEIGHT;

Vec3 camPos(50, 45, 290.6);
Raytracer tracer(SCREEN_WIDTH, SCREEN_HEIGHT, 2, 6, 4 * 1024, 128, 15.0f);

// main entry point
int main(int argc, char* argv[])
{
	// initialize SDL
	if (SDL_Init(SDL_INIT_VIDEO) != 0)
	{
		return 1;
	}

	SDL_WindowFlags window_flags = (SDL_WindowFlags)(SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);
	SDL_Window* window = SDL_CreateWindow("PathTracer Sandbox", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, SCREEN_WIDTH, SCREEN_HEIGHT, window_flags);
	// SDL_GLContext gl_context = SDL_GL_CreateContext(window);

	// SDL_Window* window = SDL_CreateWindow("PathTracer Sandbox", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, SCREEN_WIDTH, SCREEN_HEIGHT, 0);
	SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
	SDL_Texture* texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, SCREEN_WIDTH, SCREEN_HEIGHT);

	glewInit();

	float timeMs = 0.0f;
	float raysPerSecMax = 0.0f;

	std::stringstream stream;

	bool isDone = false;
	// run the main loop
	while (!isDone)
	{ 
		stream.str("");
		stream.clear();

		// poll SDL events, listen for keyboard event
		SDL_Event event = {};
		while (SDL_PollEvent(&event))
		{
			switch (event.type)
			{
			case SDL_KEYDOWN:
				break;
			case SDL_KEYUP:
				switch (event.key.keysym.sym)
				{
				case SDLK_F1:
					// TODO:
					// SaveToFile(SCREEN_WIDTH, SCREEN_HEIGHT, rgbaData, PIXEL_COUNT * 4);
					break;
				case SDLK_F2:
					tracer.SetRenderType(RenderType::Diffuse);
					break;
				case SDLK_F3:
					tracer.SetRenderType(RenderType::AmbientOcclusion);
					break;
				case SDLK_F4:
					tracer.SetRenderType(RenderType::ShadowRays);
					break;
				case SDLK_F5:
					tracer.SetRenderMode(RenderMode::SubPixelSamplingAccumulation);
					break;
				case SDLK_F6:
					tracer.SetRenderMode(RenderMode::SubPixelSampling);
					break;
				case SDLK_F7:
					tracer.SetRenderMode(RenderMode::MultiSampling);
					break;
				case SDLK_c:
					tracer.Clear();
					break;
				}
				break;
			case SDL_QUIT:
				isDone = true;
				break;
			default: break;
			}
		}


		Ray cameraRay(
			camPos,
#if USE_OPTIMIZED_VEC 
			normalize_vector(Vec3(0, 0, -1))
#else
			Vec3(0, 0, -1).normalize()
#endif
		);

		auto frameStartTime = std::chrono::high_resolution_clock::now();

		tracer.Trace(cameraRay);

		auto timeSpan = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - frameStartTime);

		float frameTimeMs = static_cast<float>(timeSpan.count());

		unsigned int rayCount = tracer.GetRayCount();
		const float agvSpp = tracer.GetAvgSamplesPerPixel();

		timeMs += frameTimeMs;
		const float frameTimeS = frameTimeMs * 0.001f;
		const float raysPerSec = static_cast<float>(rayCount / frameTimeS) * 0.000001f;
		if (raysPerSec > raysPerSecMax)
		{
			raysPerSecMax = raysPerSec;
		}

		stream << std::fixed << std::setprecision(2) << "total time: " << timeMs * 0.001f << " s";
		stream << std::fixed << std::setprecision(2) << " - render time: " << frameTimeS << " s";
		stream << std::fixed << std::setprecision(2) << " - M Rays/s: " << raysPerSec;
		stream << std::fixed << std::setprecision(2) << " - Max M Rays/s: " << raysPerSecMax;
		stream << std::fixed << std::setprecision(2) << " - Avg SPP: " << agvSpp;

		std::string s = stream.str();

		std::cout << s << std::endl;

		SDL_SetWindowTitle(window, s.c_str());

		// begin render.
		// SDL_GL_MakeCurrent(window, gl_context);

		SDL_UpdateTexture(texture, NULL, tracer.GetRGBData().data(), SCREEN_WIDTH * 4);
		SDL_RenderCopy(renderer, texture, NULL, NULL);
		
		SDL_RenderPresent(renderer);

		// end render.
		// SDL_GL_SwapWindow(window);
	}

	// clean up SDL
	// SDL_GL_DeleteContext(gl_context);

	SDL_DestroyTexture(texture);
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();

	return 0;
}
