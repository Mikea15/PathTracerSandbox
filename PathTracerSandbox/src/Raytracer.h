#pragma once

#include "Ray.h"
#include "Scene.h"

#include <vector>
#include <thread>
#include <atomic>

enum class RenderMode
{
	MultiSampling,
	SubPixelSampling,
	SubPixelSamplingAccumulation,
};

enum class RenderType
{
	Diffuse,
	AmbientOcclusion,
	ShadowRays
};

class Raytracer
{
public:
	struct Settings
	{

	};

	Raytracer() = default;
	Raytracer(int width, int height, int subSamples, int samplesPerPixel, int raysPerFrame);
	~Raytracer();

	void SetSettings(Settings settings) { m_settings = settings; }

	void Trace(const Ray& cameraRay, unsigned int rayCount, unsigned int width, unsigned int height, Vec3* colorData, unsigned char* rgbaData);
	void Clear();

	void GetImage(unsigned char* rgbData) const;

	float GetAvgSamplesPerPixel() const;
	unsigned int GetRayCount() const { return m_rayCount; }
	void SetRenderMode(RenderMode mode);
	void SetRenderType(RenderType type);

	void ToggleLightSampling();

private:
	void TraceSampleAccumulation(const Ray& cameraRay, unsigned int rayCount, unsigned int width, unsigned int height, Vec3* colorData, unsigned char* rgbaData);

	Vec3 TraceRay(const Ray& r, int depth, int E = 1);
	Vec3 TraceShadowRay(const Ray& r);
	Vec3 TraceAO(const Ray& r);
	
private:
	int m_imgWidth;
	int m_imgHeight;
	unsigned int m_pixelCount;

	float m_invWidth;
	float m_invHeight;

	int m_subSamples;
	float m_invSubSamples;

	int m_samplesPerPixel;
	float m_invSamplesPerPixel;

	int m_raysPerFrame;

	Scene scene;

	RenderMode m_renderMode;
	RenderType m_renderType;
	Settings m_settings;

	bool m_enableLightSampling;

	std::vector<Vec3> m_cols;
	std::vector<unsigned int> m_sampleAccumulation;

	std::atomic<unsigned int> m_rayCount;
};