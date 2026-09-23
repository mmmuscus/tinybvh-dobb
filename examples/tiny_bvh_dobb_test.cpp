#define FENSTER_APP_IMPLEMENTATION
#define SCRWIDTH 400 // 1280
#define SCRHEIGHT 400 // 720
#define TILESIZE 20
#include "external/fenster.h" // https://github.com/zserge/fenster

#define TINYBVH_IMPLEMENTATION
#include "tiny_bvh.h"
using namespace tinybvh;

// Other includes
#include <cstdlib>
#include <cstdio>

// Application variables
static BVH baseBvh;
static bvhvec4* tris = 0;
static int triCount = 0, frameIdx = 0, spp = 0;
static bvhvec3 accumulator[SCRWIDTH * SCRHEIGHT];
static std::atomic<int> tileIdx(0);

// ---------------------- RENDERING ----------------------
// Multi-threading
static unsigned threadCount = std::thread::hardware_concurrency();

// Setup view pyramid for a pinhole camera:
// eye, p1 (top-left), p2 (top-right), and p3 (bottom-left)
static bvhvec3 eye(0, 30, 0), p1, p2, p3;
static bvhvec3 view = tinybvh_normalize(bvhvec3(-1, 0, 0));

// Xor32 RNG
static unsigned RandomUInt(unsigned& seed) { seed ^= seed << 13, seed ^= seed >> 17, seed ^= seed << 5; return seed; }
static float RandomFloat(unsigned& seed) { return RandomUInt(seed) * 2.3283064365387e-10f; }

// Ray tracing math
bvhvec3 DiffuseReflection(const bvhvec3 N, unsigned& seed)
{
	bvhvec3 R;
	do
	{
		R = bvhvec3(RandomFloat(seed) * 2 - 1, RandomFloat(seed) * 2 - 1, RandomFloat(seed) * 2 - 1);
	} while (tinybvh_dot(R, R) > 1);
	return tinybvh_normalize(tinybvh_dot(R, N) < 0 ? R : -R);
}
bvhvec3 CosWeightedDiffReflection(const bvhvec3 N, unsigned& seed)
{
	bvhvec3 R = DiffuseReflection(N, seed);
	return tinybvh_normalize(N + R);
}

// Color conversion
bvhvec3 rgb32_to_vec3(const unsigned c)
{
	return bvhvec3((float)(c >> 16), (float)((c >> 8) & 255), (float)(c & 255)) * (1 / 255.f);
}

// Geometry access
bvhvec3 TriangleColor(const unsigned idx) { return rgb32_to_vec3(*(unsigned*)&tris[idx * 3].w); }
bvhvec3 TriangleNormal(const unsigned idx)
{
	bvhvec3 a = tris[idx * 3], b = tris[idx * 3 + 1], c = tris[idx * 3 + 2];
	return tinybvh_normalize(tinybvh_cross(b - a, a - c));
}

// Light transport calculation - Basic recursive Path Tracer with IS and Next Event Estimation
bvhvec3 Trace(BVH& bvh, Ray ray, unsigned& seed, unsigned depth = 0)
{
	// find primary intersection
	bvh.Intersect(ray);
	// shade
	if (ray.hit.t == 1e30f) return bvhvec3(0.6f, 0.7f, 1); // hit nothing
	bvhvec3 I = ray.O + ray.hit.t * ray.D;
	bvhvec3 N = TriangleNormal(ray.hit.prim);
	if (tinybvh_dot(N, ray.D) > 0) N = -N;
	bvhvec3 BRDF = TriangleColor(ray.hit.prim) * (1.0f / 3.14159f);
	bvhvec3 Lpos(RandomFloat(seed) * 30 - 15, 40, RandomFloat(seed) * 6 - 3); // virtual
	float dist = tinybvh_length(Lpos - I);
	bvhvec3 L = (Lpos - I) * (1.0f / dist); // normalize
	bvhvec3 direct = {}, indirect = {};
	float NdotL = tinybvh_dot(N, L), NLdotL = fabs(tinybvh_dot(L, bvhvec3(0, 1, 0)));
	if (NdotL > 0)
		if (!bvh.IsOccluded(Ray(I + L * 0.001f, L, dist)))
			direct = BRDF * NdotL * NLdotL * bvhvec3(9, 9, 8) * 500 * (1.0f / (dist * dist));
	// random bounce
	if (depth < 2)
	{
		bvhvec3 R = CosWeightedDiffReflection(N, seed);
		float pdf = 1.0f / tinybvh_dot(N, R);
		bvhvec3 irradiance = Trace(bvh, Ray(I + R * 0.001f, R), seed, depth + 1);
		indirect = BRDF * irradiance * (1.0f / pdf);
	}
	// finalize
	return direct + indirect;
}

void TraceWorkerThread(uint32_t* buf, float scale, int threadIdx)
{
	const int xtiles = SCRWIDTH / TILESIZE, ytiles = SCRHEIGHT / TILESIZE;
	const int tiles = xtiles * ytiles;
	int tile = threadIdx;
	while (tile < tiles)
	{
		const int tx = tile % xtiles, ty = tile / xtiles;
		unsigned seed = (tile + 17) * 171717 + frameIdx * 1023;
		for (int y = 0; y < TILESIZE; y++) for (int x = 0; x < TILESIZE; x++)
		{
			const int pixel_x = tx * TILESIZE + x, pixel_y = ty * TILESIZE + y;
			const int pixelIdx = pixel_x + pixel_y * SCRWIDTH;
			// setup primary ray
			const float u = (float)pixel_x / SCRWIDTH, v = (float)pixel_y / SCRHEIGHT;
			const bvhvec3 D = tinybvh_normalize(p1 + u * (p2 - p1) + v * (p3 - p1) - eye);
			// trace
			accumulator[pixelIdx] += Trace(baseBvh, Ray(eye, D), seed);
			const bvhvec3 E = accumulator[pixelIdx] * scale;
			// visualize, with a poor man's gamma correct
			const int r = (int)tinybvh_min(255.0f, sqrtf(E.x) * 255.0f);
			const int g = (int)tinybvh_min(255.0f, sqrtf(E.y) * 255.0f);
			const int b = (int)tinybvh_min(255.0f, sqrtf(E.z) * 255.0f);
			buf[pixelIdx] = b + (g << 8) + (r << 16);
		}
		tile = tileIdx++;
	}
}

// ---------------------- LOAD DATA ----------------------

// Scene management - Append a file, with optional position, scale and color override, tinyfied
void AddMesh(const char* file, float scale = 1, bvhvec3 pos = {}, int c = 0, int N = 0)
{
	std::fstream s{ file, s.binary | s.in };
	s.read((char*)&N, 4);
	bvhvec4* data = (bvhvec4*)malloc64((N + triCount) * 48);
	if (tris) memcpy(data, tris, triCount * 48), free64(tris);
	tris = data, s.read((char*)tris + triCount * 48, N * 48), triCount += N;
	for (int* b = (int*)tris + (triCount - N) * 12, i = 0; i < N * 3; i++)
		*(bvhvec3*)b = *(bvhvec3*)b * scale + pos, b[3] = c ? c : b[3], b += 4;
}

// Load meshes we want to test with dobb bvh
void LoadMeshData() 
{
	AddMesh("./testdata/cryteksponza.bin", 1, bvhvec3(0), 0xffffff);
}

// ---------------------- BUILD BVH ----------------------

void BuildBvh() 
{
	baseBvh.Build(tris, triCount);
}

// ---------------------- MISC ----------------------

void Init() {
	LoadMeshData();
	BuildBvh();

	// load camera position / direction from file
	std::fstream t = std::fstream{ "camera.bin", t.binary | t.in };
	if (!t.is_open()) {
		t.read((char*)&eye, sizeof(eye));
		t.read((char*)&view, sizeof(view));
		t.close();
	}

	// Update camera vector once
	bvhvec3 right = tinybvh_normalize(tinybvh_cross(bvhvec3(0, 1, 0), view));
	bvhvec3 up = 0.8f * tinybvh_cross(view, right);
	bvhvec3 C = eye + 1.2f * view;
	p1 = C - right + up;   // top-left
	p2 = C + right + up;   // top-right
	p3 = C - right - up;   // bottom-left
}

void Tick(float delta_time_s, fenster& f, uint32_t* buf)
{
	frameIdx++;
	// render tiles
	const float scale = 1.0f / spp++;
	tileIdx = threadCount;
	std::vector<std::thread> threads;
	for (uint32_t i = 0; i < threadCount; i++)
		threads.emplace_back(&TraceWorkerThread, buf, scale, i);
	for (auto& thread : threads) thread.join();
	// print frame time / rate in window title
	char title[50];
	snprintf(title, sizeof(title), "tiny_bvh %.2f s %.2f Hz", delta_time_s, 1.0f / delta_time_s);
	fenster_update_title(&f, title);
}

// Application Shutdown
void Shutdown()
{
	// save camera position / direction to file
	std::fstream s = std::fstream{ "camera.bin", s.binary | s.out };
	s.write((char*)&eye, sizeof(eye));
	s.write((char*)&view, sizeof(view));
	s.close();
}