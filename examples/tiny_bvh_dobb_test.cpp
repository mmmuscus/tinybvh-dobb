#define FENSTER_APP_IMPLEMENTATION
#define SCRWIDTH 400 // 1280
#define SCRHEIGHT 400 // 720
#define TILESIZE 20
#include "external/fenster.h" // https://github.com/zserge/fenster

#define TINYBVH_IMPLEMENTATION
#include "tiny_bvh.h"
#define TINYBVH_DOBB_IMPLEMENTATION
#include "tiny_bvh_dobb.h"
using namespace tinybvh;

// Other includes
#include <cstdlib>
#include <cstdio>

#ifdef _WIN32
#include <windows.h>
#endif
#include <iostream>
#include <iomanip>

enum Phase {
	BUILD,
	RENDERMBVH,
	RENDERDOBB,
	END
};

// Helpers for experimental setup
static Phase phase = BUILD;
// TODO: organize into single array of structures
static const char* meshLocations[2] = {
	"./testdata/cryteksponza.bin",
	"./testdata/dragon.bin"
};
static const char* meshNames[2] = {
	"CRYTEKSPONZA", "DRAGON"
};
static const int meshCount = sizeof(meshLocations) / sizeof(meshLocations[0]);
static int meshIdx = 0;

// For instrumentation for experiment
struct TraversalStats {
	uint64_t rayCount = 0;
	uint64_t nodeVisits = 0;
	uint64_t boxTests = 0;
	uint64_t triTests = 0;

	void add(const TraversalStats& other) {
		rayCount += other.rayCount;
		nodeVisits += other.nodeVisits;
		boxTests += other.boxTests;
		triTests += other.triTests;
	}

	void reset() {
		rayCount = 0; nodeVisits = 0; boxTests = 0; triTests = 0;
	}

	void print(const char* statName) {
		// Guards against division by 0
		const double guardedRayCount = rayCount ? (double)rayCount : 1.0f;
		std::cout << std::fixed << std::setprecision(3) <<
			statName << rayCount << std::endl <<
			"Avg nodes visited per ray: " << nodeVisits / guardedRayCount << std::endl <<
			"Avg boxes intersected per ray: " << (boxTests / guardedRayCount) << std::endl <<
			"Avg triangles intersected per ray: " << triTests / guardedRayCount << std::endl;
	}
};

// Application variables
static TraversalStats hitStats, shadowStats;
static const int targetSpp = 25;
static const float nodeVisitScale = 1.0f / 250.0f;
static const float triTestScale = 1.0f / 1000.0f;
static BVH_DOBB dobbBvh;
static bvhvec4* tris = 0;
static int triCount = 0, spp = 0;
static bvhvec3 accumulator[SCRWIDTH * SCRHEIGHT];
static std::atomic<int> tileIdx(0);

// ---------------------- INSTRUMENTATION ----------------------

// Same signature and behaviour as BVHBase::IntersectTri (Möller–Trumbore, no opacity maps)
static void IntersectTri(Ray& ray, const uint32_t triIdx, const bvhvec4slice& verts,
	const uint32_t i0, const uint32_t i1, const uint32_t i2)
{
	const bvhvec3 vertex0 = bvhvec3(verts[i0]);
	const bvhvec3 edge1 = bvhvec3(verts[i1]) - vertex0, edge2 = bvhvec3(verts[i2]) - vertex0;
	const bvhvec3 directionCrossEdge2 = tinybvh_cross(ray.D, edge2);
	const float determinant = tinybvh_dot(edge1, directionCrossEdge2);
	if (determinant == 0.0f) return;
	const float inverseDeterminant = 1.0f / determinant;
	const bvhvec3 originOffset = ray.O - vertex0;
	const float baryU = inverseDeterminant * tinybvh_dot(originOffset, directionCrossEdge2);
	const bvhvec3 offsetCrossEdge1 = tinybvh_cross(originOffset, edge1);
	const float baryV = inverseDeterminant * tinybvh_dot(ray.D, offsetCrossEdge1);
	if (!(baryU >= 0.0f && baryV >= 0.0f && baryU + baryV <= 1.0f)) return;
	const float distance = inverseDeterminant * tinybvh_dot(edge2, offsetCrossEdge1);
	if (!(distance >= 0.0f && distance <= ray.hit.t)) return;
	ray.hit.t = distance, ray.hit.u = baryU, ray.hit.v = baryV;
	ray.SetHitPrim(triIdx);
}

// Same signature and behaviour as BVHBase::TriOccludes
static bool TriOccludes(const Ray& ray, const bvhvec4slice& verts, const uint32_t triIdx,
	const uint32_t i0, const uint32_t i1, const uint32_t i2)
{
	(void)triIdx; // only needed for opacity maps
	const bvhvec3 vertex0 = bvhvec3(verts[i0]);
	const bvhvec3 edge1 = bvhvec3(verts[i1]) - vertex0, edge2 = bvhvec3(verts[i2]) - vertex0;
	const bvhvec3 directionCrossEdge2 = tinybvh_cross(ray.D, edge2);
	const float determinant = tinybvh_dot(edge1, directionCrossEdge2);
	if (determinant == 0.0f) return false;
	const float inverseDeterminant = 1.0f / determinant;
	const bvhvec3 originOffset = ray.O - vertex0;
	const float baryU = inverseDeterminant * tinybvh_dot(originOffset, directionCrossEdge2);
	const bvhvec3 offsetCrossEdge1 = tinybvh_cross(originOffset, edge1);
	const float baryV = inverseDeterminant * tinybvh_dot(ray.D, offsetCrossEdge1);
	if (!(baryU >= 0.0f && baryV >= 0.0f && baryU + baryV <= 1.0f)) return false;
	const float distance = inverseDeterminant * tinybvh_dot(edge2, offsetCrossEdge1);
	return distance >= 0.0f && distance <= ray.hit.t;
}

// Copy of DOBB Intersect with instrumentation
void IntersectionInstrumented(const BVH_DOBB& bvh, Ray& ray, TraversalStats* stats = nullptr) {
	if (stats) stats->rayCount++;

	uint32_t nodeStack[512]; uint32_t stackIdx = 0;
	nodeStack[stackIdx++] = 0;

	while (stackIdx > 0) {
		if (stats) stats->nodeVisits++;

		const uint32_t currNodeIdx = nodeStack[--stackIdx];
		const BVH_DOBB::BVHNode& currNode = bvh.dobbNode[currNodeIdx];

		Ray localRay;
		localRay.hit.t = ray.hit.t;
		const uint32_t rotationIdx = currNode.Rotation();

		if (rotationIdx != 255) {
			const DOBB_LUT::Mat3f rotation = bvh.getLutEntry(rotationIdx);
			const bvhvec3 rotationX = bvhvec3(rotation.m[0][0], rotation.m[1][0], rotation.m[2][0]);
			const bvhvec3 rotationY = bvhvec3(rotation.m[0][1], rotation.m[1][1], rotation.m[2][1]);
			const bvhvec3 rotationZ = bvhvec3(rotation.m[0][2], rotation.m[1][2], rotation.m[2][2]);

			localRay.O = bvhvec3(tinybvh_dot(rotationX, ray.O), tinybvh_dot(rotationY, ray.O), tinybvh_dot(rotationZ, ray.O));
			localRay.D = bvhvec3(tinybvh_dot(rotationX, ray.D), tinybvh_dot(rotationY, ray.D), tinybvh_dot(rotationZ, ray.D));
		}
		else {
			localRay.O = ray.O; localRay.D = ray.D;
		}
		localRay.rD = tinybvh_rcp(localRay.D);

		float stepSize[3] = {ldexpf(1.0f, currNode.exponent[0]), ldexpf(1.0f, currNode.exponent[1]), ldexpf(1.0f, currNode.exponent[2]),};

		uint32_t hitMask = 0; uint32_t remainingSlots = currNode.OccupiedMask();

		while (remainingSlots) {
			if (stats) stats->boxTests++;

			// Determine which is the next slot to process
			uint32_t slot = (uint32_t)std::countr_zero(remainingSlots);
			remainingSlots &= remainingSlots - 1;

			// Determine lower and higher extents of the BB
			bvhvec3 slotMinPoint = currNode.nodeMinPoint + bvhvec3(currNode.qLo[0][slot] * stepSize[0], currNode.qLo[1][slot] * stepSize[1], currNode.qLo[2][slot] * stepSize[2]);
			bvhvec3 slotMaxPoint = currNode.nodeMinPoint + bvhvec3(currNode.qHi[0][slot] * stepSize[0], currNode.qHi[1][slot] * stepSize[1], currNode.qHi[2][slot] * stepSize[2]);

			float entryDistance = tinybvh_intersect_aabb(localRay, slotMinPoint, slotMaxPoint);
			if (entryDistance < BVH_FAR)
				hitMask |= 1 << slot;
		}

		uint8_t leafHitMask = hitMask & ~currNode.iMask;
		uint8_t internalHitMask = hitMask & currNode.iMask;

		// Process leaf hits
		while (leafHitMask) {
			const uint32_t slot = (uint32_t)std::countr_zero(leafHitMask);
			leafHitMask &= leafHitMask - 1;

			// Intersect triangles
			const uint32_t firstTri = currNode.triBaseIndex + (currNode.meta[slot] & 31);
			const uint32_t triCount = __popc(currNode.meta[slot] >> 5);

			uint32_t i0, i1, i2;
			for (uint32_t triIdx = 0; triIdx < triCount; triIdx++) {
				if (stats) stats->triTests++;

				const uint32_t primIdx = bvh.dobbPrimIdx[firstTri + triIdx];
				GET_PRIM_INDICES_I0_I1_I2(bvh.bvh8.bvh, primIdx);
				IntersectTri(ray, primIdx, bvh.bvh8.bvh.verts, i0, i1, i2);
			}
		}

		// Process internal hits
		while (internalHitMask) {
			const uint32_t slot = (uint32_t)std::countr_zero(internalHitMask);
			internalHitMask &= internalHitMask - 1;

			// Push internal node
			nodeStack[stackIdx++] = currNode.ChildBaseIndex() +
				__popc(currNode.iMask & ((1u << slot) - 1));
		}
	}
}

// Copy of DOBB isOccluded with instrumentation
bool isOccludedInstrumented(const BVH_DOBB& bvh, const Ray& ray, TraversalStats* stats) {
	if (stats) stats->rayCount++;
	
	uint32_t nodeStack[512];
	uint32_t stackIdx = 0;
	nodeStack[stackIdx++] = 0;

	while (stackIdx > 0) {
		if (stats) stats->nodeVisits++;

		const uint32_t currNodeIdx = nodeStack[--stackIdx];
		const BVH_DOBB::BVHNode& currNode = bvh.dobbNode[currNodeIdx];

		Ray localRay;
		localRay.hit.t = ray.hit.t;
		const uint32_t rotationIdx = currNode.Rotation();

		if (rotationIdx != 255) {
			const DOBB_LUT::Mat3f& rotation = bvh.getLutEntry(rotationIdx);

			const bvhvec3 rotationX = bvhvec3(rotation.m[0][0], rotation.m[1][0], rotation.m[2][0]);
			const bvhvec3 rotationY = bvhvec3(rotation.m[0][1], rotation.m[1][1], rotation.m[2][1]);
			const bvhvec3 rotationZ = bvhvec3(rotation.m[0][2], rotation.m[1][2], rotation.m[2][2]);

			localRay.O = bvhvec3(tinybvh_dot(rotationX, ray.O), tinybvh_dot(rotationY, ray.O), tinybvh_dot(rotationZ, ray.O));
			localRay.D = bvhvec3(tinybvh_dot(rotationX, ray.D), tinybvh_dot(rotationY, ray.D), tinybvh_dot(rotationZ, ray.D));
		}
		else {
			localRay.O = ray.O; localRay.D = ray.D;
		}
		localRay.rD = tinybvh_rcp(localRay.D);

		float stepSize[3] = {ldexpf(1.0f, currNode.exponent[0]), ldexpf(1.0f, currNode.exponent[1]), ldexpf(1.0f, currNode.exponent[2])};

		uint32_t hitMask = 0;
		uint32_t remainingSlots = currNode.OccupiedMask();

		while (remainingSlots) {
			if (stats) stats->boxTests++;

			// Determine which is the next slot to process
			uint32_t slot = (uint32_t)std::countr_zero(remainingSlots);
			remainingSlots &= remainingSlots - 1;

			// Determine lower and higher extents of the BB
			bvhvec3 slotMinPoint = currNode.nodeMinPoint + bvhvec3(currNode.qLo[0][slot] * stepSize[0], currNode.qLo[1][slot] * stepSize[1], currNode.qLo[2][slot] * stepSize[2]);
			bvhvec3 slotMaxPoint = currNode.nodeMinPoint + bvhvec3(currNode.qHi[0][slot] * stepSize[0], currNode.qHi[1][slot] * stepSize[1], currNode.qHi[2][slot] * stepSize[2]);

			float entryDistance = tinybvh_intersect_aabb(localRay, slotMinPoint, slotMaxPoint);
			if (entryDistance < BVH_FAR) 
				hitMask |= 1 << slot;
		}

		uint8_t leafHitMask = hitMask & ~currNode.iMask;
		uint8_t internalHitMask = hitMask & currNode.iMask;

		// Process leaf hits
		while (leafHitMask) {
			const uint32_t slot = (uint32_t)std::countr_zero(leafHitMask);
			leafHitMask &= leafHitMask - 1;

			// Intersect triangles
			const uint32_t firstTri = currNode.triBaseIndex + (currNode.meta[slot] & 31);
			const uint32_t triCount = __popc(currNode.meta[slot] >> 5);

			uint32_t i0, i1, i2;
			for (uint32_t triIdx = 0; triIdx < triCount; triIdx++) {
				if (stats) stats->triTests++;

				const uint32_t primIdx = bvh.dobbPrimIdx[firstTri + triIdx];
				GET_PRIM_INDICES_I0_I1_I2(bvh.bvh8.bvh, primIdx);
				if (TriOccludes(ray, bvh.bvh8.bvh.verts, primIdx, i0, i1, i2))
					return true;
			}
		}

		// Process internal hits
		while (internalHitMask) {
			const uint32_t slot = (uint32_t)std::countr_zero(internalHitMask);
			internalHitMask &= internalHitMask - 1;

			// Push internal node
			nodeStack[stackIdx++] = currNode.ChildBaseIndex() +
				__popc(currNode.iMask & ((1u << slot) - 1));
		}
	}

	return false;
}

// MBVH Intersect matching with DOBB processing order
void IntersectionInstrumented(const MBVH<8>& bvh8, Ray& ray, TraversalStats* stats) {
	if (stats) stats->rayCount++;

	uint32_t nodeStack[512];
	uint32_t stackIdx = 0;
	nodeStack[stackIdx++] = 0;

	while (stackIdx > 0) {
		const MBVH<8>::MBVHNode& currNode = bvh8.mbvhNode[nodeStack[--stackIdx]];

		if (currNode.isLeaf()) {
			uint32_t i0, i1, i2;
			for (uint32_t triOffset = 0; triOffset < currNode.triCount; triOffset++) {
				if (stats) stats->triTests++;

				const uint32_t primIdx = bvh8.bvh.primIdx[currNode.firstTri + triOffset];
				GET_PRIM_INDICES_I0_I1_I2(bvh8.bvh, primIdx);
				IntersectTri(ray, primIdx, bvh8.bvh.verts, i0, i1, i2);
			}
		}
		else {
			// Leaves are not nodes in DOBB
			if (stats) stats->nodeVisits++;

			// Make sure nodes are traversed w the same logic as DOBB
			bool childHit[8] = {};
			for (uint32_t slot = 0; slot < 8; slot++) {
				if (currNode.child[slot] != 0) {
					if (stats) stats->boxTests++;

					const MBVH<8>::MBVHNode& child = bvh8.mbvhNode[currNode.child[slot]];
					childHit[slot] = tinybvh_intersect_aabb(ray, child.aabbMin, child.aabbMax) < BVH_FAR;;
				}
			}

			for (uint32_t slot = 0; slot < 8; slot++)
				if (childHit[slot] && !bvh8.mbvhNode[currNode.child[slot]].isLeaf())
					nodeStack[stackIdx++] = currNode.child[slot];

			for (uint32_t slot = 0; slot < 8; slot++)
				if (childHit[slot] && bvh8.mbvhNode[currNode.child[slot]].isLeaf())
					nodeStack[stackIdx++] = currNode.child[slot];
		}
	}
}

// MBVH isOccluded matching with DOBB processing order
bool isOccludedInstrumented(const MBVH <8>& bvh8, const Ray& ray, TraversalStats* stats) {
	if (stats) stats->rayCount++;

	Ray boxRay = ray;

	uint32_t nodeStack[512];
	uint32_t stackIdx = 0;
	nodeStack[stackIdx++] = 0;

	while (stackIdx > 0) {
		const MBVH<8>::MBVHNode& currNode = bvh8.mbvhNode[nodeStack[--stackIdx]];

		if (currNode.isLeaf()) {
			uint32_t i0, i1, i2;
			for (uint32_t triOffset = 0; triOffset < currNode.triCount; triOffset++) {
				if (stats) stats->triTests++;

				const uint32_t primIdx = bvh8.bvh.primIdx[currNode.firstTri + triOffset];
				GET_PRIM_INDICES_I0_I1_I2(bvh8.bvh, primIdx);
				if (TriOccludes(ray, bvh8.bvh.verts, primIdx, i0, i1, i2))
					return true;
			}
		}
		else {
			// Leaves are not nodes in DOBB
			if (stats) stats->nodeVisits++;

			// Make sure nodes are traversed w the same logic as DOBB
			bool childHit[8] = {};
			for (uint32_t slot = 0; slot < 8; slot++) {
				if (currNode.child[slot] != 0) {
					if (stats) stats->boxTests++;

					const MBVH<8>::MBVHNode& child = bvh8.mbvhNode[currNode.child[slot]];
					childHit[slot] = tinybvh_intersect_aabb(boxRay, child.aabbMin, child.aabbMax) < BVH_FAR;;
				}
			}

			for (uint32_t slot = 0; slot < 8; slot++)
				if (childHit[slot] && !bvh8.mbvhNode[currNode.child[slot]].isLeaf())
					nodeStack[stackIdx++] = currNode.child[slot];

			for (uint32_t slot = 0; slot < 8; slot++)
				if (childHit[slot] && bvh8.mbvhNode[currNode.child[slot]].isLeaf())
					nodeStack[stackIdx++] = currNode.child[slot];
		}
	}

	return false;
}

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
template <class BVHType>
bvhvec3 Trace(const BVHType& bvh, Ray ray, unsigned& seed, uint64_t& traceCost, 
	TraversalStats& threadHitStats, TraversalStats& threadShadowStats, unsigned depth = 0)
{
	// find primary intersection
	IntersectionInstrumented(bvh, ray, &threadHitStats);
	// shade
	if (ray.hit.t == 1e30f)
		return bvhvec3(0.6f, 0.7f, 1); // hit nothing
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
		if (!isOccludedInstrumented(bvh, Ray(I + L * 0.001f, L, dist), &threadShadowStats))
			direct = BRDF * NdotL * NLdotL * bvhvec3(9, 9, 8) * 500 * (1.0f / (dist * dist));
	// random bounce
	if (depth < 2)
	{
		bvhvec3 R = CosWeightedDiffReflection(N, seed);
		float pdf = 1.0f / tinybvh_dot(N, R);
		bvhvec3 irradiance = Trace(bvh, Ray(I + R * 0.001f, R), seed, traceCost, threadHitStats, threadShadowStats, depth + 1);
		indirect = BRDF * irradiance * (1.0f / pdf);
	}
	// finalize
	return direct + indirect;
}

template <class BVHType>
void TraceWorkerThread(const BVHType& bvh, uint32_t* buf, float scale, int threadIdx,
	TraversalStats& threadHitStats, TraversalStats& threadShadowStats)
{
	uint64_t localCost = 0;
	const int xtiles = SCRWIDTH / TILESIZE, ytiles = SCRHEIGHT / TILESIZE;
	const int tiles = xtiles * ytiles;
	int tile = threadIdx;
	while (tile < tiles)
	{
		const int tx = tile % xtiles, ty = tile / xtiles;
		unsigned seed = (tile + 17) * 171717 + spp * 1023;
		for (int y = 0; y < TILESIZE; y++) for (int x = 0; x < TILESIZE; x++)
		{
			const int pixel_x = tx * TILESIZE + x, pixel_y = ty * TILESIZE + y;
			const int pixelIdx = pixel_x + pixel_y * SCRWIDTH;
			// setup primary ray
			const float u = (float)pixel_x / SCRWIDTH, v = (float)pixel_y / SCRHEIGHT;
			const bvhvec3 D = tinybvh_normalize(p1 + u * (p2 - p1) + v * (p3 - p1) - eye);
			// trace
			const TraversalStats hitBefore = threadHitStats, shadowBefore = threadShadowStats;
			Trace(bvh, Ray(eye, D), seed, localCost, threadHitStats, threadShadowStats);
			// red: total node visits; green: total triangle tests
			accumulator[pixelIdx] += bvhvec3(
				nodeVisitScale * ((threadHitStats.nodeVisits - hitBefore.nodeVisits) + (threadShadowStats.nodeVisits - shadowBefore.nodeVisits)),
				triTestScale * ((threadHitStats.triTests - hitBefore.triTests) + (threadShadowStats.triTests - shadowBefore.triTests)),
				0.0f
			);
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
void AddMesh(const char* file, float scale = 1, bvhvec3 pos = {}, int c = 0,
	bvhvec3 rotationAxis = bvhvec3(0, 1, 0), float rotationDegrees = 0.0f)
{
	std::fstream s{ file, s.binary | s.in };
	int N = 0;
	s.read((char*)&N, 4);
	bvhvec4* data = (bvhvec4*)malloc64((N + triCount) * 48);
	if (tris) memcpy(data, tris, triCount * 48), free64(tris);
	tris = data, s.read((char*)tris + triCount * 48, N * 48), triCount += N;
	
	rotationAxis = tinybvh_normalize(rotationAxis);
	const float angle = rotationDegrees * (3.14159265f / 180.0f);
	const float cosAngle = cosf(angle);
	const float sinAngle = sinf(angle);

	for (int vertexIdx = (triCount - N) * 3; vertexIdx < triCount * 3; vertexIdx++) {
		bvhvec4& vertex = tris[vertexIdx];
		// Scale vertex position
		bvhvec3 vertPos = bvhvec3(vertex) * scale;
		// Rotate vertex position 
		vertPos = vertPos * cosAngle
			+ tinybvh_cross(rotationAxis, vertPos) * sinAngle
			+ rotationAxis * (tinybvh_dot(rotationAxis, vertPos) * (1.0f - cosAngle));
		// Add pos, to get world position
		vertPos = vertPos + pos;

		vertex.x = vertPos.x;
		vertex.y = vertPos.y;
		vertex.z = vertPos.z;
		if (c) *(int*)&vertex.w = c;
	}
}

void ClearMesh()
{
	if (tris)
	{
		free64(tris);
		tris = 0;
	}

	triCount = 0;
}

// ---------------------- BUILD BVH ----------------------

// Build BVHs and then set their c_trav and c_int values to 1.0
// to measure BVH quality
void InitBvhs() 
{
	dobbBvh.Build(tris, triCount);
}

// ---------------------- MISC ----------------------

void Init() {
#ifdef _WIN32
	AllocConsole();
	FILE* dummy;
	freopen_s(&dummy, "CONOUT$", "w", stdout);
	freopen_s(&dummy, "CONOUT$", "w", stderr);
	std::cout.clear();
	std::cerr.clear();
#endif
	// load camera position / direction from file
	std::fstream t = std::fstream{ "camera.bin", t.binary | t.in };
	if (t.is_open()) {
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

template <class BVHType>
void RenderTick(const BVHType& bvh, float delta_time_s, fenster& f, uint32_t* buf)
{
	// render tiles
	const float scale = 1.0f / ++spp;
	tileIdx = threadCount;
	std::vector<TraversalStats> threadHitStats(threadCount), threadShadowStats(threadCount);
	std::vector<std::thread> threads;
	for (uint32_t i = 0; i < threadCount; i++)
		threads.emplace_back(TraceWorkerThread<BVHType>, std::cref(bvh), buf, scale, i,
			std::ref(threadHitStats[i]), std::ref(threadShadowStats[i]));
	for (auto& thread : threads) thread.join();
	for (uint32_t i = 0; i < threadCount; i++)
		hitStats.add(threadHitStats[i]), shadowStats.add(threadShadowStats[i]);
	// print frame time / rate in window title
	char title[50];
	snprintf(
		title, sizeof(title), 
		"tiny_bvh %s, spp: %d", 
		meshNames[meshIdx], spp
	);
	fenster_update_title(&f, title);
}

void resetRender() {
	memset(accumulator, 0, sizeof(accumulator));
	spp = 0;

	hitStats.reset(); shadowStats.reset();
}

void Tick(float delta_time_s, fenster& f, uint32_t* buf)
{
	switch (phase) {
		case BUILD:
			// Reset and load next mesh to be tested
			std::cout << "Loading mesh: " << meshNames[meshIdx] << std::endl;
			ClearMesh();
			AddMesh(
				meshLocations[meshIdx], 
				1, bvhvec3(0), 0xffffff,
				bvhvec3(0.0f, 1.0f, 0.0f), 45.0f 
			);
			std::cout << "Building BVHs" << std::endl;
			// Build BVHs 
			InitBvhs();
			
			// Swap to the next phase
			phase = RENDERMBVH;
			resetRender();
			break;
		case RENDERMBVH:
			RenderTick(dobbBvh.bvh8, delta_time_s, f, buf);
			if (spp >= targetSpp)
			{
				std::cout << "MBVH8 BVH, mesh: " << meshNames[meshIdx] << std::endl;
				hitStats.print("Closest-hit rays");
				shadowStats.print("Shadow rays");
				hitStats.add(shadowStats);
				hitStats.print("All rays");
				phase = RENDERDOBB;
				resetRender();
			}
			break;

		case RENDERDOBB:
			RenderTick(dobbBvh, delta_time_s, f, buf);
			if (spp >= targetSpp)
			{
				std::cout << "DOBB BVH, mesh: " << meshNames[meshIdx] << std::endl;
				hitStats.print("Closest-hit rays");
				shadowStats.print("Shadow rays");
				hitStats.add(shadowStats);
				hitStats.print("All rays");
				phase = BUILD;
				meshIdx++;
				if (meshIdx >= meshCount) 
				{
					phase = END;
					std::cout << "END" << std::endl;
				}
			}
			break;
		case END: break;
		default: break;
	}
}

// Application Shutdown
void Shutdown()
{
	// save camera position / direction to file
	std::fstream s = std::fstream{ "camera.bin", s.binary | s.out };
	s.write((char*)&eye, sizeof(eye));
	s.write((char*)&view, sizeof(view));
	s.close();
#ifdef _WIN32
	FreeConsole();
#endif
}