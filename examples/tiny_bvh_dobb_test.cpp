#define TINYBVH_IMPLEMENTATION
#include "tiny_bvh.h"
#include <cstdlib>
#include <cstdio>

// Adding functionality from minimal example to test
static constexpr uint32_t TRIANGLE_COUNT = 8192;
static tinybvh::bvhvec4 vertices[TRIANGLE_COUNT * 3]; // must be 16 byte!

float uniform_rand() { return (float)rand() / (float)RAND_MAX; }

int main()
{
	printf("Hello tinybvh!\n");

	// Adding functionality from minimal example to test
	
	// create a scene consisting of some random small triangles
	for (uint32_t i = 0; i < TRIANGLE_COUNT; i++)
	{
		// create a random triangle
		tinybvh::bvhvec4& v0 = vertices[i * 3 + 0];
		tinybvh::bvhvec4& v1 = vertices[i * 3 + 1];
		tinybvh::bvhvec4& v2 = vertices[i * 3 + 2];
		// triangle position, x/y/z = 0..1
		float x = uniform_rand();
		float y = uniform_rand();
		float z = uniform_rand();
		// set first vertex
		v0.x = x + 0.1f * uniform_rand();
		v0.y = y + 0.1f * uniform_rand();
		v0.z = z + 0.1f * uniform_rand();
		// set second vertex
		v1.x = x + 0.1f * uniform_rand();
		v1.y = y + 0.1f * uniform_rand();
		v1.z = z + 0.1f * uniform_rand();
		// set third vertex
		v2.x = x + 0.1f * uniform_rand();
		v2.y = y + 0.1f * uniform_rand();
		v2.z = z + 0.1f * uniform_rand();
	}

	// construct a ray
	tinybvh::bvhvec3 O(0.5f, 0.5f, -1);
	tinybvh::bvhvec3 D(0.1f, 0, 2);
	tinybvh::Ray ray(O, D);

	// build a BVH over the scene
	tinybvh::BVH bvh;
	bvh.Build(vertices, TRIANGLE_COUNT);

	// from here: play with the BVH!
	int steps = bvh.Intersect(ray);
	if (ray.hit.t < BVH_FAR)
		printf("hit prim %u at t=%f (%i steps)\n", ray.hit.prim, ray.hit.t, steps);
	else
		printf("no hit (%i steps)\n", steps);
}