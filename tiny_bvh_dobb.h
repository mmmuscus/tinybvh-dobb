#ifndef TINY_BVH_DOBB_H_
#define TINY_BVH_DOBB_H_

#ifndef TINY_BVH_H_
#include "tiny_bvh.h"
#endif

#define AXESNO 13
#define SMALLM 4

namespace tinybvh {
	static const bvhvec3 proxyKDopAxes[AXESNO] = {
		// 3 Euclidian axes
		bvhvec3(1.0f, 0.0f, 0.0f),
		bvhvec3(0.0f, 1.0f, 0.0f),
		bvhvec3(0.0f, 0.0f, 1.0f),
		// 6 diagonal axes
		tinybvh_normalize(bvhvec3(1.0f, 1.0f, 0.0f)),
		tinybvh_normalize(bvhvec3(1.0f, -1.0f, 0.0f)),
		tinybvh_normalize(bvhvec3(1.0f, 0.0f, 1.0f)),
		tinybvh_normalize(bvhvec3(1.0f, 0.0f, -1.0f)),
		tinybvh_normalize(bvhvec3(0.0f, 1.0f, 1.0f)),
		tinybvh_normalize(bvhvec3(0.0f, 1.0f, -1.0f)),
		// 4 "axes in 3D space"
		// pointing towards the vertices of the basis cube
		tinybvh_normalize(bvhvec3(1.0f, 1.0f, 1.0f)),
		tinybvh_normalize(bvhvec3(1.0f, 1.0f, -1.0f)),
		tinybvh_normalize(bvhvec3(1.0f, -1.0f, 1.0f)),
		tinybvh_normalize(bvhvec3(1.0f, -1.0f, -1.0f)),
	};

	// TODO: indirection of LUT
	class DOBB_LUT {
	public: 
		struct Mat3 { float m[3][3]; };

		Mat3 rotations[AXESNO * SMALLM * 2];

		DOBB_LUT() { buildLUT(); }

		const float similarity(const bvhvec3 (&capitalB)[3], const uint32_t rotIdx) {
			const Mat3 rotMat = rotations[rotIdx];
			// Using dot product, set up matrix that we can use to compare
			// all possible permutations of bases
			Mat3 dotMat;

			for (uint32_t i = 0; i < 3; i++) {
				const bvhvec3 r(rotMat.m[0][i], rotMat.m[1][i], rotMat.m[2][i]);
				
				for (uint32_t j = 0; j < 3; j++)
					dotMat.m[i][j] = fabsf(tinybvh_dot(r, capitalB[j]));
			}

			static const int Perms[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
			float maxSim = 0.0f;
			for (uint32_t i = 0; i < 6; i++)
				maxSim = tinybvh_max(
					maxSim,
					dotMat.m[0][Perms[i][0]] +
					dotMat.m[1][Perms[i][1]] +
					dotMat.m[2][Perms[i][2]]
				);

			return maxSim;
		}

	private:
		// https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
		// axis must be normalized
		static Mat3 rodriguesRotMat(bvhvec3 axis, float angle) {
			const float cs = cosf(angle);
			const float oneMinusCos = 1 - cs;
			const float sn = sinf(angle);
			Mat3 ret;

			ret.m[0][0] = cs + oneMinusCos * axis.x * axis.x;
			ret.m[0][1] = oneMinusCos * axis.x * axis.y - sn * axis.z;
			ret.m[0][2] = oneMinusCos * axis.x * axis.z + sn * axis.y;

			ret.m[1][0] = oneMinusCos * axis.x * axis.y + sn * axis.z;
			ret.m[1][1] = cs + oneMinusCos * axis.y * axis.y;
			ret.m[1][2] = oneMinusCos * axis.y * axis.z - sn * axis.x;

			ret.m[2][0] = oneMinusCos * axis.x * axis.z - sn * axis.y;
			ret.m[2][1] = oneMinusCos * axis.y * axis.z + sn * axis.x;
			ret.m[2][2] = cs + oneMinusCos * axis.z * axis.z;

			return ret;
		}

		void buildLUT() {
			const float delta = 3.14159265358979f / (2.0f * SMALLM);
			uint32_t rotationsIdx = 0;

			for (uint32_t axisIdx = 0; axisIdx < AXESNO; axisIdx++) 
				for (int i = 1; i <= SMALLM; i++) {
					rotations[rotationsIdx++] = rodriguesRotMat(
						proxyKDopAxes[axisIdx], i * delta);

					rotations[rotationsIdx++] = rodriguesRotMat(
						proxyKDopAxes[axisIdx], i * -delta);
				}
		}
	};

	class BVH_DOBB : public BVHBase {
	public:
		// TODO: do compressed BVHNode
		struct BVHNode
		{
			uint8_t rotation;
		};

		// TODO: move outside of class
		// Only implementing |K| = 13 case
		struct kDop
		{
			// Index corresponds to axes
			bvhvec2 extents[AXESNO];
		};

		BVH_DOBB(BVHContext ctx = {}) { context = ctx; } // TODO: add layout
		~BVH_DOBB();
		void ConvertFrom(const MBVH<8>& original);

		// BVH data
		BVHNode* dobbNode = 0;
		MBVH<8> bvh8;
		bool ownBVH8 = true;

		kDop* proxyKDop = 0;
		DOBB_LUT lut;
};

} // namespace tinybvh
#endif // TINY_BVH_DOBB_H_


#if defined( TINYBVH_DOBB_IMPLEMENTATION ) && !defined( TINY_BVH_DOBB_IMPL_DONE )
#define TINY_BVH_DOBB_IMPL_DONE

namespace tinybvh {
	BVH_DOBB::~BVH_DOBB()
	{
		if (!ownBVH8) bvh8.ReleaseOwnership();
		AlignedFree(dobbNode);

		AlignedFree(proxyKDop);
	}

	void BVH_DOBB::ConvertFrom(const MBVH<8>& original)
	{
		// copy source
		if (&original != &bvh8)
		{
			if (ownBVH8) bvh8 = MBVH<8>(context); 
			else bvh8.DropReference(context);
			ownBVH8 = false;
		}
		bvh8.ReferenceFrom(original);
		// get base properties
		CopyBasePropertiesFrom(bvh8);
		// allocate if need be
		const uint32_t nodesNeeded = bvh8.usedNodes;
		if (allocatedNodes < nodesNeeded)
		{
			// Allocate nodes
			AlignedFree(dobbNode);
			dobbNode = (BVHNode*)AlignedAlloc(nodesNeeded * sizeof(BVHNode));
			allocatedNodes = nodesNeeded;

			// Allocate proxy k-DOPs
			AlignedFree(proxyKDop);
			proxyKDop = (kDop*)AlignedAlloc(nodesNeeded * sizeof(kDop));
		}
		usedNodes = nodesNeeded;

		// Generate LUT
		lut = DOBB_LUT();

		// Iterate through nodes backwards (children are processed
		// implicitly before parents)
		for (int nodeIdx = usedNodes - 1; nodeIdx >= 0; nodeIdx--)
		{
			// Reset proxy k-DOP for min/max selection
			for (int axisIdx = 0; axisIdx < AXESNO; axisIdx++)
			{
				proxyKDop[nodeIdx].extents[axisIdx] = bvhvec2(BVH_FAR, -BVH_FAR);
			}

			const auto& currNode = bvh8.mbvhNode[nodeIdx];

			// Decide if leaf or not
			if (currNode.isLeaf())
			{
				// Process leaf nodes
				bvhvec3 leafVerts[6];
				uint32_t noLeafVerts = 0;
				// used for GET_PRIM_INDICES_I0_I1_I2 macro
				uint32_t i0, i1, i2;

				// For determining B;
				bvhvec3 a0, a1, a2;

				if (currNode.triCount == 1) {
					noLeafVerts = 3;

					// Gather vertices in leaf
					const uint32_t prim = bvh8.bvh.primIdx[currNode.firstTri];
					GET_PRIM_INDICES_I0_I1_I2(bvh8.bvh, prim);

					leafVerts[0] = bvh8.bvh.verts[i0];
					leafVerts[1] = bvh8.bvh.verts[i1];
					leafVerts[2] = bvh8.bvh.verts[i2];

					// TODO: Create helper function to copy in the non-quad triCount == 2 case
					// Determine a0 and a1 for B
					// Determine longest edge of triangle => find 
					// "middle point" of shortest edges
					const bvhvec3 edges[3] = {
						leafVerts[1] - leafVerts[0], // (v0, v1) = e0
						leafVerts[2] - leafVerts[1], // (v1, v2) = e1
						leafVerts[0] - leafVerts[2]  // (v2, v0) = e2
					};
					// Get lengths
					const float lengthsSqd[3] = {
						tinybvh_dot(edges[0], edges[0]),
						tinybvh_dot(edges[1], edges[1]),
						tinybvh_dot(edges[2], edges[2])
					};

					uint32_t longestIdx = lengthsSqd[0] > lengthsSqd[1] ? 0 : 1;
					if (lengthsSqd[2] > lengthsSqd[longestIdx]) longestIdx = 2;

					uint32_t midPointIdx = (longestIdx + 2) % 3;

					if (lengthsSqd[midPointIdx] < lengthsSqd[(midPointIdx + 2) % 3]) {
						a0 = edges[midPointIdx]; 
						a1 = -edges[(midPointIdx + 2) % 3];
					}
					else { 
						a0 = -edges[(midPointIdx + 2) % 3]; 
						a1 = edges[midPointIdx];
					}
				}
				// Implicitly, leaf has 2 triangles
				else
				{
					const uint32_t prim0 = bvh8.bvh.primIdx[currNode.firstTri];
					GET_PRIM_INDICES_I0_I1_I2(bvh8.bvh, prim0);

					leafVerts[0] = bvh8.bvh.verts[i0];
					leafVerts[1] = bvh8.bvh.verts[i1];
					leafVerts[2] = bvh8.bvh.verts[i2];

					const uint32_t prim1 = bvh8.bvh.primIdx[currNode.firstTri + 1];
					GET_PRIM_INDICES_I0_I1_I2(bvh8.bvh, prim1);

					leafVerts[3] = bvh8.bvh.verts[i0];
					leafVerts[4] = bvh8.bvh.verts[i1];
					leafVerts[5] = bvh8.bvh.verts[i2];

					// Determine shared vertices
					uint32_t sharedVert1[2] = { 6, 6 };
					uint32_t sharedVert2[2] = { 6, 6 };

					for (uint32_t tri1VertIdx = 0; tri1VertIdx < 3; tri1VertIdx++) {
						for (uint32_t tri2VertIdx = 3; tri2VertIdx < 6; tri2VertIdx++) {
							if (leafVerts[tri1VertIdx].x == leafVerts[tri2VertIdx].x &&
								leafVerts[tri1VertIdx].y == leafVerts[tri2VertIdx].y &&
								leafVerts[tri1VertIdx].z == leafVerts[tri2VertIdx].z)
							{
								if (sharedVert1[0] == 6) {
									sharedVert1[0] = tri1VertIdx;
									sharedVert1[1] = tri2VertIdx;
								}
								else {
									sharedVert2[0] = tri1VertIdx;
									sharedVert2[1] = tri2VertIdx;
								}
							}

						}
					}

					// Check if triangles share an edge 
					if (sharedVert2[0] == 6) {
						// TODO: handle case where triangles share exactly 1 vertex
						// This case will be selected if the triangles match
						// in at most one vertex, thus here there is a wasted
						// projection
						noLeafVerts = 6;

						// Determine a0 and a1 for B
						// Determine longest edge of triangle => find 
						// "middle point" of shortest edges
						const bvhvec3 edges[3] = {
							leafVerts[1] - leafVerts[0], // (v0, v1) = e0
							leafVerts[2] - leafVerts[1], // (v1, v2) = e1
							leafVerts[0] - leafVerts[2]  // (v2, v0) = e2
						};
						// Get lengths
						const float lengthsSqd[3] = {
							tinybvh_dot(edges[0], edges[0]),
							tinybvh_dot(edges[1], edges[1]),
							tinybvh_dot(edges[2], edges[2])
						};

						uint32_t longestIdx = lengthsSqd[0] > lengthsSqd[1] ? 0 : 1;
						if (lengthsSqd[2] > lengthsSqd[longestIdx]) longestIdx = 2;

						uint32_t midPointIdx = (longestIdx + 2) % 3;

						if (lengthsSqd[midPointIdx] < lengthsSqd[(midPointIdx + 2) % 3]) {
							a0 = edges[midPointIdx];
							a1 = -edges[(midPointIdx + 2) % 3];
						}
						else {
							a0 = -edges[(midPointIdx + 2) % 3];
							a1 = edges[midPointIdx];
						}
					}
					else {
						noLeafVerts = 4;

						// Init tmpleagVerts for swaps
						bvhvec3 tmpleafVerts[6];
						for (uint32_t i = 0; i < 6; i++)
							tmpleafVerts[i] = leafVerts[i];

						// Insert first three vertex, making sure the shared
						// edge is in slots 1 and 2
						for (uint32_t i = 0; i < 3; i++) {
							if (i != sharedVert1[0] && i != sharedVert2[0])
								leafVerts[0] = tmpleafVerts[i];
							if (i == sharedVert1[0])
								leafVerts[1] = tmpleafVerts[i];
							if (i == sharedVert2[0])
								leafVerts[2] = tmpleafVerts[i];
						}

						// Insert final vertex
						for (uint32_t i = 3; i < 6; i++)
							if (i != sharedVert1[1] && i != sharedVert2[1])
								leafVerts[3] = tmpleafVerts[i];

						// Determine a0 and a1 for B
						// Length does not matter in this case, only that
						// the "middle" vertex is in the shared edge
						a0 = leafVerts[0] - leafVerts[1];
						a1 = leafVerts[3] - leafVerts[1];
					}
				}

				// Calculate B
				a0 = tinybvh_normalize(a0);
				a2 = tinybvh_normalize(tinybvh_cross(a0, a1));
				a1 = tinybvh_normalize(tinybvh_cross(a2, a0)); // normalized for safety
				
				// TODO: Implement axis azimuth version instead of brute force
				// find best rotation candidate for leaf node with B
				float bestSim = 0.0f;
				for (uint8_t rotationIdx = 0; rotationIdx < AXESNO * SMALLM * 2; rotationIdx++) {
					float sim = lut.similarity({a0, a1, a2}, rotationIdx);
					if (sim > bestSim) {
						bestSim = sim;
						dobbNode[nodeIdx].rotation = rotationIdx;
					}
				}

				// Create proxy KDop to propagate up
				for (uint32_t vertIdx = 0; vertIdx < noLeafVerts; vertIdx++) {
					for (uint32_t axisIdx = 0; axisIdx < AXESNO; axisIdx++) {
						proxyKDop[nodeIdx].extents[axisIdx].x = tinybvh_min(
							tinybvh_dot(proxyKDopAxes[axisIdx], leafVerts[vertIdx]),
							proxyKDop[nodeIdx].extents[axisIdx].x);

						proxyKDop[nodeIdx].extents[axisIdx].y = tinybvh_max(
							tinybvh_dot(proxyKDopAxes[axisIdx], leafVerts[vertIdx]),
							proxyKDop[nodeIdx].extents[axisIdx].y);
					}
				}
			}
			else
			{
				float AABBSAMax = 0.0f;
				uint8_t candidateRotationIdx = 127;	// invalid value
				float AABBSASum = 0.0f;

				// Process internal nodes
				for (uint32_t childIdx = 0; childIdx < currNode.childCount; childIdx++)
				{
					uint32_t currChildIdx = currNode.child[childIdx];

					// Gather extents from euclidean axes from axes collection
					float a = fabsf(proxyKDop[currChildIdx].extents[0].x - proxyKDop[currChildIdx].extents[0].y);
					float b = fabsf(proxyKDop[currChildIdx].extents[1].x - proxyKDop[currChildIdx].extents[1].y);
					float c = fabsf(proxyKDop[currChildIdx].extents[2].x - proxyKDop[currChildIdx].extents[2].y);

					// Select rotation by maxselecting on AABB surface area
					float childAABBSA = 2.0f * (a * b + b * c + c * a);
					if (childAABBSA > AABBSAMax) {
						AABBSAMax = childAABBSA;
						candidateRotationIdx = dobbNode[currChildIdx].rotation;
					}

					// Determine k-DOP extents
					for (uint32_t axisIdx = 0; axisIdx < AXESNO; axisIdx++)
					{
						proxyKDop[nodeIdx].extents[axisIdx].x = tinybvh_min(
							proxyKDop[currChildIdx].extents[axisIdx].x,
							proxyKDop[nodeIdx].extents[axisIdx].x);

						proxyKDop[nodeIdx].extents[axisIdx].y = tinybvh_max(
							proxyKDop[currChildIdx].extents[axisIdx].y,
							proxyKDop[nodeIdx].extents[axisIdx].y);
					}
				}
			}
		}
	}

} // namespace tinybvh
#endif // TINYBVH_DOBB_IMPLEMENTATION