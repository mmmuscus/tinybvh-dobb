#ifndef TINY_BVH_DOBB_H_
#define TINY_BVH_DOBB_H_

#ifndef TINY_BVH_H_
#include "tiny_bvh.h"
#endif

#define AXESNO 13
#define SMALLM 4

#ifndef FATAL_ERROR_IF
#define FATAL_ERROR_IF(condition, message) if (condition) { fprintf( stderr, "Fatal error in tiny_bvh_dobb.h, line %i:\n%s\n", __LINE__, message ); exit( 1 ); }
#endif

// TODO: General refactor:
//	- use const consistently
//  - remove uneccessary declaration
//	- you can index bvhvec3-s, use it where its useful!

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

		struct reprojectionEntry {
			uint8_t axisIdx[3];
			int sign[3];
		};

		Mat3 rotations[AXESNO * SMALLM * 2];
		reprojectionEntry reprojections[AXESNO * SMALLM * 2][6];

		DOBB_LUT() { 
			buildLUT(); 
			buildReprojectionLUT();
		}

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

		void buildReprojectionLUT() {
			for (uint32_t rotationIdx = 0; rotationIdx < AXESNO * SMALLM * 2; rotationIdx++) {
				const Mat3& rotation = rotations[rotationIdx];
				for (uint32_t candidateRotAxisIdx = 0; candidateRotAxisIdx < 3; candidateRotAxisIdx++) {
					const bvhvec3 candidateAxis(
						rotation.m[0][candidateRotAxisIdx],
						rotation.m[1][candidateRotAxisIdx],
						rotation.m[2][candidateRotAxisIdx]
					);

					findQuadSphereTriangle(
						candidateAxis, 
						reprojections[rotationIdx][candidateRotAxisIdx * 2]
					);
					findQuadSphereTriangle(
						-candidateAxis,
						reprojections[rotationIdx][candidateRotAxisIdx * 2 + 1]
					);
				}
			}
		}

		// Each triangle of the quad sphere for out purposes is 
		// corresponds to a face center,an edge middle point and a 
		// corner (of the unit cube)
		void findQuadSphereTriangle(const bvhvec3& direction, reprojectionEntry& ret) {
			const float components[3] = { direction.x, direction.y, direction.z };

			// Find face center
			float faceCenter[3] = {0.0f, 0.0f, 0.0f};

			uint8_t majorAxis = 0;
			if (fabsf(components[1]) > fabsf(components[majorAxis])) 
				majorAxis = 1;
			if (fabsf(components[2]) > fabsf(components[majorAxis]))
				majorAxis = 2;

			faceCenter[majorAxis] = components[majorAxis] >= 0.0f ? 1.0f : -1.0f;
			quadSpherePointToProxyAxis(faceCenter, ret.axisIdx[0], ret.sign[0]);

			// Find cube corner
			float cubeCorner[3] = {faceCenter[0], faceCenter[1], faceCenter[2]};

			int minorAxisFirst = (majorAxis + 1) % 3;
			int minorAxisSecond = (majorAxis + 2) % 3;
			const float minorAxisFirstSign = components[minorAxisFirst] >= 0.0f ? 1.0f : -1.0f;
			const float minorAxisSecondSign = components[minorAxisSecond] >= 0.0f ? 1.0f : -1.0f;
			
			cubeCorner[minorAxisFirst] = minorAxisFirstSign;
			cubeCorner[minorAxisSecond] = minorAxisSecondSign;
			quadSpherePointToProxyAxis(cubeCorner, ret.axisIdx[1], ret.sign[1]);

			// Find edge midpoint
			float edgeMidPoint[3] = {faceCenter[0], faceCenter[1], faceCenter[2]};
			
			const bool minorAxisFirstIsLarger = fabsf(components[minorAxisFirst]) >= fabsf(components[minorAxisSecond]);
			if (minorAxisFirstIsLarger)
				edgeMidPoint[minorAxisFirst] = minorAxisFirstSign;
			else
				edgeMidPoint[minorAxisSecond] = minorAxisSecondSign;
			quadSpherePointToProxyAxis(edgeMidPoint, ret.axisIdx[2], ret.sign[2]);
		}

		void quadSpherePointToProxyAxis(const float quadSpherePoint[3], uint8_t& axisIdx, int& sign) {
			bvhvec3 quadSpherePointVec = {quadSpherePoint[0], quadSpherePoint[1], quadSpherePoint[2] };
			quadSpherePointVec = tinybvh_normalize(quadSpherePointVec);
			
			for (uint8_t proxyAxisIdx = 0; proxyAxisIdx < AXESNO; proxyAxisIdx++) {
				float alignment = tinybvh_dot(
						proxyKDopAxes[proxyAxisIdx],
						quadSpherePointVec);

				if (alignment > 0.999f) { axisIdx = proxyAxisIdx; sign = 1; return; }
				if (alignment < -0.999f) { axisIdx = proxyAxisIdx; sign = -1; return; }
			}
		}
	};

	class BVH_DOBB : public BVHBase {
	public:
		// TODO: do compressed BVHNode
		struct ALIGNED(16) BVHNode
		{
			// Block 0
			bvhvec3 nodeMinPoint;
			int8_t exponent[3];
			uint8_t iMask;
			// Block 1
			uint32_t childBaseIndexRotation;	// 0 - 23: childBaseIndex, 24 - 31: rotation
			uint32_t triBaseIndex;
			uint8_t meta[8];
			// Block 2 - 4
			uint8_t qLo[3][8];
			uint8_t qHi[3][8];

			static constexpr uint32_t CHILD_INDEX_BITS = 24;
			static constexpr uint32_t CHILD_INDEX_MASK = (1u << CHILD_INDEX_BITS) - 1u;

			uint32_t ChildBaseIndex() const { return childBaseIndexRotation & CHILD_INDEX_MASK; }
			uint8_t Rotation() const { return (uint8_t)(childBaseIndexRotation >> CHILD_INDEX_BITS); }
			void SetChildBaseIndexRotation(const uint32_t childBaseIndex, const uint8_t rotationIndex) {
				childBaseIndexRotation = (childBaseIndex & CHILD_INDEX_MASK) | ((uint32_t)rotationIndex << CHILD_INDEX_BITS);
			}
		};

		// TODO: move outside of class
		// Only implementing |K| = 13 case
		struct sidebandData
		{
			// Index corresponds to axes
			bvhvec2 kDopExtents[AXESNO];
			uint8_t rotation;
			float AABBSASum;
		};

		BVH_DOBB(BVHContext ctx = {}) { context = ctx; } // TODO: add layout
		~BVH_DOBB();

		void Build(const bvhvec4* vertices, const uint32_t primCount);
		void Build(const bvhvec4slice& vertices);
		void Build(const bvhvec4* vertices, const uint32_t* indices, const uint32_t primCount);
		void Build(const bvhvec4slice& vertices, const uint32_t* indices, const uint32_t primCount);

		// Leaves can have at most two triangles in them!
		void ConvertFrom(const MBVH<8>& original);

		// BVH data
		BVHNode* dobbNode = 0;
		uint32_t* dobbPrimIdx = 0;
		uint32_t allocatedPrimIndices = 0;
		MBVH<8> bvh8;
		bool ownBVH8 = true;

		sidebandData* sideband = 0;
		static DOBB_LUT lut;
};

} // namespace tinybvh
#endif // TINY_BVH_DOBB_H_


#if defined( TINYBVH_DOBB_IMPLEMENTATION ) && !defined( TINY_BVH_DOBB_IMPL_DONE )
#define TINY_BVH_DOBB_IMPL_DONE

namespace tinybvh {
	DOBB_LUT BVH_DOBB::lut;

	BVH_DOBB::~BVH_DOBB()
	{
		if (!ownBVH8) bvh8.ReleaseOwnership();
		AlignedFree(dobbNode);
		AlignedFree(dobbPrimIdx);
		
		AlignedFree(sideband);
	}

	void BVH_DOBB::Build(const bvhvec4* vertices, const uint32_t primCount) {
		Build(bvhvec4slice(vertices, primCount * 3, sizeof(bvhvec4)), 0, 0);
	}

	void BVH_DOBB::Build(const bvhvec4slice& vertices) {
		Build(vertices, 0, 0);
	}

	void BVH_DOBB::Build(const bvhvec4* vertices, const uint32_t* indices, const uint32_t primCount) {
		Build(bvhvec4slice(vertices, primCount * 3, sizeof(bvhvec4)), indices, primCount);
	}

	void BVH_DOBB::Build(const bvhvec4slice& vertices, const uint32_t* indices, const uint32_t primCount) {
		// make sure we own bvh8 (a previous ConvertFrom may have left it referencing someone else's)
		if (!ownBVH8) bvh8.DropReference(context), ownBVH8 = true;
		bvh8.bvh.context = context;
		// Build binary BVH
		bvh8.bvh.Build(vertices, indices, primCount);
		// DOBB leaves can hold only up to 2 triangles
		bvh8.bvh.SplitLeafs(2);
		// Create 8 wide BVH
		bvh8.ConvertFrom(bvh8.bvh, true);
		// Create oriented tree
		ConvertFrom(bvh8);
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
		// ConvertFrom does not handle cases where there are 2 or less primitives
		FATAL_ERROR_IF(bvh8.usedNodes == 0, "BVH_DOBB::ConvertFrom( .. ), source MBVH is empty.");
		FATAL_ERROR_IF(bvh8.mbvhNode[0].isLeaf(), "BVH_DOBB::ConvertFrom( .. ), root is a leaf (scene has 2 triangles or fewer).");
		bvh8.ReferenceFrom(original);
		// get base properties
		CopyBasePropertiesFrom(bvh8);
		// Allocate proxy k-DOPs and rotation sideband data
		AlignedFree(sideband);
		sideband = (sidebandData*)AlignedAlloc(bvh8.usedNodes * sizeof(sidebandData));

		// Count the number of interior nodes in bvh8, DOBB BVH needs this
		// amount of nodes in total
		uint32_t nodesNeeded = 0;

		// Iterate through nodes backwards (children are processed
		// implicitly before parents)
		for (int nodeIdx = bvh8.usedNodes - 1; nodeIdx >= 0; nodeIdx--)
		{
			// Reset proxy k-DOP for min/max selection
			for (int axisIdx = 0; axisIdx < AXESNO; axisIdx++)
			{
				sideband[nodeIdx].kDopExtents[axisIdx] = bvhvec2(BVH_FAR, -BVH_FAR);
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
				float bestSim = -1.0f;
				sideband[nodeIdx].rotation = 0;
				for (uint8_t rotationIdx = 0; rotationIdx < AXESNO * SMALLM * 2; rotationIdx++) {
					float sim = lut.similarity({ a0, a1, a2 }, rotationIdx);
					if (sim > bestSim) {
						bestSim = sim;
						sideband[nodeIdx].rotation = rotationIdx;
					}
				}

				// Create proxy KDop to propagate up
				for (uint32_t vertIdx = 0; vertIdx < noLeafVerts; vertIdx++) {
					for (uint32_t axisIdx = 0; axisIdx < AXESNO; axisIdx++) {
						sideband[nodeIdx].kDopExtents[axisIdx].x = tinybvh_min(
							tinybvh_dot(proxyKDopAxes[axisIdx], leafVerts[vertIdx]),
							sideband[nodeIdx].kDopExtents[axisIdx].x);

						sideband[nodeIdx].kDopExtents[axisIdx].y = tinybvh_max(
							tinybvh_dot(proxyKDopAxes[axisIdx], leafVerts[vertIdx]),
							sideband[nodeIdx].kDopExtents[axisIdx].y);
					}
				}
			}
			else
			{
				// Increment nodes needed for DOBB_BVH
				nodesNeeded++;

				// Process internal nodes
				for (uint32_t childIdx = 0; childIdx < currNode.childCount; childIdx++)
				{
					uint32_t currChildIdx = currNode.child[childIdx];

					// Determine k-DOP extents
					for (uint32_t axisIdx = 0; axisIdx < AXESNO; axisIdx++)
					{
						sideband[nodeIdx].kDopExtents[axisIdx].x = tinybvh_min(
							sideband[currChildIdx].kDopExtents[axisIdx].x,
							sideband[nodeIdx].kDopExtents[axisIdx].x);

						sideband[nodeIdx].kDopExtents[axisIdx].y = tinybvh_max(
							sideband[currChildIdx].kDopExtents[axisIdx].y,
							sideband[nodeIdx].kDopExtents[axisIdx].y);
					}
				}

				float AABBSAMax = 0.0f;
				uint8_t candidateRotationIdx = 0;
				sideband[nodeIdx].AABBSASum = 0.0f;

				// Process internal nodes
				for (uint32_t childIdx = 0; childIdx < currNode.childCount; childIdx++)
				{
					uint32_t currChildIdx = currNode.child[childIdx];

					// Gather extents from euclidean axes from axes collection
					float a = fabsf(sideband[currChildIdx].kDopExtents[0].x - sideband[currChildIdx].kDopExtents[0].y);
					float b = fabsf(sideband[currChildIdx].kDopExtents[1].x - sideband[currChildIdx].kDopExtents[1].y);
					float c = fabsf(sideband[currChildIdx].kDopExtents[2].x - sideband[currChildIdx].kDopExtents[2].y);

					// Select rotation by maxselecting on AABB surface area
					float childAABBSA = 2.0f * (a * b + b * c + c * a);
					sideband[nodeIdx].AABBSASum += childAABBSA;

					if (childAABBSA > AABBSAMax) {
						AABBSAMax = childAABBSA;
						candidateRotationIdx = sideband[currChildIdx].rotation;
					}
				}

				sideband[nodeIdx].rotation = candidateRotationIdx;
			}
		}

		// Allocate nodes and tri indices for DOBB_BVH
		if (allocatedNodes < nodesNeeded) {
			AlignedFree(dobbNode);
			dobbNode = (BVHNode*)AlignedAlloc(nodesNeeded * sizeof(BVHNode));

			allocatedNodes = nodesNeeded;
		}

		if (allocatedPrimIndices < idxCount) {
			AlignedFree(dobbPrimIdx);
			dobbPrimIdx = (uint32_t*)AlignedAlloc(idxCount * sizeof(uint32_t));
		
			allocatedPrimIndices = idxCount;
		}
	
		uint32_t stackBvhNodes[256], stackDobbNodes[256];
		stackBvhNodes[0] = 0; stackDobbNodes[0] = 0;
		uint32_t stackIdx = 1;
		uint32_t dobbNodeIdx = 1, dobbTriIdx = 0;

		while (stackIdx > 0) {
			// Pop top member of stack
			stackIdx--;
			uint32_t currBvhIdx = stackBvhNodes[stackIdx];
			uint32_t currDobbIdx = stackDobbNodes[stackIdx];
			const auto& currBvhNode = bvh8.mbvhNode[currBvhIdx];
			BVHNode& currDobbNode = dobbNode[currDobbIdx];
			memset(&currDobbNode, 0, sizeof(BVHNode));

			uint8_t candidateRotationIdx = sideband[currBvhIdx].rotation;
			DOBB_LUT::Mat3 candidateRotation = lut.rotations[candidateRotationIdx];

			// Collect basis vectors of current rotation
			bvhvec3 rotatedAxes[6] = {
				bvhvec3(candidateRotation.m[0][0], candidateRotation.m[1][0], candidateRotation.m[2][0]),
				bvhvec3(-candidateRotation.m[0][0], -candidateRotation.m[1][0], -candidateRotation.m[2][0]),
				bvhvec3(candidateRotation.m[0][1], candidateRotation.m[1][1], candidateRotation.m[2][1]),
				bvhvec3(-candidateRotation.m[0][1], -candidateRotation.m[1][1], -candidateRotation.m[2][1]),
				bvhvec3(candidateRotation.m[0][2], candidateRotation.m[1][2], candidateRotation.m[2][2]),
				bvhvec3(-candidateRotation.m[0][2], -candidateRotation.m[1][2], -candidateRotation.m[2][2])
			};

			float rotatedAxesExtents[6];

			// Values for e
			float DOBBSASum = 0.0f;
			bvhvec3 childMinPoint[8], childMaxPoint[8];

			for (uint32_t childIdx = 0; childIdx < currBvhNode.childCount; childIdx++)
			{
				uint32_t currChildIdx = currBvhNode.child[childIdx];

				// Do reprojection based on the selected rotation
				for (uint32_t rotatedAxesIdx = 0; rotatedAxesIdx < 6; rotatedAxesIdx++)
				{
					bvhvec3 apexPoint;

					// Do reprojection for the given axis
					// Reprojection entry corresponding to current rotation and axis
					const DOBB_LUT::reprojectionEntry& reprojection = lut.reprojections[candidateRotationIdx][rotatedAxesIdx];
					// Normals of selected planes
					bvhvec3 v0 = proxyKDopAxes[reprojection.axisIdx[0]] * (float)reprojection.sign[0];
					bvhvec3 v1 = proxyKDopAxes[reprojection.axisIdx[1]] * (float)reprojection.sign[1];
					bvhvec3 v2 = proxyKDopAxes[reprojection.axisIdx[2]] * (float)reprojection.sign[2];
					// Extents of selected planes
					float h0, h1, h2;
					if (reprojection.sign[0] == 1)
						h0 = sideband[currChildIdx].kDopExtents[reprojection.axisIdx[0]].y;
					else h0 = -sideband[currChildIdx].kDopExtents[reprojection.axisIdx[0]].x;
					if (reprojection.sign[1] == 1)
						h1 = sideband[currChildIdx].kDopExtents[reprojection.axisIdx[1]].y;
					else h1 = -sideband[currChildIdx].kDopExtents[reprojection.axisIdx[1]].x;
					if (reprojection.sign[2] == 1)
						h2 = sideband[currChildIdx].kDopExtents[reprojection.axisIdx[2]].y;
					else h2 = -sideband[currChildIdx].kDopExtents[reprojection.axisIdx[2]].x;
					// TODO: Remove unnecessary variable declarations after confirming it works correctly
					// Determine apex point by finding intersection of three planes
					// Find vectors colinear with the intersection of pairwise planes
					bvhvec3 c0 = tinybvh_cross(v1, v2);
					bvhvec3 c1 = tinybvh_cross(v2, v0);
					bvhvec3 c2 = tinybvh_cross(v0, v1);
					apexPoint = (c0 * h0 + c1 * h1 + c2 * h2) / tinybvh_dot(v0, c0);
					// Determine axis extent by projecting apex point onto the axis
					rotatedAxesExtents[rotatedAxesIdx] = tinybvh_dot(rotatedAxes[rotatedAxesIdx], apexPoint);
				}

				childMinPoint[childIdx].x = -rotatedAxesExtents[1];
				childMinPoint[childIdx].y = -rotatedAxesExtents[3];
				childMinPoint[childIdx].z = -rotatedAxesExtents[5];

				childMaxPoint[childIdx].x = rotatedAxesExtents[0];
				childMaxPoint[childIdx].y = rotatedAxesExtents[2];
				childMaxPoint[childIdx].z = rotatedAxesExtents[4];

				float a = childMaxPoint[childIdx].x - childMinPoint[childIdx].x;
				float b = childMaxPoint[childIdx].y - childMinPoint[childIdx].y;
				float c = childMaxPoint[childIdx].z - childMinPoint[childIdx].z;
				DOBBSASum += 2 * (a * b + a * c + b * c);
			}

			// Determine chosen rotation, record extents if no rotation
			uint8_t chosenRotation;

			if (DOBBSASum < sideband[currBvhIdx].AABBSASum) {
				chosenRotation = candidateRotationIdx;
			}
			else {
				chosenRotation = 255;

				for (uint32_t childIdx = 0; childIdx < currBvhNode.childCount; childIdx++) {
					uint32_t currChildIdx = currBvhNode.child[childIdx];

					childMinPoint[childIdx] = bvhvec3(
						sideband[currChildIdx].kDopExtents[0].x,
						sideband[currChildIdx].kDopExtents[1].x,
						sideband[currChildIdx].kDopExtents[2].x);
					childMaxPoint[childIdx] = bvhvec3(
						sideband[currChildIdx].kDopExtents[0].y,
						sideband[currChildIdx].kDopExtents[1].y,
						sideband[currChildIdx].kDopExtents[2].y);
				}
			}

			// Determine min and max points of bounding boxes
			bvhvec3 minPoint = { BVH_FAR, BVH_FAR, BVH_FAR };
			bvhvec3 maxPoint = { -BVH_FAR, -BVH_FAR, -BVH_FAR };

			for (uint32_t childIdx = 0; childIdx < currBvhNode.childCount; childIdx++) {
				minPoint = tinybvh_min(minPoint, childMinPoint[childIdx]);
				maxPoint = tinybvh_max(maxPoint, childMaxPoint[childIdx]);
			}

			currDobbNode.nodeMinPoint = minPoint;

			// To be used to calculate qLo-s and qHi-s
			float inverseScale[3];

			// Quantize extent of current bounding box
			for (uint32_t axisIdx = 0; axisIdx < 3; axisIdx++) {
				const float extent = maxPoint[axisIdx] - minPoint[axisIdx];
				// guard against in-plane BBs
				int32_t exponent = extent > 0.0f ? (int32_t)ceilf(log2f(extent / 255.0f)) : -126;
				// clamp into int8_t range
				exponent = tinybvh_clamp(exponent, -126, 127);
				// guard against ceil rounding "extent" below the real extent
				// indirectly by giving exponent too small of a value
				if (exponent < 127 && ldexpf(255.0f, exponent) < extent)
					exponent++;
				currDobbNode.exponent[axisIdx] = (int8_t)exponent;
				// calculate inverseScale here, once
				inverseScale[axisIdx] = ldexpf(1.0f, -exponent);
			}

			uint32_t currDobbNodeChildBaseIdx = 0;
			uint32_t currDobbNodeInternalChildCount = 0;
			uint32_t currDobbNodeTriCount = 0;

			// Fill rest of the data for the compressed nodes
			for (uint32_t childIdx = 0; childIdx < currBvhNode.childCount; childIdx++) {
				// Fill qLo and qHi
				for (uint32_t axisIdx = 0; axisIdx < 3; axisIdx++) {
					const float currQLo = floorf((childMinPoint[childIdx][axisIdx] - minPoint[axisIdx]) * inverseScale[axisIdx]);
					currDobbNode.qLo[axisIdx][childIdx] = (uint8_t)tinybvh_clamp(currQLo, 0.0f, 255.0f);
				
					const float currQHi = ceilf((childMaxPoint[childIdx][axisIdx] - minPoint[axisIdx]) * inverseScale[axisIdx]);
					currDobbNode.qHi[axisIdx][childIdx] = (uint8_t)tinybvh_clamp(currQHi, 0.0f, 255.0f);
				}

				uint32_t currChildIdx = currBvhNode.child[childIdx];
				const auto& currChild = bvh8.mbvhNode[currChildIdx];

				// Pack internal child nodes
				if (!currChild.isLeaf())
				{
					// Set base index if its the first internal child
					// dobbNodeIdx holds the next available dobb node
					if (currDobbNodeInternalChildCount++ == 0)
						currDobbNodeChildBaseIdx = dobbNodeIdx;

					// Add next index pair to stack
					stackBvhNodes[stackIdx] = currChildIdx;
					stackDobbNodes[stackIdx++] = dobbNodeIdx++;

					currDobbNode.iMask |= (uint8_t)(1 << childIdx);
					currDobbNode.meta[childIdx] = (uint8_t)(0b00100000 | (24 + childIdx));
				}
				// Pack leaf child nodes
				else
				{
					// Set base index if it is the first leaf child
					// dobbTriIdx holds the next available triangle
					if (currDobbNodeTriCount == 0)
						currDobbNode.triBaseIndex = dobbTriIdx;

					// Pack meta
					const uint8_t unaryTriCount = currChild.triCount == 1 ? 0b001 : 0b011;
					currDobbNode.meta[childIdx] = (uint8_t)((unaryTriCount << 5) | currDobbNodeTriCount);

					// Increment tri count
					currDobbNodeTriCount += currChild.triCount;

					// Copy triangle indices into dobbPrimIdx in the correct
					// order 
					for (uint32_t triIdx = 0; triIdx < currChild.triCount; triIdx++)
						dobbPrimIdx[dobbTriIdx++] = bvh8.bvh.primIdx[currChild.firstTri + triIdx];
				}
			}

			// Set child base index and rotation in dobb node
			currDobbNode.SetChildBaseIndexRotation(currDobbNodeChildBaseIdx, chosenRotation);
		}

		usedNodes = dobbNodeIdx;

		AlignedFree(sideband);
		sideband = 0;
	}

} // namespace tinybvh
#endif // TINYBVH_DOBB_IMPLEMENTATION