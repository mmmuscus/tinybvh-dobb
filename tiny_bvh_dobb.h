#ifndef TINY_BVH_DOBB_H_
#define TINY_BVH_DOBB_H_

#ifndef TINY_BVH_H_
#include "tiny_bvh.h"
#endif

#define AXESNO 13

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

	class BVH_DOBB : public BVHBase {
	public:
		struct BVHNode
		{

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
				// Create proxyKDop to propagate up
				for (uint32_t triIdx = 0; triIdx < currNode.triCount; triIdx++) {
					// Find vertices of current triangle
					const uint32_t prim = bvh8.bvh.primIdx[currNode.firstTri + triIdx];
					uint32_t i0, i1, i2;

					// TODO: Consider using GET_PRIM_INDICES_I0_I1_I2 
					// Indexed
					if (bvh8.bvh.vertIdx) {
						i0 = bvh8.bvh.vertIdx[prim * 3];
						i1 = bvh8.bvh.vertIdx[prim * 3 + 1];
						i2 = bvh8.bvh.vertIdx[prim * 3 + 2];
					}
					else {
						i0 = prim * 3;
						i1 = i0 + 1;
						i2 = i0 + 2;
					}

					const bvhvec3 v0 = bvh8.bvh.verts[i0];
					const bvhvec3 v1 = bvh8.bvh.verts[i1];
					const bvhvec3 v2 = bvh8.bvh.verts[i2];
					
					for (int axisIdx = 0; axisIdx < AXESNO; axisIdx++)
					{
						// Project each vertex of triangle onto current axis
						// float tinybvh_dot( const bvhvec3& a, const bvhvec3& b ); 
						float v0Extent = tinybvh_dot(proxyKDopAxes[axisIdx], v0);
						float v1Extent = tinybvh_dot(proxyKDopAxes[axisIdx], v1);
						float v2Extent = tinybvh_dot(proxyKDopAxes[axisIdx], v2);

						proxyKDop[nodeIdx].extents[axisIdx].x = tinybvh_min(
							v0Extent, proxyKDop[nodeIdx].extents[axisIdx].x);
						proxyKDop[nodeIdx].extents[axisIdx].x = tinybvh_min(
							v1Extent, proxyKDop[nodeIdx].extents[axisIdx].x);
						proxyKDop[nodeIdx].extents[axisIdx].x = tinybvh_min(
							v2Extent, proxyKDop[nodeIdx].extents[axisIdx].x);

						proxyKDop[nodeIdx].extents[axisIdx].y = tinybvh_max(
							v0Extent, proxyKDop[nodeIdx].extents[axisIdx].y);
						proxyKDop[nodeIdx].extents[axisIdx].y = tinybvh_max(
							v1Extent, proxyKDop[nodeIdx].extents[axisIdx].y);
						proxyKDop[nodeIdx].extents[axisIdx].y = tinybvh_max(
							v2Extent, proxyKDop[nodeIdx].extents[axisIdx].y);
					}
				}
			}
			else
			{
				// Process internal nodes
				// Iterate through children
				for (uint32_t childIdx = 0; childIdx < currNode.childCount; childIdx++)
				{
					for (uint32_t axisIdx = 0; axisIdx < AXESNO; axisIdx++)
					{
						proxyKDop[nodeIdx].extents[axisIdx].x = tinybvh_min(
							proxyKDop[currNode.child[childIdx]].extents[axisIdx].x,
							proxyKDop[nodeIdx].extents[axisIdx].x);

						proxyKDop[nodeIdx].extents[axisIdx].y = tinybvh_max(
							proxyKDop[currNode.child[childIdx]].extents[axisIdx].y,
							proxyKDop[nodeIdx].extents[axisIdx].y);
					}
				}
			}
		}
	}

} // namespace tinybvh
#endif // TINYBVH_DOBB_IMPLEMENTATION