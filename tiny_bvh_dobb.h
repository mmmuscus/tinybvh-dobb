#ifndef TINY_BVH_DOBB_H_
#define TINY_BVH_DOBB_H_

#ifndef TINY_BVH_H_
#include "tiny_bvh.h"
#endif

namespace tinybvh {

	class BVH_DOBB : public BVHBase {
	public:
		struct BVHNode
		{

		};

		BVH_DOBB(BVHContext ctx = {}) { context = ctx; } // TODO: add layout
		~BVH_DOBB();
		void ConvertFrom(const MBVH<8>& original);

		// BVH data
		BVHNode* dobbNode = 0;
		MBVH<8> bvh8;
		bool ownBVH8 = true;
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
			AlignedFree(dobbNode);
			dobbNode = (BVHNode*)AlignedAlloc(nodesNeeded * sizeof(BVHNode));
			allocatedNodes = nodesNeeded;
		}
		usedNodes = nodesNeeded;
		// TODO: Implement DOBB conversion
	}

} // namespace tinybvh
#endif // TINYBVH_DOBB_IMPLEMENTATION