#ifndef TINY_BVH_DOBB_H_
#define TINY_BVH_DOBB_H_

#ifndef TINY_BVH_H_
#include "tiny_bvh.h"
#endif

namespace tinybvh {

class BVH_DOBB : public BVHBase {
	struct BVHNode
	{
		
	};
};


} // namespace tinybvh
#endif // TINY_BVH_DOBB_H_


#if defined( TINYBVH_DOBB_IMPLEMENTATION ) && !defined( TINY_BVH_DOBB_IMPL_DONE )
#define TINY_BVH_DOBB_IMPL_DONE

namespace tinybvh {


} // namespace tinybvh
#endif // TINYBVH_DOBB_IMPLEMENTATION