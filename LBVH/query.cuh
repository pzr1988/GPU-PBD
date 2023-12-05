#ifndef LBVH_QUERY_CUH
#define LBVH_QUERY_CUH
#include "predicator.cuh"

namespace lbvh {
// query object indices that potentially overlaps with query aabb.
//
// requirements:
//
template<typename Real, typename Objects, bool IsConst,
         template<typename> class QueryObjects,
         template<typename> class Collision>
__device__
unsigned int query_device(
  const detail::basic_device_bvh<Real, Objects, IsConst>& bvh,
  const query_overlap<QueryObjects, Real> q,
  size_t maxCollisionsPerNode,
  Collision<Real>* localMemory) noexcept {
  using bvh_type   = detail::basic_device_bvh<Real, Objects, IsConst>;
  using index_type = typename bvh_type::index_type;
  using aabb_type  = typename bvh_type::aabb_type;
  using node_type  = typename bvh_type::node_type;

  index_type  stack[64]; // is it okay?
  index_type* stack_ptr = stack;
  *stack_ptr++ = 0; // root node is always 0

  size_t num_found = 0;
  do {
    const index_type node  = *--stack_ptr;
    const index_type L_idx = bvh.nodes[node].left_idx;
    const index_type R_idx = bvh.nodes[node].right_idx;

    if(intersects(q.bbox, bvh.aabbs[L_idx])) { //包围盒是否碰撞
      const auto obj_idx = bvh.nodes[L_idx].object_idx;
      if(obj_idx != 0xFFFFFFFF) {
        // 胶囊体碰撞
        if(obj_idx > q.object_idx && num_found < maxCollisionsPerNode) {
          int numCollision = narrowPhaseCollision(
                               q.origin, q.object_idx,
                               bvh.objects[obj_idx], obj_idx,
                               localMemory, num_found, maxCollisionsPerNode);
        }
      } else { // the node is not a leaf.
        *stack_ptr++ = L_idx;
      }
    }
    if(intersects(q.bbox, bvh.aabbs[R_idx])) { //包围盒是否碰撞
      const auto obj_idx = bvh.nodes[R_idx].object_idx;
      if(obj_idx != 0xFFFFFFFF) {
        // 胶囊体碰撞
        if(obj_idx > q.object_idx && num_found < maxCollisionsPerNode) {
          int numCollision = narrowPhaseCollision(
                              q.origin, q.object_idx,
                              bvh.objects[obj_idx], obj_idx,
                              localMemory, num_found, maxCollisionsPerNode);
        }
      } else { // the node is not a leaf.
        *stack_ptr++ = R_idx;
      }
    }
  } while (stack < stack_ptr);
  return num_found;
}
} // lbvh
#endif// LBVH_QUERY_CUH
