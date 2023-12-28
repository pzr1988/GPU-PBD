#include "Collision.h"
#include "ContactGenerator.h"

namespace GPUPBD {
template <typename T>
CollisionDetector<T>::CollisionDetector(std::shared_ptr<Geometry<T>> geometry):_geometry(geometry) {}
template <typename T>
void CollisionDetector<T>::detectCollisions() {
  // bvh 内部都是float类型
  auto bIter = _geometry->getCapsules().cbegin();
  auto eIter = _geometry->getCapsules().cend();
  // TODO yeti, lbvh里的host和device中的capsule是没必要创建的。
  lbvh::bvh<float, GPUPBD::Capsule<T>, lbvh::aabb_getter<GPUPBD::Capsule, T>> bvh(bIter, eIter, true);
  const auto bvh_dev = bvh.get_device_repr();
  const int maxCollisionsPerNode = 10;
  std::size_t numCapsules = _geometry->getCapsules().size();
  thrust::device_vector<unsigned int> numCollisionPerNode(numCapsules+1, 0);
  unsigned int* d_numCollisionPerNode = thrust::raw_pointer_cast(numCollisionPerNode.data());

  // First compute the number of total collision, to resize the globalMemory
  thrust::for_each(thrust::device,
                   thrust::make_counting_iterator<std::size_t>(0),
                   thrust::make_counting_iterator<std::size_t>(numCapsules),
  [bvh_dev, maxCollisionsPerNode, d_numCollisionPerNode] __device__ (std::size_t idx) {
    unsigned int buffer[maxCollisionsPerNode];
    const auto & self = bvh_dev.objects[idx];
    // broad phase
    const auto num_found = lbvh::query_device(
                             bvh_dev, lbvh::overlaps(lbvh::aabb_getter<GPUPBD::Capsule, T>()(self)),
                             buffer, maxCollisionsPerNode);

    // narrow phase
    Collision<T> localMemory[maxCollisionsPerNode];
    typename ContactGenerator<T>::ContactManifold contactM(&self, idx, localMemory);
    for (size_t i = 0; i < num_found; i++) {
      if(idx < buffer[i]) {
        const auto &rhs = bvh_dev.objects[buffer[i]];
        contactM.UpdateRhs(&rhs, buffer[i]);
        if(contactM._numCollision < maxCollisionsPerNode) {
          ContactGenerator<T>::narrowPhaseCollision(
            contactM, maxCollisionsPerNode);
        }
      }
    }
    d_numCollisionPerNode[idx] = contactM._numCollision;
    return;
  });

  cudaError_t cudaStatus = cudaDeviceSynchronize();
  thrust::exclusive_scan(numCollisionPerNode.begin(), numCollisionPerNode.end(), numCollisionPerNode.begin());
  _collisions.resize(numCollisionPerNode[numCapsules]);

  // Second fill collisions to the correct position
  Collision<T>* d_collisions = thrust::raw_pointer_cast(_collisions.data());
  thrust::for_each(thrust::device,
                   thrust::make_counting_iterator<std::size_t>(0),
                   thrust::make_counting_iterator<std::size_t>(numCapsules),
  [bvh_dev, maxCollisionsPerNode, d_numCollisionPerNode, d_collisions] __device__ (std::size_t idx) {
    unsigned int buffer[maxCollisionsPerNode];
    const auto & self = bvh_dev.objects[idx];
    // broad phase
    const auto num_found = lbvh::query_device(
                             bvh_dev, lbvh::overlaps(lbvh::aabb_getter<GPUPBD::Capsule, T>()(self)),
                             buffer, maxCollisionsPerNode);

    // narrow phase
    Collision<T> localMemory[maxCollisionsPerNode];
    typename ContactGenerator<T>::ContactManifold contactM(&self, idx, localMemory);
    for (size_t i = 0; i < num_found; i++) {
      if(idx < buffer[i]) {
        const auto & rhs = bvh_dev.objects[buffer[i]];
        contactM.UpdateRhs(&rhs, buffer[i]);
        if(contactM._numCollision < maxCollisionsPerNode) {
          ContactGenerator<T>::narrowPhaseCollision(
            contactM, maxCollisionsPerNode);
        }
      }
    }
    // fill
    for (size_t i = 0; i < contactM._numCollision; i++) {
      d_collisions[d_numCollisionPerNode[idx]+i] = contactM._localMemory[i];
    }
    return;
  });
  cudaStatus = cudaDeviceSynchronize();
}
template <typename T>
Collision<T> CollisionDetector<T>::operator[](int id) {
  return _collisions[id];
}
template <typename T>
int CollisionDetector<T>::size() const {
  return _collisions.size();
}
template <typename T>
const thrust::device_vector<Collision<T>>& CollisionDetector<T>::getCollisions() const {
  return _collisions;
}

//declare instance
template struct CollisionDetector<LSCALAR>;
}
