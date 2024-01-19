#include "Collision.h"
#include "NarrowPhase.h"

namespace GPUPBD {
template <typename T>
CollisionDetector<T>::CollisionDetector(std::shared_ptr<Geometry<T>> geometry)
  :_geometry(geometry) {}
template <typename T>
void CollisionDetector<T>::detectCollisions() {
  // TODO yeti, lbvh里的host和device中的capsule是没必要创建的。
  lbvh::bvh<float, Capsule<T>, AABBGetter<Capsule, T>> bvh(_geometry->getCapsules().cbegin(), _geometry->getCapsules().cend(), true);
  const auto bvh_dev = bvh.get_device_repr();
  std::size_t numCapsules = _geometry->getCapsules().size();
  if(_collisionsTemporary.size() < numCapsules * maxCollisionPerObject)
    _collisionsTemporary.resize(numCapsules * maxCollisionPerObject);
  Collision<T>* d_collisionsTemporary = thrust::raw_pointer_cast(_collisionsTemporary.data());

  // First fill all the collisions, including invalid ones
  thrust::for_each(thrust::device,
                   thrust::make_counting_iterator<std::size_t>(0),
                   thrust::make_counting_iterator<std::size_t>(numCapsules),
  [bvh_dev, d_collisionsTemporary] __device__ (std::size_t idx) {
    unsigned int buffer[maxCollisionPerObject];
    const auto& self = bvh_dev.objects[idx];
    // broad phase
    const auto num_found = lbvh::query_device(bvh_dev, lbvh::overlaps(AABBGetter<Capsule, T>()(self)), buffer, maxCollisionPerObject);
    // narrow phase
    Collision<T>* localMemory = d_collisionsTemporary + idx * maxCollisionPerObject;
    ContactManifold<T> contactM(&self, idx, localMemory);
    for (size_t i = 0; i < num_found; i++)
      if(idx < buffer[i]) {
        const auto &rhs = bvh_dev.objects[buffer[i]];
        contactM.UpdateRhs(&rhs, buffer[i]);
        if(contactM._numCollision < maxCollisionPerObject)
          NarrowPhase<T>::narrowPhaseCollision(contactM, maxCollisionPerObject);
      }
    // fill invalid
    assert(contactM._numCollision < maxCollisionPerObject);
    for(size_t i = contactM._numCollision; i < maxCollisionPerObject; i++)
      localMemory[i]._isValid = false;
  });
  cudaDeviceSynchronize();

  // Then remove invalid ones
  if(_collisions.size() < numCapsules * maxCollisionPerObject)
    _collisions.resize(numCapsules * maxCollisionPerObject);
  auto ret = thrust::copy_if(_collisionsTemporary.begin(),
                  _collisionsTemporary.end(),
                  _collisions.begin(),
                  [] __device__(const Collision<T>& c) {
    return c._isValid;
  });
  int nCollision = std::distance(_collisions.begin(), ret);
  _collisions.resize(nCollision);
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
