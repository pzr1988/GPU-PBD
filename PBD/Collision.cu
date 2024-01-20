#include "Collision.h"
#include "NarrowPhase.h"

namespace GPUPBD {
template <typename T>
CollisionDetector<T>::CollisionDetector(std::shared_ptr<Geometry<T>> geometry)
  :_geometry(geometry) {}
template <typename T>
void CollisionDetector<T>::detectCollisions() {
  // TODO yeti, no need to create Capsule in lbvh
  lbvh::bvh<float, Capsule<T>, AABBGetter<Capsule, T>> bvh(_geometry->begin(), _geometry->end(), true);
  const auto bvh_dev = bvh.get_device_repr();
  std::size_t numCapsules = _geometry->size();
  if(_collisionsTemporary.size() < numCapsules * maxCollisionPerObject)
    _collisionsTemporary.resize(numCapsules * maxCollisionPerObject);
  Constraint<T>* d_collisionsTemporary = thrust::raw_pointer_cast(_collisionsTemporary.data());

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
    Constraint<T>* localMemory = d_collisionsTemporary + idx * maxCollisionPerObject;
    ContactManifold<T> contactM(&self, idx, localMemory);
    for (size_t i = 0; i < num_found; i++)
      if(idx < buffer[i]) {
        const auto &rhs = bvh_dev.objects[buffer[i]];
        contactM.updateRhs(&rhs, buffer[i]);
        if(contactM.canCollide())
          NarrowPhase<T>::narrowPhaseCollision(contactM, maxCollisionPerObject);
      }
    // fill invalid
    assert(contactM._numCollision < maxCollisionPerObject);
    for(size_t i = 0; i < maxCollisionPerObject; i++)
      localMemory[i]._type = Collision;
    for(size_t i = contactM._numCollision; i < maxCollisionPerObject; i++)
      localMemory[i]._isValid = false;
  });
  cudaDeviceSynchronize();

  // Second remove invalid ones
  if(_collisions.size() < numCapsules * maxCollisionPerObject)
    _collisions.resize(numCapsules * maxCollisionPerObject);
  auto collisionsEnd = thrust::copy_if(_collisionsTemporary.begin(), _collisionsTemporary.end(), _collisions.begin(), [] __device__(const Constraint<T>& c) {return c._isValid;});
  _size = std::distance(_collisions.begin(), collisionsEnd);
}
template <typename T>
Constraint<T> CollisionDetector<T>::operator[](int id) {
  return _collisions[id];
}
template <typename T>
size_t CollisionDetector<T>::size() const {
  return _size;
}
template <typename T>
typename thrust::device_vector<Constraint<T>>::iterator CollisionDetector<T>::begin() {
  return _collisions.begin();
}
template <typename T>
typename thrust::device_vector<Constraint<T>>::iterator CollisionDetector<T>::end() {
  return _collisions.begin()+_size;
}
template <typename T>
typename thrust::device_vector<Constraint<T>>::const_iterator CollisionDetector<T>::begin() const {
  return _collisions.begin();
}
template <typename T>
typename thrust::device_vector<Constraint<T>>::const_iterator CollisionDetector<T>::end() const {
  return _collisions.begin()+_size;
}
template <typename T>
typename thrust::device_ptr<const Constraint<T>> CollisionDetector<T>::getCollisions() const {
  return _collisions.data();
}
template <typename T>
typename thrust::device_ptr<Constraint<T>> CollisionDetector<T>::getCollisions() {
  return _collisions.data();
}

//declare instance
template struct CollisionDetector<LSCALAR>;
}
