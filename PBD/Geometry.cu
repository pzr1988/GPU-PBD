#include "Geometry.h"

namespace GPUPBD {
template <typename T>
void Geometry<T>::resize(size_t nrCapsule) {
  _nrCapsule=nrCapsule;
  if(_nrCapsule<_capsules.size())
    _capsules.resize(_nrCapsule);
}
template <typename T>
void Geometry<T>::reserve(size_t nrCapsule) {
  _capsules.resize(std::max(nrCapsule,_nrCapsule));
}
template <typename T>
void Geometry<T>::setCapsule(const std::vector<Capsule<T>>& c) {
  assert(c.size()==_nrCapsule);
  thrust::copy(c.begin(),c.end(),_capsules.begin());
}
template <typename T>
void Geometry<T>::setCapsule(size_t id, const Capsule<T>& c) {
  _capsules[id]=c;
}
template <typename T>
Capsule<T> Geometry<T>::operator[](size_t id) const {
  return _capsules[id];
}
template <typename T>
const thrust::device_vector<Capsule<T>>& Geometry<T>::getCapsules() const {
  return _capsules;
}
template <typename T>
thrust::device_vector<Capsule<T>>& Geometry<T>::getMutableCapsules() {
  return _capsules;
}

//declare instance
template struct Geometry<LSCALAR>;
}
