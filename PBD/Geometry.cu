#include "Geometry.h"

namespace GPUPBD {
template <typename T>
void Capsule<T>::initInertiaTensor(T rho) {
  //https://gamedev.net/tutorials/programming/math-and-physics/capsule-inertia-tensor-r3856/
  _mass = 3.14f * (_radius*_radius*_len) + 3.14f * (4.f/3.f)*_radius*_radius*_radius;
  _mass *= rho;
  T mCy = _mass * (_radius*_radius*_len)/(_radius*_radius*_len + (4.f/3.f)*_radius*_radius*_radius);
  T mHSph = _mass * ((2.f/3.f)*_radius*_radius*_radius)/(_radius*_radius*_len + (4.f/3.f)*_radius*_radius*_radius);
  T Ixx = mCy*(_radius*_radius/2.f)
          + 2.f*mHSph*(2.f*_radius*_radius/5.f);
  T Iyy = mCy*(_radius*_radius/4.f + _len*_len/12.f)
          + 2.f*mHSph*(2.f*_radius*_radius/5.f+_len*_len/2.f+3.f*_len*_radius/8);
  _Ibody.setIdentity();
  _Ibody.diagonal()=Vec3T(Ixx,Iyy,Iyy);
  _Ibodyinv=_Ibody.inverse();
}

template <typename T>
size_t Geometry<T>::size() const {
  return _nrCapsule;
}
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
typename thrust::device_vector<Capsule<T>>::iterator Geometry<T>::begin() {
  return _capsules.begin();
}
template <typename T>
typename thrust::device_vector<Capsule<T>>::iterator Geometry<T>::end() {
  return _capsules.begin()+_nrCapsule;
}
template <typename T>
typename thrust::device_vector<Capsule<T>>::const_iterator Geometry<T>::begin() const {
  return _capsules.begin();
}
template <typename T>
typename thrust::device_vector<Capsule<T>>::const_iterator Geometry<T>::end() const {
  return _capsules.begin()+_nrCapsule;
}
template <typename T>
typename thrust::device_ptr<const Capsule<T>> Geometry<T>::getCapsules() const {
  return _capsules.data();
}
template <typename T>
typename thrust::device_ptr<Capsule<T>> Geometry<T>::getCapsules() {
  return _capsules.data();
}

//declare instance
template struct Capsule<LSCALAR>;
template struct Geometry<LSCALAR>;
}
