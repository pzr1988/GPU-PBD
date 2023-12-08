#include "Geometry.h"

namespace GPUPBD
{
template <typename T>
void Geometry<T>::resize(int nrCapsule)
{
  _nrCapsule = nrCapsule;
  if (_nrCapsule < (int)_capsules.size())
    _capsules.resize(_nrCapsule);
}
template <typename T>
void Geometry<T>::reserve(int nrCapsule)
{
  _capsules.resize(std::max<int>(nrCapsule, _nrCapsule));
}
template <typename T>
void Geometry<T>::setCapsule(const std::vector<Capsule<T>> &c)
{
  // ASSERT_MSG((int)c.size()==_nrCapsule,"Incorrect number of capsules on host!")
  thrust::copy(c.begin(), c.end(), _capsules.begin());
}
template <typename T>
void Geometry<T>::setCapsule(int id, const Capsule<T> &c)
{
  _capsules[id] = c;
}
template <typename T>
Capsule<T> Geometry<T>::operator[](int id) const
{
  return _capsules[id];
}
template <typename T>
const thrust::device_vector<Capsule<T>> &Geometry<T>::getCapsules() const
{
  return _capsules;
}

// declare instance
template struct Geometry<LSCALAR>;
}
