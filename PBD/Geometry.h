#ifndef GEOMETRY_CUH
#define GEOMETRY_CUH
#include "Utils.h"
#include <thrust/device_vector.h>
#include "Capsule.h"

namespace GPUPBD
{
// The geometry stores a vector of capsules
// The vector is over-sized and pre-allocated
// The number of capsules in use is stored in _nrCapsule
template <typename T>
struct Geometry
{
  DECL_MAT_VEC_MAP_TYPES_T
  // Set the actual number of capsules used
  void resize(int nrCapsule);
  // Set the pre-allocated capsule list
  void reserve(int nrCapsule);
  // CPU->GPU transfer: copying the list of capsules to GPU
  void setCapsule(const std::vector<Capsule<T>> &c);
  // Set the id-th capsule
  void setCapsule(int id, const Capsule<T> &c);
  // Get the id-th capsule
  Capsule<T> operator[](int id) const;
  // Get All Capsules
  const thrust::device_vector<Capsule<T>> &getCapsules() const;

 protected:
  thrust::device_vector<Capsule<T>> _capsules;
  int _nrCapsule = 0;
};
}

#endif
