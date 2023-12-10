#ifndef CAPSULE_H
#define CAPSULE_H
#include "Utils.h"

namespace GPUPBD
{
// This struct represents a capsule-shaped object
// The capsule's centerline extends from (-_len/2,0,0) to (_len/2,0,0)
// The radius of the capsule is stored in _radius
// The local to global transformation is stored in _trans
// The 3x3 rotational matrix is: _trans.template block<3,3>(0,0) (you can also use macro: ROT(_trans))
// The 3x1 translational vector is: _trans.template block<3,1>(0,3) (you can also use macro: CTR(_trans))
template <typename T>
struct Capsule
{
  DECL_MAT_VEC_MAP_TYPES_T
  T _len, _radius;
  Mat3X4T _trans;
  DEVICE_HOST Vec3T minCorner() const
  {
    return Vec3T(-_len / 2, 0, 0);
  };
  DEVICE_HOST Vec3T maxCorner() const
  {
    return Vec3T(_len / 2, 0, 0);
  };
  DEVICE_HOST Vec3T absoluteMinCorner() const
  {
    return ROT(_trans) * minCorner() + CTR(_trans);
  };
  DEVICE_HOST Vec3T absoluteMaxCorner() const
  {
    return ROT(_trans) * maxCorner() + CTR(_trans);
  };
};
}
#endif