#ifndef GEOMETRY_H
#define GEOMETRY_H
#include "Pragma.h"
#include "LBVH/aabb.cuh"
#include <thrust/device_vector.h>

namespace GPUPBD {
//This struct represents a capsule-shaped object
//The capsule's centerline extends from (-_len/2,0,0) to (_len/2,0,0)
//The radius of the capsule is stored in _radius
//The local to global transformation is stored in _x and _q
//The 3x3 rotational matrix is: _q.toRotationMatrix()
//The 3x1 translational vector is: _x
template <typename T>
struct Capsule {
  DECL_MAT_VEC_MAP_TYPES_T
  /* Constant quantities */
  T _len,_radius;
  T _mass;
  Mat3T _Ibody, _Ibodyinv; //body space inertia tensor
  bool _isDynamic;
  /* State variables */
  Vec3T _x;
  QuatT _q;
  Vec3T _xPrev; //Auxilliary variable for XPBD
  QuatT _qPrev;
  Vec3T _v; //linear velocity
  Vec3T _w; //angular velocity
  /* Derived quantities (auxiliary variables) */
  Mat3T _Iinv; // inverse of inertia tensor
  Mat3T _R; //rotation matrix
  /* Computed quantities */
  Vec3T _force;
  Vec3T _torque;
  DEVICE_HOST Vec3T minCorner() const {
    return Vec3T(-_len / 2, 0, 0);
  }
  DEVICE_HOST Vec3T maxCorner() const {
    return Vec3T(_len / 2, 0, 0);
  }
  DEVICE_HOST Vec3T globalMinCorner() const {
    return _q.toRotationMatrix()*minCorner()+_x;
  }
  DEVICE_HOST Vec3T globalMaxCorner() const {
    return _q.toRotationMatrix()*maxCorner()+_x;
  }
  DEVICE_HOST void initInertiaTensor(T rho=1) {
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
  DEVICE_HOST Mat3T getInertiaTensorInv() const {
    auto R = _q.toRotationMatrix();
    return R*_Ibodyinv*R.transpose();
  }
};

//The geometry stores a vector of capsules
//The vector is over-sized and pre-allocated
//The number of capsules in use is stored in _nrCapsule
template <typename T>
struct Geometry {
  DECL_MAT_VEC_MAP_TYPES_T
  //Set the actual number of capsules used
  void resize(size_t nrCapsule);
  //Set the pre-allocated capsule list
  void reserve(size_t nrCapsule);
  //CPU->GPU transfer: copying the list of capsules to GPU
  void setCapsule(const std::vector<Capsule<T>>& c);
  //Set the id-th capsule
  void setCapsule(size_t id, const Capsule<T>& c);
  //Get the id-th capsule
  Capsule<T> operator[](size_t id) const;
  //Get All Capsules
  const thrust::device_vector<Capsule<T>>& getCapsules() const;
  // Get mutable capsules for integration, don't change the number of capsuels.
  thrust::device_vector<Capsule<T>>& getMutableCapsules();
 protected:
  thrust::device_vector<Capsule<T>> _capsules;
  size_t _nrCapsule=0;
};

// 获得物体的bounding box
template <template<typename> class Geometry, typename T>
struct AABBGetter;
// 获得胶囊体的bounding box
template <typename T>
struct AABBGetter<Capsule, T> {
  DEVICE_HOST lbvh::aabb<T> operator()(const Capsule<T> &c) const noexcept {
    lbvh::aabb<T> retval;
    auto transformedEnd1 = c.globalMaxCorner();
    auto transformedEnd2 = c.globalMinCorner();
    auto upper = transformedEnd1.cwiseMax(transformedEnd2);
    T radius = static_cast<T>(c._radius);
    retval.upper.x = upper.x() + radius;
    retval.upper.y = upper.y() + radius;
    retval.upper.z = upper.z() + radius;
    auto lower = transformedEnd1.cwiseMin(transformedEnd2);
    retval.lower.x = lower.x() - radius;
    retval.lower.y = lower.y() - radius;
    retval.lower.z = lower.z() - radius;
    return retval;
  }
};
}
#endif
