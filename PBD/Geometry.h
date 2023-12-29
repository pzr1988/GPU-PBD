#ifndef GEOMETRY_H
#define GEOMETRY_H
#include "Utils.h"
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
  /* State variables */
  Vec3T _x;
  Eigen::Quaternion<T> _q;
  Vec3T _xNext; //Tentative transformation at the next timestep
  Eigen::Quaternion<T> _qNext;
  /* Derived quantities (auxiliary variables) */
  // Mat3T _Iinv; // inverse of inertia tensor
  // Mat3T _R; //rotation matrix
  // Vec3T _v; //linear velocity
  // Vec3T _w; //angular velocity
  /* Computed quantities */
  // Vec3T _force;
  // Vec3T _torque;
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
  DEVICE_HOST void initInertiaTensor() {
    //https://gamedev.net/tutorials/programming/math-and-physics/capsule-inertia-tensor-r3856/
    T mCy = _mass * (_radius*_radius*_len)/(_radius*_radius*_len + (4./3.)*_radius*_radius*_radius);
    T mHSph = _mass * ((2./3.)*_radius*_radius*_radius)/(_radius*_radius*_len + (4./3.)*_radius*_radius*_radius);
    T Ixx = mCy*(_radius*_radius/2.)
            + 2.*mHSph*(2.*_radius*_radius/5.);
    T Iyy = mCy*(_radius*_radius/4. + _len*_len/12.)
            + 2.*mHSph*(2.*_radius*_radius/5.+_len*_len/2.+3.*_len*_radius/8);
    _Ibody << Ixx, 0, 0,
           0, Iyy, 0,
           0, 0, Iyy;
    _Ibodyinv << 1./Ixx, 0, 0,
              0, 1./Iyy, 0,
              0, 0, 1./Iyy;
  }
};

//This structure represents a collision between a pair of capsules
//The index of first/second capsule is store in _capsuleIdA/_capsuleIdB
//The points of contact in global coordinates are stored in _gobalPointA/_gobalPointB
//The separating direction of contact is stored in _globalNormal, extending from A to B and having unit norm
template <typename T>
struct Collision {
  DECL_MAT_VEC_MAP_TYPES_T
  int _capsuleIdA;
  int _capsuleIdB;
  Vec3T _gobalPointA;
  Vec3T _gobalPointB;
  Vec3T _globalNormal;
  bool _isValid;
  DEVICE_HOST Collision():_capsuleIdA(-1),_capsuleIdB(-1),_gobalPointA(Vec3T()),_gobalPointB(Vec3T()),_globalNormal(Vec3T()),_isValid(false) {}
  DEVICE_HOST T depth() const {
    return (_gobalPointA-_gobalPointB).dot(_globalNormal);
  }
};

//The geometry stores a vector of capsules
//The vector is over-sized and pre-allocated
//The number of capsules in use is stored in _nrCapsule
template <typename T>
struct Geometry {
  DECL_MAT_VEC_MAP_TYPES_T
  //Set the actual number of capsules used
  void resize(int nrCapsule);
  //Set the pre-allocated capsule list
  void reserve(int nrCapsule);
  //CPU->GPU transfer: copying the list of capsules to GPU
  void setCapsule(const std::vector<Capsule<T>>& c);
  //Set the id-th capsule
  void setCapsule(int id, const Capsule<T>& c);
  //Get the id-th capsule
  Capsule<T> operator[](int id) const;
  //Get All Capsules
  const thrust::device_vector<Capsule<T>>& getCapsules() const;
 protected:
  thrust::device_vector<Capsule<T>> _capsules;
  int _nrCapsule=0;
};
}


namespace lbvh {
// 获得物体的bounding box
template<template<typename> class Geometry, typename T>
struct aabb_getter {
  __device__
  lbvh::aabb<float> operator()(const Geometry<T>& c) const noexcept {
    lbvh::aabb<float> retval;
    return retval;
  }
};
// 获得胶囊体的bounding box
template<>
struct lbvh::aabb_getter<GPUPBD::Capsule, float> {
  __device__
  lbvh::aabb<float> operator()(const GPUPBD::Capsule<float> &c) const noexcept {
    lbvh::aabb<float> retval;
    auto transformedEnd1 = c.globalMaxCorner();
    auto transformedEnd2 = c.globalMinCorner();
    auto upper = transformedEnd1.cwiseMax(transformedEnd2);
    float radius = static_cast<float>(c._radius);
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
