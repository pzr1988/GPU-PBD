#ifndef GEOMETRY_H
#define GEOMETRY_H
#include "Utils.h"
#include "LBVH/aabb.cuh"
#include <thrust/device_vector.h>

namespace GPUPBD {
//This struct represents a capsule-shaped object
//The capsule's centerline extends from (-_len/2,0,0) to (_len/2,0,0)
//The radius of the capsule is stored in _radius
//The local to global transformation is stored in _trans
//The 3x3 rotational matrix is: _trans.template block<3,3>(0,0) (you can also use macro: ROT(_trans))
//The 3x1 translational vector is: _trans.template block<3,1>(0,3) (you can also use macro: CTR(_trans))
template <typename T>
struct Capsule {
  DECL_MAT_VEC_MAP_TYPES_T
  T _len,_radius;
  Mat3X4T _trans;
  Mat3X4T _transNext;   //Tentative transformation at the next timestep
  Vec3T _v; //linear velocity
  Vec3T _w; //angular velocity
  DEVICE_HOST Vec3T minCorner() const {
    return Vec3T(-_len / 2, 0, 0);
  }
  DEVICE_HOST Vec3T maxCorner() const {
    return Vec3T(_len / 2, 0, 0);
  }
  DEVICE_HOST Vec3T absoluteMinCorner() const {
    return ROT(_trans)*minCorner()+CTR(_trans);
  }
  DEVICE_HOST Vec3T absoluteMaxCorner() const {
    return ROT(_trans)*maxCorner()+CTR(_trans);
  }
};

//This structure represents a collision between a pair of capsules
//The index of first/second capsule is store in _capsuleIdA/_capsuleIdB
//The points of contact in local coordinates are stored in _localPointA/_localPointB
//The separating direction of contact is stored in _globalNormal, extending from A to B and having unit norm
template <typename T>
struct Collision {
  DECL_MAT_VEC_MAP_TYPES_T
  int _capsuleIdA;
  int _capsuleIdB;
  Vec3T _localPointA;
  Vec3T _localPointB;
  Vec3T _globalNormal;
  bool _isValid;
  DEVICE_HOST Collision():_capsuleIdA(-1),_capsuleIdB(-1),_localPointA(Vec3T()),_localPointB(Vec3T()),_globalNormal(Vec3T()),_isValid(false) {}
  DEVICE_HOST T depth() const {
    return (_localPointA-_localPointB).dot(_globalNormal);
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
    Eigen::Matrix<float, 4, 1> end1(static_cast<float>(c._len)/2.0, 0, 0, 1); // 第一个端点
    Eigen::Matrix<float, 4, 1> end2(-static_cast<float>(c._len)/2.0, 0, 0, 1); // 第二个端点

    Eigen::Matrix<float, 3, 1> transformedEnd1 = c._trans * end1;
    Eigen::Matrix<float, 3, 1> transformedEnd2 = c._trans * end2;

    Eigen::Matrix<float, 3, 1> upper = transformedEnd1.head<3>().cwiseMax(transformedEnd2.head<3>());
    float radius = static_cast<float>(c._radius);
    retval.upper.x = upper.x() + radius;
    retval.upper.y = upper.y() + radius;
    retval.upper.z = upper.z() + radius;
    Eigen::Matrix<float, 3, 1> lower = transformedEnd1.head<3>().cwiseMin(transformedEnd2.head<3>()) ;
    retval.lower.x = lower.x() - radius;
    retval.lower.y = lower.y() - radius;
    retval.lower.z = lower.z() - radius;
    return retval;
  }
};

}
#endif
