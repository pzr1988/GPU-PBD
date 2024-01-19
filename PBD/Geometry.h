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
  void initInertiaTensor(T rho=1);
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
  //Get/Set the actual number of capsules used
  size_t size() const;
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
  typename thrust::device_vector<Capsule<T>>::iterator begin();
  typename thrust::device_vector<Capsule<T>>::iterator end();
  typename thrust::device_vector<Capsule<T>>::const_iterator begin() const;
  typename thrust::device_vector<Capsule<T>>::const_iterator end() const;
  typename thrust::device_ptr<const Capsule<T>> getCapsules() const;
  typename thrust::device_ptr<Capsule<T>> getCapsules();
 protected:
  thrust::device_vector<Capsule<T>> _capsules;
  size_t _nrCapsule=0;
};

// General function to get AABB
template <template<typename> class Geometry, typename T>
struct AABBGetter;
// AABBGetter for Capsule<T>
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
