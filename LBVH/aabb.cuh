#ifndef LBVH_AABB_CUH
#define LBVH_AABB_CUH
#include <vector_types.h>
#include <thrust/swap.h>
#include <cmath>

namespace lbvh {

template<typename T> struct vector_of;
template<> struct vector_of<float>  {
  using type = float4;
};
template<> struct vector_of<double> {
  using type = double4;
};

template<typename T>
struct aabb {
  typename vector_of<T>::type upper;
  typename vector_of<T>::type lower;
};

template<typename T>
__device__ __host__
inline bool intersects(const aabb<T>& lhs, const aabb<T>& rhs) noexcept {
  if(lhs.upper.x < rhs.lower.x || rhs.upper.x < lhs.lower.x) {
    return false;
  }
  if(lhs.upper.y < rhs.lower.y || rhs.upper.y < lhs.lower.y) {
    return false;
  }
  if(lhs.upper.z < rhs.lower.z || rhs.upper.z < lhs.lower.z) {
    return false;
  }
  return true;
}

__device__ __host__
inline aabb<double> merge(const aabb<double>& lhs, const aabb<double>& rhs) noexcept {
  aabb<double> merged;
  merged.upper.x = ::fmax(lhs.upper.x, rhs.upper.x);
  merged.upper.y = ::fmax(lhs.upper.y, rhs.upper.y);
  merged.upper.z = ::fmax(lhs.upper.z, rhs.upper.z);
  merged.lower.x = ::fmin(lhs.lower.x, rhs.lower.x);
  merged.lower.y = ::fmin(lhs.lower.y, rhs.lower.y);
  merged.lower.z = ::fmin(lhs.lower.z, rhs.lower.z);
  return merged;
}

__device__ __host__
inline aabb<float> merge(const aabb<float>& lhs, const aabb<float>& rhs) noexcept {
  aabb<float> merged;
  merged.upper.x = ::fmaxf(lhs.upper.x, rhs.upper.x);
  merged.upper.y = ::fmaxf(lhs.upper.y, rhs.upper.y);
  merged.upper.z = ::fmaxf(lhs.upper.z, rhs.upper.z);
  merged.lower.x = ::fminf(lhs.lower.x, rhs.lower.x);
  merged.lower.y = ::fminf(lhs.lower.y, rhs.lower.y);
  merged.lower.z = ::fminf(lhs.lower.z, rhs.lower.z);
  return merged;
}



template<typename T>
__device__ __host__
inline typename vector_of<T>::type centroid(const aabb<T>& box) noexcept {
  typename vector_of<T>::type c;
  c.x = (box.upper.x + box.lower.x) * 0.5;
  c.y = (box.upper.y + box.lower.y) * 0.5;
  c.z = (box.upper.z + box.lower.z) * 0.5;
  return c;
}

// 获得物体的bounding box
template<template<typename> class Geometry, typename T>
struct aabb_getter {
  __device__
  lbvh::aabb<float> operator()(const Geometry<T>& c) const noexcept {
    lbvh::aabb<float> retval;
    return retval;
  }
};

} // lbvh
#endif// LBVH_AABB_CUH
