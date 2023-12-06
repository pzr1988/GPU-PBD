#ifndef LBVH_PREDICATOR_CUH
#define LBVH_PREDICATOR_CUH
#include "aabb.cuh"

namespace lbvh
{

template<template<typename> class Geometry, typename Real>
struct query_overlap
{
  __device__ __host__
  query_overlap(const Geometry<Real>& tgt, std::uint32_t index):origin(tgt),
    bbox(aabb_getter<Geometry, Real>()(tgt)), object_idx(index) {}

  query_overlap()  = default;
  ~query_overlap() = default;
  query_overlap(const query_overlap&) = default;
  query_overlap(query_overlap&&)      = default;
  query_overlap& operator=(const query_overlap&) = default;
  query_overlap& operator=(query_overlap&&)      = default;

  __device__ __host__
  inline bool operator()(const aabb<float>& box) noexcept
  {
    return intersects(box, bbox);
  }
  Geometry<Real> origin;
  aabb<float> bbox;
  std::uint32_t object_idx;
};

template<template<typename> class Geometry, typename Real>
__device__ __host__
query_overlap<Geometry, Real> overlaps(const Geometry<Real>& region, std::uint32_t idx) noexcept
{
  return query_overlap<Geometry, Real>(region, idx);
}

} // lbvh
#endif// LBVH_PREDICATOR_CUH
