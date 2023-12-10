#include "Collision.h"
#include "Capsule.h"
#include "LBVH/lbvh.cuh"
#include "ContactGenerator.h"

namespace GPUPBD
{
template <typename T>
CollisionDetector<T>::CollisionDetector(std::shared_ptr<Geometry<T>> geometry) : _geometry(geometry) {}
template <typename T>
void CollisionDetector<T>::detectCollisions()
{
  // bvh 内部都是float类型
  auto bIter = _geometry->getCapsules().cbegin();
  auto eIter = _geometry->getCapsules().cend();
  // TODO yeti, lbvh里的host和device中的capsule是没必要创建的。
  lbvh::bvh<float, GPUPBD::Capsule<T>, lbvh::aabb_getter<GPUPBD::Capsule, T>> bvh(bIter, eIter, true);
  const auto bvh_dev = bvh.get_device_repr();
  const int maxCollisionsPerNode = 4;
  std::size_t numCapsules = _geometry->getCapsules().size();
  std::size_t maxCollisions = maxCollisionsPerNode * numCapsules;
  _globalMemory.clear();
  _globalMemory.resize(maxCollisions);
  thrust::fill(thrust::device, _globalMemory.begin(), _globalMemory.end(), Collision<T>());
  Collision<T> *collisionsPtr = thrust::raw_pointer_cast(_globalMemory.data());
  thrust::for_each(thrust::device,
                   thrust::make_counting_iterator<std::size_t>(0),
                   thrust::make_counting_iterator<std::size_t>(numCapsules),
                   [bvh_dev, maxCollisionsPerNode, numCapsules, collisionsPtr] __device__(std::size_t idx)
  {
    Collision<T> localMemory[maxCollisionsPerNode];
    const auto self = bvh_dev.objects[idx];
    const auto num_found = lbvh::query_device(
                             bvh_dev, lbvh::overlaps(self, idx),
                             maxCollisionsPerNode,
                             localMemory);
    for (size_t i = 0; i < num_found; i++)
    {
      size_t globalIndex = idx * maxCollisionsPerNode + i;
      collisionsPtr[globalIndex] = localMemory[i];
    }
    return;
  });
  // 等待设备上所有先前的 kernel 完成
  // TODO yeti 耗时问题
  cudaError_t cudaStatus = cudaDeviceSynchronize();
  if (cudaStatus != cudaSuccess)
  {
    // TODO yeti 处理错误
    std::cerr << "CUDA Synchronization error: " << cudaGetErrorString(cudaStatus) << std::endl;
  }
  size_t numValidCollisions = thrust::count_if(thrust::device, _globalMemory.begin(), _globalMemory.end(),
                              [] __device__(const Collision<T> &collision)
  {
    return collision._isValid;
  });
  _collisions.resize(numValidCollisions);
  thrust::copy_if(thrust::device, _globalMemory.begin(), _globalMemory.end(),
                  _collisions.begin(),
                  [] __device__(const Collision<T> &collision)
  {
    return collision._isValid;
  });
}

template <typename T>
void CollisionDetector<T>::detectCollisionsDebug()
{
  _collisions.clear();
  size_t numCapsules = _geometry->getCapsules().size();
  const int maxCollisionsPerNode = 4;
  for (int i = 0; i < numCapsules; i++)
  {
    Collision<T> localMemory[maxCollisionsPerNode];
    auto lhs = (*_geometry)[i];
    typename ContactGenerator<T>::ContactManifold contactM(&lhs, i, localMemory);
    for (int j = i + 1; j < numCapsules; j++)
    {
      auto rhs = (*_geometry)[j];
      contactM.UpdateRhs(&rhs, j);
      if(contactM._numCollision < maxCollisionsPerNode) {
        int numCollision = ContactGenerator<T>::narrowPhaseCollision(
                             contactM, maxCollisionsPerNode);
      }
    }
    for (size_t k = 0; k < contactM._numCollision; k++)
    {
      _collisions.push_back(contactM._localMemory[k]);
    }
  }
}

template <typename T>
Collision<T> CollisionDetector<T>::operator[](int id)
{
  return _collisions[id];
}
template <typename T>
int CollisionDetector<T>::size() const
{
  return _collisions.size();
}

// declare instance
template struct CollisionDetector<LSCALAR>;
}
