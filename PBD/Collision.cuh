#ifndef COLLISION_CUH
#define COLLISION_CUH
#include "Geometry.cuh"
#include <LBVH/lbvh.cuh>

namespace GPUPBD {
//The collisionDetector has the capability of detecting all pairs of collisions between all pairs of capsules
//The main function (detectCollisions()) is supposed to fill up the vector _collisions with active collision constraints
template <typename T>
class CollisionDetector {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  //The contructor takes the set of capsules
  CollisionDetector(std::shared_ptr<Geometry<T>> geometry);
  //The main detection function that fill up the vector _collisions
  void detectCollisions();
  //Fetch the id-th collision
  Collision<T> operator[](int id);
  //Return the number of detected collisions
  int size() const;
 protected:
  std::shared_ptr<Geometry<T>> _geometry;
  thrust::device_vector<Collision<T>> _collisions;
 private:
  // 在中间过程使用的，碰撞点集合 
  thrust::device_vector<Collision<T>> _globalMemory;
};

template <typename T>
CollisionDetector<T>::CollisionDetector(std::shared_ptr<Geometry<T>> geometry)
  :_geometry(geometry) {}
template <typename T>
void CollisionDetector<T>::detectCollisions() {
  // bvh 内部都是float类型
  auto bIter = _geometry->getCapsules().cbegin();
  auto eIter = _geometry->getCapsules().cend();
  // TODO yeti, lbvh里的host和device中的capsule是没必要创建的。
  lbvh::bvh<float, GPUPBD::Capsule<T>, lbvh::aabb_getter<GPUPBD::Capsule, T>> bvh(bIter, eIter, true);
  const auto bvh_dev = bvh.get_device_repr();
  std::cout << "testing query_device:overlap ...\n";
  const int maxCollisionsPerNode = 4;
  std::size_t numCapsules = _geometry->getCapsules().size();
  std::size_t maxCollisions = maxCollisionsPerNode * numCapsules;
  _globalMemory.clear();
  _globalMemory.resize(maxCollisions);
  thrust::fill(thrust::device, _globalMemory.begin(), _globalMemory.end(), Collision<T>());
  Collision<T>* collisionsPtr = thrust::raw_pointer_cast(_globalMemory.data());
  thrust::for_each(thrust::device,
      thrust::make_counting_iterator<std::size_t>(0),
      thrust::make_counting_iterator<std::size_t>(numCapsules),
      [bvh_dev, maxCollisionsPerNode, collisionsPtr] __device__ (std::size_t idx) {
      Collision<T> localMemory[maxCollisionsPerNode];
      const auto self = bvh_dev.objects[idx];
      const auto num_found = lbvh::query_device(
          bvh_dev, lbvh::overlaps(self, idx),
          maxCollisionsPerNode,
          localMemory);
      for (size_t i = 0; i < num_found; i++) {
          size_t globalIndex = idx * maxCollisionsPerNode + i;
          collisionsPtr[globalIndex] = localMemory[i];
      }
      return;
  });
  // 等待设备上所有先前的 kernel 完成
  // TODO yeti 耗时问题
  cudaError_t cudaStatus = cudaDeviceSynchronize();
  if (cudaStatus != cudaSuccess) {
    // TODO yeti 处理错误
    std::cerr << "CUDA Synchronization error: " << cudaGetErrorString(cudaStatus) << std::endl;
  }
  size_t numValidCollisions = thrust::count_if(thrust::device, _globalMemory.begin(), _globalMemory.end(),
    [] __device__ (const Collision<T>& collision) {
    return collision._isValid;
  });
  _collisions.resize(numValidCollisions);
  thrust::copy_if(thrust::device, _globalMemory.begin(), _globalMemory.end(),
      _collisions.begin(),
      [] __device__ (const Collision<T>& collision) {
      return collision._isValid;
  });
}

template <typename T>
Collision<T> CollisionDetector<T>::operator[](int id) {
  return _collisions[id];
}
template <typename T>
int CollisionDetector<T>::size() const {
  return _collisions.size();
}

}
#endif
