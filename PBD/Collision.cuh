#ifndef COLLISION_CUH
#define COLLISION_CUH
#include "Geometry.cuh"
#include <LBVH/lbvh.cuh>

namespace GPUPBD {
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
};
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
};

template <typename T>
CollisionDetector<T>::CollisionDetector(std::shared_ptr<Geometry<T>> geometry)
  :_geometry(geometry) {}
template <typename T>
void CollisionDetector<T>::detectCollisions() {
  //To see how to implement this function, refer to the tutorial:
  //http://media.steampowered.com/apps/valve/2015/DirkGregorius_Contacts.pdf
  //To see a sample code of how this can be implemented on CPU, see:
  //https://bitbucket.org/runningblade/libdifferentiable/src/ConvexHullPBAD/Environment/ContactGenerator.h

  // bvh 内部都是float类型
  auto bIter = _geometry->getCapsules().cbegin();
  auto eIter = _geometry->getCapsules().cend();
  // TODO yeti, lbvh里的host和device中的capsule是没必要创建的。
  lbvh::bvh<float, GPUPBD::Capsule<T>, aabb_getter<T>> bvh(bIter, eIter, true);
  const auto bvh_dev = bvh.get_device_repr();
  int N = 10;
  std::cout << "testing query_device:overlap ...\n";
  thrust::for_each(thrust::device,
      thrust::make_counting_iterator<std::size_t>(0),
      thrust::make_counting_iterator<std::size_t>(N),
      [bvh_dev] __device__ (std::size_t idx) {
      unsigned int buffer[10];
      const auto self = bvh_dev.objects[idx];
      const float  dr = 0.1f;
      for(std::size_t i=1; i<10; ++i)
      {
          for(unsigned int j=0; j<10; ++j)
          {
              buffer[j] = 0xFFFFFFFF;
          }
          const float r = dr * i;
          lbvh::aabb<float> query_box;
          query_box.lower = make_float4(self._trans(0,3)-r, self._trans(1,3)-r, self._trans(2,3)-r, 0);
          query_box.upper = make_float4(self._trans(0,3)+r, self._trans(1,3)+r, self._trans(2,3)+r, 0);
          const auto num_found = lbvh::query_device(
              bvh_dev, lbvh::overlaps(query_box), buffer, 10);
      }
      return ;
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
