#ifndef COLLISION_CUH
#define COLLISION_CUH
#include "Geometry.h"
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

}
#endif
