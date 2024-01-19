#ifndef COLLISION_H
#define COLLISION_H

#include "Geometry.h"
#include <LBVH/lbvh.cuh>

namespace GPUPBD {

//This structure represents a collision between a pair of capsules
//The index of first/second capsule is store in _capsuleIdA/_capsuleIdB
//The points of contact in local coordinates are stored in _localPointA/_localPointB
//The separating direction of contact is stored in _globalNormal, extending from A to B and having unit norm
template <typename T>
struct Collision {
  DECL_MAT_VEC_MAP_TYPES_T
  bool _isValid;
  int _capsuleIdA;
  int _capsuleIdB;
  Vec3T _localPointA;
  Vec3T _localPointB;
  Vec3T _globalNormal;
  T _alpha; // compliance of the constraint
  DEVICE_HOST Collision():_capsuleIdA(-1),_capsuleIdB(-1),
    _localPointA(Vec3T()),_localPointB(Vec3T()),_globalNormal(Vec3T()),
    _alpha(.0001), _isValid(false) {}
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
  //Get All Collisions
  const thrust::device_vector<Collision<T>>& getCollisions() const;
 protected:
  std::shared_ptr<Geometry<T>> _geometry;
  thrust::device_vector<Collision<T>> _collisions;
};

}
#endif
