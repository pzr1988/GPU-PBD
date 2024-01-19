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
  DEVICE_HOST Collision():_capsuleIdA(-1),_capsuleIdB(-1), _localPointA(Vec3T()),_localPointB(Vec3T()),_globalNormal(Vec3T()), _alpha(.0001), _isValid(false) {}
};

template <typename T>
struct ContactManifold {
  DECL_MAT_VEC_MAP_TYPES_T
  unsigned int _lhsId, _rhsId;
  const Capsule<T> *_lhs,*_rhs;
  Vec3T _lhsMinCorner,_lhsMaxCorner;
  Vec3T _rhsMinCorner,_rhsMaxCorner;
  Collision<T> *_localMemory;
  unsigned int _numCollision;
  DEVICE_HOST ContactManifold(const Capsule<T> *lhs, const Capsule<T> *rhs, unsigned int lhsId, unsigned int rhsId, Collision<T> *localMemory)
    :_lhsId(lhsId), _rhsId(rhsId), _lhs(lhs), _rhs(rhs), _localMemory(localMemory), _numCollision(0) {
    _lhsMinCorner = _lhs->globalMinCorner().template cast<T>();
    _lhsMaxCorner = _lhs->globalMaxCorner().template cast<T>();
    _rhsMinCorner = _rhs->globalMinCorner().template cast<T>();
    _rhsMaxCorner = _rhs->globalMaxCorner().template cast<T>();
  }
  DEVICE_HOST ContactManifold(const Capsule<T> *lhs, unsigned int lhsId, Collision<T> *localMemory)
    :_lhsId(lhsId), _lhs(lhs), _localMemory(localMemory), _numCollision(0) {
    _numCollision = 0;
    _lhsMinCorner = _lhs->globalMinCorner().template cast<T>();
    _lhsMaxCorner = _lhs->globalMaxCorner().template cast<T>();
  }
  DEVICE_HOST void UpdateRhs(const Capsule<T> *rhs, unsigned int rhsId) {
    _rhsId = rhsId;
    _rhs = rhs;
    _rhsMinCorner = _rhs->globalMinCorner().template cast<T>();
    _rhsMaxCorner = _rhs->globalMaxCorner().template cast<T>();
  }
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
  thrust::device_vector<Collision<T>> _collisionsTemporary;
};

}
#endif
