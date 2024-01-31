#ifndef COLLISION_H
#define COLLISION_H

#include "Geometry.h"
#include <LBVH/lbvh.cuh>

namespace GPUPBD {

//This structure represents a constraint between a pair of capsules
//The index of first/second capsule is store in _capsuleIdA/_capsuleIdB
//The points of contact in local coordinates are stored in _localPointA/_localPointB
//The separating direction of contact is stored in _globalNormal, extending from A to B and having unit norm
enum ConstraintType {
  Collision,
  JointPosition,
  JointAngular,
  Unknown,
};
template <typename T>
struct Constraint {
  DECL_MAT_VEC_MAP_TYPES_T
  ConstraintType _type=Unknown;
  bool _isValid;
  int _capsuleIdA;
  int _capsuleIdB;
  Vec3T _localPointA;
  Vec3T _localPointB;
  Vec3T _globalNormal;
  T _alpha; // compliance of the constraint
  DEVICE_HOST Constraint():_capsuleIdA(-1),_capsuleIdB(-1), _localPointA(Vec3T()),_localPointB(Vec3T()),_globalNormal(Vec3T()), _alpha(.0001f), _isValid(false) {}
};

template <typename T>
struct ContactManifold {
  DECL_MAT_VEC_MAP_TYPES_T
  unsigned int _lhsId, _rhsId;
  const Capsule<T> *_lhs,*_rhs;
  Vec3T _lhsMinCorner,_lhsMaxCorner;
  Vec3T _rhsMinCorner,_rhsMaxCorner;
  Constraint<T> *_localMemory;
  unsigned int _numCollision;
  DEVICE_HOST ContactManifold(const Capsule<T> *lhs, const Capsule<T> *rhs, unsigned int lhsId, unsigned int rhsId, Constraint<T> *localMemory)
    :_lhsId(lhsId), _rhsId(rhsId), _lhs(lhs), _rhs(rhs), _localMemory(localMemory), _numCollision(0) {
    _lhsMinCorner = _lhs->globalMinCorner().template cast<T>();
    _lhsMaxCorner = _lhs->globalMaxCorner().template cast<T>();
    _rhsMinCorner = _rhs->globalMinCorner().template cast<T>();
    _rhsMaxCorner = _rhs->globalMaxCorner().template cast<T>();
  }
  DEVICE_HOST ContactManifold(const Capsule<T> *lhs, unsigned int lhsId, Constraint<T> *localMemory)
    :_lhsId(lhsId), _lhs(lhs), _localMemory(localMemory), _numCollision(0) {
    _lhsMinCorner = _lhs->globalMinCorner().template cast<T>();
    _lhsMaxCorner = _lhs->globalMaxCorner().template cast<T>();
  }
  DEVICE_HOST void updateRhs(const Capsule<T> *rhs, unsigned int rhsId) {
    _rhsId = rhsId;
    _rhs = rhs;
    _rhsMinCorner = _rhs->globalMinCorner().template cast<T>();
    _rhsMaxCorner = _rhs->globalMaxCorner().template cast<T>();
  }
  DEVICE_HOST bool canCollide() const {
    return _numCollision < maxCollisionPerObject && _lhs->_parent != _rhsId && _rhs->_parent != _lhsId;
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
  Constraint<T> operator[](int id);
  //Return the number of detected collisions
  size_t size() const;
  //Get All Collisions
  typename thrust::device_vector<Constraint<T>>::iterator begin();
  typename thrust::device_vector<Constraint<T>>::iterator end();
  typename thrust::device_vector<Constraint<T>>::const_iterator begin() const;
  typename thrust::device_vector<Constraint<T>>::const_iterator end() const;
  typename thrust::device_ptr<const Constraint<T>> getCollisions() const;
  typename thrust::device_ptr<Constraint<T>> getCollisions();
 protected:
  size_t _size=0;
  std::shared_ptr<Geometry<T>> _geometry;
  thrust::device_vector<Constraint<T>> _collisions;
  thrust::device_vector<Constraint<T>> _collisionsTemporary;
};

}
#endif
