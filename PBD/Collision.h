#ifndef COLLISION_H
#define COLLISION_H

#include "Geometry.h"
#include <LBVH/lbvh.cuh>

namespace GPUPBD {

//This structure represents a constraint between a pair of shapes
//The index of first/second shape is store in _shapeIdA/_shapeIdB
//The points of contact in local coordinates are stored in _localPointA/_localPointB
//The separating direction of contact is stored in _globalNormal, extending from A to B and having unit norm
enum class ConstraintType {
  Collision,
  JointPosition,
  JointAngular,
  Unknown,
};
template <typename T>
struct Constraint {
  DECL_MAT_VEC_MAP_TYPES_T
  ConstraintType _type=ConstraintType::Unknown;
  bool _isValid;
  int _shapeIdA;
  int _shapeIdB; // B is the parent of A
  Vec3T _localPointA;
  Vec3T _localPointB;
  Vec3T _globalNormal;
  QuatT _targetQ; // the rotation from B to A
  QuatT _aQ; // the rotation of A in it's world
  QuatT _bQ; // the rotation of B in it's world
  Vec3T _axis; // rotate axis of angular constraint
  T _theta; // angular of rotation.
  T _alpha; // compliance of the constraint
  DEVICE_HOST Constraint():_shapeIdA(-1),_shapeIdB(-1), _localPointA(Vec3T()),_localPointB(Vec3T()),_globalNormal(Vec3T()),_targetQ(QuatT::Identity()),_aQ(QuatT::Identity()),_bQ(QuatT::Identity()),_alpha(.0001f), _isValid(false) {}
};

template <typename T>
struct ContactManifold {
  DECL_MAT_VEC_MAP_TYPES_T
  unsigned int _lhsId, _rhsId;
  const Shape<T> *_lhs,*_rhs;
  Vec3T _lhsMinCorner,_lhsMaxCorner;
  Vec3T _rhsMinCorner,_rhsMaxCorner;
  Constraint<T> *_localMemory;
  unsigned int _numCollision;
  DEVICE_HOST ContactManifold(const Shape<T> *lhs, const Shape<T> *rhs, unsigned int lhsId, unsigned int rhsId, Constraint<T> *localMemory)
    :_lhsId(lhsId), _rhsId(rhsId), _lhs(lhs), _rhs(rhs), _localMemory(localMemory), _numCollision(0) {
    _lhsMinCorner = _lhs->globalMinCorner().template cast<T>();
    _lhsMaxCorner = _lhs->globalMaxCorner().template cast<T>();
    _rhsMinCorner = _rhs->globalMinCorner().template cast<T>();
    _rhsMaxCorner = _rhs->globalMaxCorner().template cast<T>();
  }
  DEVICE_HOST ContactManifold(const Shape<T> *lhs, unsigned int lhsId, Constraint<T> *localMemory)
    :_lhsId(lhsId), _lhs(lhs), _localMemory(localMemory), _numCollision(0) {
    _lhsMinCorner = _lhs->globalMinCorner().template cast<T>();
    _lhsMaxCorner = _lhs->globalMaxCorner().template cast<T>();
  }
  DEVICE_HOST void updateRhs(const Shape<T> *rhs, unsigned int rhsId) {
    _rhsId = rhsId;
    _rhs = rhs;
    _rhsMinCorner = _rhs->globalMinCorner().template cast<T>();
    _rhsMaxCorner = _rhs->globalMaxCorner().template cast<T>();
  }
  DEVICE_HOST bool canCollide(bool closeSelfCollision=false) const {
    if(!closeSelfCollision)
      return _numCollision < maxCollisionPerObject && _lhs->_parent != _rhsId && _rhs->_parent != _lhsId;
    else
      return _numCollision < maxCollisionPerObject && _lhs->_groupid != _rhs->_groupid;
  }
  DEVICE_HOST void swap() {
    unsigned int tmpId = _lhsId;
    _lhsId=_rhsId, _rhsId=tmpId;
    const Shape<T> *tmpP = _lhs;
    _lhs=_rhs, _rhs=tmpP;
    Vec3T tmpV = _lhsMinCorner;
    _lhsMinCorner=_rhsMinCorner, _rhsMinCorner=tmpV;
    tmpV = _lhsMaxCorner;
    _lhsMaxCorner=_rhsMaxCorner, _rhsMaxCorner=tmpV;
  }
};

//The collisionDetector has the capability of detecting all pairs of collisions between all pairs of shapes
//The main function (detectCollisions()) is supposed to fill up the vector _collisions with active collision constraints
template <typename T>
class CollisionDetector {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  //The contructor takes the set of shapes
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
