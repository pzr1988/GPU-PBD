#ifndef CONTACT_GENERATOR_H
#define CONTACT_GENERATOR_H
#include "Pragma.h"
#include "Collision.h"
#include "Capsule.h"

namespace GPUPBD {

template<typename T>
class ContactGenerator {
 public:

  DECL_MAT_VEC_MAP_TYPES_T
  struct ContactManifold {
    uint32_t _lhsId, _rhsId;
    const Capsule<T> *_lhs, *_rhs;
    T _lhsRadius, _rhsRadius;
    Vec3T _lhsMinCorner, _lhsMaxCorner, _rhsMinCorner, _rhsMaxCorner;
    Collision<T> *_localMemory;
    size_t _numCollision;
    DEVICE_HOST ContactManifold(const Capsule<T> *lhs, const Capsule<T> *rhs, uint32_t lhsId, uint32_t rhsId, Collision<T> *localMemory):
      _lhsId(lhsId), _rhsId(rhsId), _lhs(lhs), _rhs(rhs), _lhsRadius(lhs->_radius), _rhsRadius(rhs->_radius),
      _localMemory(localMemory), _numCollision(0) {
      _lhsMinCorner = _lhs->absoluteMinCorner().template cast<T>();
      _lhsMaxCorner = _lhs->absoluteMaxCorner().template cast<T>();
      _rhsMinCorner = _rhs->absoluteMinCorner().template cast<T>();
      _rhsMaxCorner = _rhs->absoluteMaxCorner().template cast<T>();
    }
    DEVICE_HOST ContactManifold(const Capsule<T> *lhs,  uint32_t lhsId, Collision<T> *localMemory):
      _lhsId(lhsId), _lhs(lhs), _lhsRadius(lhs->_radius),
      _localMemory(localMemory), _numCollision(0) {
      _numCollision = 0;
      _lhsRadius = _lhs->_radius;
      _lhsMinCorner = _lhs->absoluteMinCorner().template cast<T>();
      _lhsMaxCorner = _lhs->absoluteMaxCorner().template cast<T>();
    }
    void DEVICE_HOST UpdateRhs(const Capsule<T> *rhs, uint32_t rhsId) {
      _rhsId = rhsId;
      _rhs = rhs;
      _rhsRadius = _rhs->_radius;
      _rhsMinCorner = _rhs->absoluteMinCorner().template cast<T>();
      _rhsMaxCorner = _rhs->absoluteMaxCorner().template cast<T>();
    }
  };


  // 胶囊体碰撞
  // maxCollisionsPerNode是单个胶囊体碰撞的最大数目
  static DEVICE_HOST int narrowPhaseCollision(ContactManifold& contactM, size_t maxCollisionsPerNode) noexcept;

 private:
  // 胶囊体碰撞内部使用
  static DEVICE_HOST int generateManifoldSphereSphereInternal(ContactManifold& contactM, size_t maxCollisionsPerNode);
  static DEVICE_HOST void generateManifoldSphereCapsuleInternal(GPUPBD::Collision<T> &originCollision,
      const Vec3T &cA, T cARadius, int lhs_idx,
      const Vec3T &cB1, const Vec3T &cB2, T cBRadius, int rhs_idx);
};
}

#endif
