#ifndef CONTACT_GENERATOR_H
#define CONTACT_GENERATOR_H
#include "Pragma.h"
#include "Geometry.h"

namespace GPUPBD {

template<typename T>
class ContactGenerator {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  struct ContactManifold {
    uint32_t _lhsId, _rhsId;
    const Capsule<T> *_lhs,*_rhs;
    T _lhsRadius,_rhsRadius;
    Vec3T _lhsMinCorner,_lhsMaxCorner,_rhsMinCorner,_rhsMaxCorner;
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
  static DEVICE_HOST int narrowPhaseCollision(ContactManifold& contactM,size_t maxCollisionsPerNode) noexcept {
    typedef Eigen::Matrix<T, 3, 2> Mat3X2T;
    int numCollision = 0;
    Vec3T &cA1 = contactM._lhsMinCorner;
    Vec3T &cA2 = contactM._lhsMaxCorner;
    Vec3T &cB1 = contactM._rhsMinCorner;
    Vec3T &cB2 = contactM._rhsMaxCorner;
    Vec3T nA = cA2 - cA1, nB = cB2 - cB1;
    T nLenASqr = nA.squaredNorm(), nLenA = sqrt((double)nLenASqr);
    T nLenBSqr = nB.squaredNorm(), nLenB = sqrt((double)nLenBSqr);
    nA /= nLenA;
    nB /= nLenB;

    if (abs(nA.dot(nB)) > 1 - epsDir) {
      // nearly parallel 产生1~2个碰撞点
      T dB1 = (cB1 - cA1).dot(nA);
      T dB2 = (cB2 - cA1).dot(nA);
      if (dB1 <= 0 && dB2 <= 0) {
        // sphere-sphere 产生0/1个碰撞点
        generateManifoldSphereSphereInternal(contactM, maxCollisionsPerNode);
      } else if (dB1 >= nLenA && dB2 >= nLenA) {
        // sphere-sphere 产生0/1个碰撞点
        generateManifoldSphereSphereInternal(contactM, maxCollisionsPerNode);
      } else if (maxCollisionsPerNode - contactM._numCollision > 2) {
        // 假如localMemory不夠就忽略了
        // range 产生0/2个碰撞点
        Vec3T dir = cB1 - cA1 - dB1 * nA;
        T distSqr = dir.squaredNorm(), dist = 0;
        T sumRad = contactM._lhsRadius + contactM._rhsRadius, sumRadSqr = sumRad * sumRad;
        // not in contact
        if (distSqr > sumRadSqr)
          return numCollision;
        // in contact
        Vec3T nA2B;
        if (distSqr > epsDist * epsDist) {
          nA2B = dir;
          nA2B /= dist = sqrt((double)distSqr);
        } else {
          typename Vec3T::Index id;
          nA.cwiseAbs().minCoeff(&id);
          nA2B = nA.cross(Vec3T::Unit(id));
          nA2B /= nA2B.template cast<double>().norm();
        }
        // two contacts
        Vec2T range(std::max<T>(0, std::min(dB1, dB2)), std::min(nLenA, std::max(dB1, dB2)));
        for (const T &r : range) {
          contactM._localMemory[contactM._numCollision]._capsuleIdA = contactM._lhsId;
          contactM._localMemory[contactM._numCollision]._capsuleIdB = contactM._rhsId;
          contactM._localMemory[contactM._numCollision]._gobalPointA = cA1 + nA * r + nA2B * contactM._lhsRadius;
          contactM._localMemory[contactM._numCollision]._gobalPointB = cA1 + nA * r + dir - nA2B * contactM._rhsRadius;
          contactM._localMemory[contactM._numCollision]._globalNormal = nA2B;
          contactM._localMemory[contactM._numCollision]._isValid = true;
          contactM._numCollision++;
          numCollision++;
        }
      }
    } else {
      // not parallel
      Mat3X2T LHS;
      LHS.col(0) = -(cA2 - cA1);
      LHS.col(1) = (cB2 - cB1);
      Vec2T bary = (LHS.transpose() * LHS).inverse() * (LHS.transpose() * (cA1 - cB1));
      if ((bary.array() >= 0).all() && (bary.array() <= 1).all()) {
        Vec3T cA = cA1 * (1 - bary[0]) + cA2 * bary[0];
        Vec3T cB = cB1 * (1 - bary[1]) + cB2 * bary[1];
        T distSqr = (cA - cB).squaredNorm();
        T sumRad = contactM._lhsRadius + contactM._rhsRadius, sumRadSqr = sumRad * sumRad;
        if (distSqr > sumRadSqr)
          return numCollision;
        Vec3T nA2B = nA.cross(nB);
        nA2B /= nA2B.template cast<double>().norm();
        if ((cB - cA).dot(nA2B) < 0)
          nA2B *= -1;
        contactM._localMemory[contactM._numCollision]._capsuleIdA = contactM._lhsId;
        contactM._localMemory[contactM._numCollision]._capsuleIdB = contactM._rhsId;
        contactM._localMemory[contactM._numCollision]._gobalPointA = cA + nA2B * contactM._lhsRadius;
        contactM._localMemory[contactM._numCollision]._gobalPointB = cB - nA2B * contactM._rhsRadius;
        contactM._localMemory[contactM._numCollision]._globalNormal = nA2B;
        contactM._localMemory[contactM._numCollision]._isValid = true;
        contactM._numCollision++;
        numCollision++;
      } else {
        GPUPBD::Collision<T> collision;
        generateManifoldSphereCapsuleInternal(collision,
                                              cA1, contactM._lhsRadius, contactM._lhsId,
                                              cB1, cB2, contactM._rhsRadius, contactM._rhsId);
        generateManifoldSphereCapsuleInternal(collision,
                                              cA2, contactM._lhsRadius, contactM._lhsId,
                                              cB1, cB2, contactM._rhsRadius, contactM._rhsId);
        generateManifoldSphereCapsuleInternal(collision,
                                              cB1, contactM._rhsRadius, contactM._rhsId,
                                              cA1, cA2, contactM._lhsRadius, contactM._lhsId);
        generateManifoldSphereCapsuleInternal(collision,
                                              cB2, contactM._rhsRadius, contactM._rhsId,
                                              cA1, cA2, contactM._lhsRadius, contactM._lhsId);
        if (collision._isValid) {
          contactM._localMemory[contactM._numCollision]._capsuleIdA = collision._capsuleIdA;
          contactM._localMemory[contactM._numCollision]._capsuleIdB = collision._capsuleIdB;
          contactM._localMemory[contactM._numCollision]._gobalPointA = collision._gobalPointA;
          contactM._localMemory[contactM._numCollision]._gobalPointB = collision._gobalPointB;
          contactM._localMemory[contactM._numCollision]._globalNormal = collision._globalNormal;
          contactM._localMemory[contactM._numCollision]._isValid = true;
          contactM._numCollision++;
          numCollision++;
        }
      }
    }
    return numCollision;
  }
 private:
  // 胶囊体碰撞内部使用
  static DEVICE_HOST int generateManifoldSphereSphereInternal(ContactManifold& contactM,size_t maxCollisionsPerNode) {
    Vec3T &cA1 = contactM._lhsMinCorner;
    Vec3T &cA2 = contactM._lhsMaxCorner;
    Vec3T &cB1 = contactM._rhsMinCorner;
    Vec3T &cB2 = contactM._rhsMaxCorner;

    typedef Eigen::Matrix<T, 4, 1> Vec4T;
    Vec4T distSqrs((cA1 - cB1).squaredNorm(), (cA2 - cB1).squaredNorm(), (cA1 - cB2).squaredNorm(), (cA2 - cB2).squaredNorm());
    typename Vec4T::Index id;
    T distSqr = distSqrs.minCoeff(&id), dist = 0;
    T sumRad = contactM._lhsRadius + contactM._rhsRadius, sumRadSqr = sumRad * sumRad;
    const Vec3T &cA = (id % 2 == 0) ? cA1 : cA2;
    const Vec3T &cB = (id < 2) ? cB1 : cB2;
    // not in contact
    if (distSqr > sumRadSqr)
      return 0;
    // in contact
    Vec3T nA2B;
    if (distSqr > epsDist * epsDist) {
      // a single contact point
      nA2B = cB - cA;
      nA2B /= dist = sqrt((double)distSqr);
    } else {
      // overlapping degenerate case
      distSqrs.maxCoeff(&id);
      const Vec3T &cAn = (id % 2 == 0) ? cA1 : cA2;
      const Vec3T &cBn = (id < 2) ? cB1 : cB2;
      nA2B = cBn - cAn;
      nA2B /= nA2B.template cast<double>().norm();
    }
    contactM._localMemory[contactM._numCollision]._capsuleIdA = contactM._lhsId;
    contactM._localMemory[contactM._numCollision]._capsuleIdB = contactM._rhsId;
    contactM._localMemory[contactM._numCollision]._gobalPointA = cA + nA2B * contactM._lhsRadius;
    contactM._localMemory[contactM._numCollision]._gobalPointB = cB - nA2B * contactM._rhsRadius;
    contactM._localMemory[contactM._numCollision]._globalNormal = nA2B;
    contactM._localMemory[contactM._numCollision]._isValid = true;
    contactM._numCollision++;
    return 1;
  }
  static DEVICE_HOST void generateManifoldSphereCapsuleInternal(GPUPBD::Collision<T>& originCollision,
      const Vec3T &cA,T cARadius,int lhs_idx,
      const Vec3T &cB1,const Vec3T &cB2,T cBRadius,int rhs_idx) {

    GPUPBD::Collision<T> collision;
    Vec3T n = cB2 - cB1;
    T nLenSqr = n.squaredNorm(), nLen = sqrt((double)nLenSqr);
    n /= nLen;
    T d = (cA - cB1).dot(n);
    T sumRad = cARadius + cBRadius, sumRadSqr = sumRad * sumRad;
    // three cases
    if (d <= 0) {
      T distSqr = (cA - cB1).squaredNorm(), dist = 0;
      // not in contact
      if (distSqr > sumRadSqr)
        return ;
      // in contact
      if (distSqr > epsDist * epsDist) {
        // a single contact point
        collision._globalNormal = cB1 - cA;
        collision._globalNormal /= dist = sqrt((double)distSqr);
      } else {
        collision._globalNormal = n;
      }
      collision._gobalPointA = cA + collision._globalNormal * cARadius;
      collision._gobalPointB = cB1 - collision._globalNormal * cBRadius;
    } else if (d >= nLen) {
      T distSqr = (cA - cB2).squaredNorm(), dist = 0;
      // not in contact
      if (distSqr > sumRadSqr)
        return;
      // in contact
      if (distSqr > epsDist * epsDist) {
        // a single contact point
        collision._globalNormal = cB2 - cA;
        collision._globalNormal /= dist = sqrt((double)distSqr);
      } else {
        collision._globalNormal = -n;
      }
      collision._gobalPointA = cA + collision._globalNormal * cARadius;
      collision._gobalPointB = cB2 - collision._globalNormal * cBRadius;
    } else if (d > 0 && d < nLen) {
      Vec3T dir = cA - cB1 - n * d;
      T distSqr = dir.squaredNorm(), dist = 0;
      // not in contact
      if (distSqr > sumRadSqr)
        return;
      // in contact
      if (distSqr > epsDist * epsDist) {
        collision._globalNormal = -dir;
        collision._globalNormal /= dist = sqrt((double)distSqr);
      } else {
        typename Vec3T::Index id;
        n.cwiseAbs().minCoeff(&id);
        collision._globalNormal = n.cross(Vec3T::Unit(id));
        collision._globalNormal /= collision._globalNormal.template cast<double>().norm();
      }
      collision._gobalPointA = cA + collision._globalNormal * cARadius;
      collision._gobalPointB = cA - dir - collision._globalNormal * cBRadius;
    }
    if (!originCollision._isValid || originCollision.depth() < collision.depth()) {
      originCollision._capsuleIdA = lhs_idx;
      originCollision._capsuleIdB = rhs_idx;
      originCollision._gobalPointA = collision._gobalPointA;
      originCollision._gobalPointB = collision._gobalPointB;
      originCollision._globalNormal = collision._globalNormal;
      originCollision._isValid = true;
    }
  }
};
}

#endif