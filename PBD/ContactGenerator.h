#ifndef CONTACT_GENERATOR_H
#define CONTACT_GENERATOR_H
#include "Pragma.h"
#include "Collision.h"

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
      _lhsMinCorner = _lhs->globalMinCorner().template cast<T>();
      _lhsMaxCorner = _lhs->globalMaxCorner().template cast<T>();
      _rhsMinCorner = _rhs->globalMinCorner().template cast<T>();
      _rhsMaxCorner = _rhs->globalMaxCorner().template cast<T>();
    }
    DEVICE_HOST ContactManifold(const Capsule<T> *lhs,  uint32_t lhsId, Collision<T> *localMemory):
      _lhsId(lhsId), _lhs(lhs), _lhsRadius(lhs->_radius),
      _localMemory(localMemory), _numCollision(0) {
      _numCollision = 0;
      _lhsRadius = _lhs->_radius;
      _lhsMinCorner = _lhs->globalMinCorner().template cast<T>();
      _lhsMaxCorner = _lhs->globalMaxCorner().template cast<T>();
    }
    void DEVICE_HOST UpdateRhs(const Capsule<T> *rhs, uint32_t rhsId) {
      _rhsId = rhsId;
      _rhs = rhs;
      _rhsRadius = _rhs->_radius;
      _rhsMinCorner = _rhs->globalMinCorner().template cast<T>();
      _rhsMaxCorner = _rhs->globalMaxCorner().template cast<T>();
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
      T dB1 = (cA1 - cB1).dot(nB);
      T dB2 = (cA2 - cB1).dot(nB);
      if (dB1 <= 0 && dB2 <= 0) {
        // sphere-sphere 产生0/1个碰撞点
        generateManifoldSphereSphereInternal(contactM, maxCollisionsPerNode);
      } else if (dB1 >= nLenB && dB2 >= nLenB) {
        // sphere-sphere 产生0/1个碰撞点
        generateManifoldSphereSphereInternal(contactM, maxCollisionsPerNode);
      } else if (maxCollisionsPerNode - contactM._numCollision > 2) {
        // 假如localMemory不夠就忽略了
        // range 产生0/2个碰撞点
        Vec3T dir = cB1 - cA1 + dB1 * nB;
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
        T range[2] = {std::max<T>(0, std::min(dB1, dB2)), std::min(nLenB, std::max(dB1, dB2))};
        for (int i=0; i < 2; i++) {
          T r = range[i];
          contactM._localMemory[contactM._numCollision]._capsuleIdA = contactM._lhsId;
          contactM._localMemory[contactM._numCollision]._capsuleIdB = contactM._rhsId;
          T f = dB1 > dB2 ? 1.0 : -1.0;
          Vec3T globalPointA = cA1 + nA * (dB1 -r) * f + nA2B * contactM._lhsRadius;
          Vec3T globalPointB = cB1 + nB * r - nA2B * contactM._rhsRadius;
          contactM._localMemory[contactM._numCollision]._localPointA =
            contactM._lhs->_q.conjugate().toRotationMatrix()*(globalPointA-contactM._lhs->_x);
          contactM._localMemory[contactM._numCollision]._localPointB =
            contactM._rhs->_q.conjugate().toRotationMatrix()*(globalPointB-contactM._rhs->_x);
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
        auto globalPointA = cA + nA2B * contactM._lhsRadius;
        auto globalPointB = cB - nA2B * contactM._rhsRadius;
        contactM._localMemory[contactM._numCollision]._localPointA =
          contactM._lhs->_q.conjugate().toRotationMatrix()*(globalPointA-contactM._lhs->_x);
        contactM._localMemory[contactM._numCollision]._localPointB =
          contactM._rhs->_q.conjugate().toRotationMatrix()*(globalPointB-contactM._rhs->_x);
        contactM._localMemory[contactM._numCollision]._globalNormal = nA2B;
        contactM._localMemory[contactM._numCollision]._isValid = true;
        contactM._numCollision++;
        numCollision++;
      } else {
        GPUPBD::Collision<T> collision;
        T collisionDepth = 0;
        generateManifoldSphereCapsuleInternal(collision, collisionDepth,
                                              cA1, contactM._lhs, contactM._lhsId,
                                              contactM._rhs, contactM._rhsId);
        generateManifoldSphereCapsuleInternal(collision,collisionDepth,
                                              cA2, contactM._lhs, contactM._lhsId,
                                              contactM._rhs, contactM._rhsId);
        generateManifoldSphereCapsuleInternal(collision, collisionDepth,
                                              cB1, contactM._rhs, contactM._rhsId,
                                              contactM._lhs, contactM._lhsId);
        generateManifoldSphereCapsuleInternal(collision, collisionDepth,
                                              cB2, contactM._rhs, contactM._rhsId,
                                              contactM._lhs, contactM._lhsId);
        if (collision._isValid) {
          contactM._localMemory[contactM._numCollision]._capsuleIdA = collision._capsuleIdA;
          contactM._localMemory[contactM._numCollision]._capsuleIdB = collision._capsuleIdB;
          contactM._localMemory[contactM._numCollision]._localPointA = collision._localPointA;
          contactM._localMemory[contactM._numCollision]._localPointB = collision._localPointB;
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
    auto globalPointA = cA + nA2B * contactM._lhsRadius;
    auto globalPointB = cB - nA2B * contactM._rhsRadius;
    contactM._localMemory[contactM._numCollision]._localPointA =
      contactM._lhs->_q.conjugate().toRotationMatrix()*(globalPointA-contactM._lhs->_x);
    contactM._localMemory[contactM._numCollision]._localPointB =
      contactM._rhs->_q.conjugate().toRotationMatrix()*(globalPointB-contactM._rhs->_x);
    contactM._localMemory[contactM._numCollision]._globalNormal = nA2B;
    contactM._localMemory[contactM._numCollision]._isValid = true;
    contactM._numCollision++;
    return 1;
  }
  static DEVICE_HOST void generateManifoldSphereCapsuleInternal(
    GPUPBD::Collision<T>& originCollision, T& originCollisionDepth,
    const Vec3T &cA,const Capsule<T>* lhs,int lhs_idx,
    const Capsule<T>* rhs,int rhs_idx) {

    GPUPBD::Collision<T> collision;
    auto cB1 = rhs->globalMinCorner().template cast<T>();
    auto cB2 = rhs->globalMaxCorner().template cast<T>();
    auto cARadius = lhs->_radius;
    auto cBRadius = rhs->_radius;
    Vec3T n = cB2 - cB1;
    T nLenSqr = n.squaredNorm(), nLen = sqrt((double)nLenSqr);
    n /= nLen;
    T d = (cA - cB1).dot(n);
    T sumRad = cARadius + cBRadius, sumRadSqr = sumRad * sumRad;
    // three cases
    Vec3T globalPointA;
    Vec3T globalPointB;
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
      globalPointA = cA + collision._globalNormal * cARadius;
      globalPointB = cB1 - collision._globalNormal * cBRadius;
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
      globalPointA = cA + collision._globalNormal * cARadius;
      globalPointB = cB2 - collision._globalNormal * cBRadius;
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
      globalPointA = cA + collision._globalNormal * cARadius;
      globalPointB = cA - dir - collision._globalNormal * cBRadius;
    }
    collision._localPointA =
      lhs->_q.conjugate().toRotationMatrix()*(globalPointA-lhs->_x);
    collision._localPointB =
      rhs->_q.conjugate().toRotationMatrix()*(globalPointB-rhs->_x);
    auto collisionDepth = (globalPointA-globalPointB).dot(collision._globalNormal);
    if (!originCollision._isValid || originCollisionDepth < collisionDepth) {
      originCollisionDepth = collisionDepth;
      originCollision._capsuleIdA = lhs_idx;
      originCollision._capsuleIdB = rhs_idx;
      originCollision._localPointA = collision._localPointA;
      originCollision._localPointB = collision._localPointB;
      originCollision._globalNormal = collision._globalNormal;
      originCollision._isValid = true;
    }
  }
};
}

#endif
