#ifndef NARROW_PHASE_H
#define NARROW_PHASE_H

#include "Geometry.h"
#include "Pragma.h"
#include "Collision.h"
#include "GJK.h"
#include "SAT.h"

namespace GPUPBD {
template<typename T>
class NarrowPhase {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  // maxCollisionsPerNode is the number of collisions per shape
  static DEVICE_HOST int narrowPhaseCollision(ContactManifold<T>& contactM,size_t maxCollisionsPerNode) noexcept {
    if(contactM._numCollision >= maxCollisionPerObject)
      return 0;
    if(ShapeType::Capsule==contactM._lhs->_type && ShapeType::Capsule==contactM._rhs->_type)
      return generateManifoldCapsuleCapsule(contactM, maxCollisionsPerNode);
    if(ShapeType::Capsule==contactM._lhs->_type && ShapeType::Box==contactM._rhs->_type)
      return generateManifoldCapsuleBox(contactM, maxCollisionsPerNode);
    if(ShapeType::Box==contactM._lhs->_type && ShapeType::Capsule==contactM._rhs->_type) {
      contactM.swap();
      int result =  generateManifoldCapsuleBox(contactM, maxCollisionsPerNode);
      contactM.swap();
      return result;
    }
    if(ShapeType::Box==contactM._lhs->_type && ShapeType::Box==contactM._rhs->_type)
      return generateManifoldBoxBox(contactM, maxCollisionsPerNode);
    return 0;
  }
 private:
  static DEVICE_HOST int generateManifoldCapsuleCapsule(ContactManifold<T>& contactM, size_t maxCollisionsPerNode) {
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
        T sumRad = contactM._lhs->_radius + contactM._rhs->_radius, sumRadSqr = sumRad * sumRad;
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
          T f = dB1 > dB2 ? 1.0 : -1.0;
          Vec3T globalPointA = cA1 + nA * (dB1 -r) * f + nA2B * contactM._lhs->_radius;
          Vec3T globalPointB = cB1 + nB * r - nA2B * contactM._rhs->_radius;
          T depth = (globalPointA-globalPointB).dot(nA2B);
          if(depth > 0) {
            contactM._localMemory[contactM._numCollision]._shapeIdA = contactM._lhsId;
            contactM._localMemory[contactM._numCollision]._shapeIdB = contactM._rhsId;
            contactM._localMemory[contactM._numCollision]._localPointA = contactM._lhs->_q.conjugate().toRotationMatrix()*(globalPointA-contactM._lhs->_x);
            contactM._localMemory[contactM._numCollision]._localPointB = contactM._rhs->_q.conjugate().toRotationMatrix()*(globalPointB-contactM._rhs->_x);
            contactM._localMemory[contactM._numCollision]._globalNormal = nA2B;
            contactM._localMemory[contactM._numCollision]._isValid = true;
            contactM._numCollision++;
            numCollision++;
          }
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
        T sumRad = contactM._lhs->_radius + contactM._rhs->_radius, sumRadSqr = sumRad * sumRad;
        if (distSqr > sumRadSqr)
          return numCollision;
        Vec3T nA2B = nA.cross(nB);
        nA2B /= nA2B.template cast<double>().norm();
        if ((cB - cA).dot(nA2B) < 0)
          nA2B *= -1;
        contactM._localMemory[contactM._numCollision]._shapeIdA = contactM._lhsId;
        contactM._localMemory[contactM._numCollision]._shapeIdB = contactM._rhsId;
        auto globalPointA = cA + nA2B * contactM._lhs->_radius;
        auto globalPointB = cB - nA2B * contactM._rhs->_radius;
        contactM._localMemory[contactM._numCollision]._localPointA = contactM._lhs->_q.conjugate().toRotationMatrix()*(globalPointA-contactM._lhs->_x);
        contactM._localMemory[contactM._numCollision]._localPointB = contactM._rhs->_q.conjugate().toRotationMatrix()*(globalPointB-contactM._rhs->_x);
        contactM._localMemory[contactM._numCollision]._globalNormal = nA2B;
        contactM._localMemory[contactM._numCollision]._isValid = true;
        contactM._numCollision++;
        numCollision++;
      } else {
        Constraint<T> collision;
        T collisionDepth = 0;
        generateManifoldSphereShapeInternal(collision, collisionDepth, cA1, contactM._lhs, contactM._lhsId, contactM._rhs, contactM._rhsId);
        generateManifoldSphereShapeInternal(collision,collisionDepth, cA2, contactM._lhs, contactM._lhsId, contactM._rhs, contactM._rhsId);
        generateManifoldSphereShapeInternal(collision, collisionDepth, cB1, contactM._rhs, contactM._rhsId, contactM._lhs, contactM._lhsId);
        generateManifoldSphereShapeInternal(collision, collisionDepth, cB2, contactM._rhs, contactM._rhsId, contactM._lhs, contactM._lhsId);
        if (collision._isValid) {
          contactM._localMemory[contactM._numCollision]._shapeIdA = collision._shapeIdA;
          contactM._localMemory[contactM._numCollision]._shapeIdB = collision._shapeIdB;
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
  // For internal use: sphere-sphere collision
  static DEVICE_HOST int generateManifoldSphereSphereInternal(ContactManifold<T>& contactM, size_t maxCollisionsPerNode) {
    Vec3T &cA1 = contactM._lhsMinCorner;
    Vec3T &cA2 = contactM._lhsMaxCorner;
    Vec3T &cB1 = contactM._rhsMinCorner;
    Vec3T &cB2 = contactM._rhsMaxCorner;
    Vec4T distSqrs((cA1 - cB1).squaredNorm(), (cA2 - cB1).squaredNorm(), (cA1 - cB2).squaredNorm(), (cA2 - cB2).squaredNorm());
    typename Vec4T::Index id;
    T distSqr = distSqrs.minCoeff(&id), dist = 0;
    T sumRad = contactM._lhs->_radius + contactM._rhs->_radius, sumRadSqr = sumRad * sumRad;
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
    contactM._localMemory[contactM._numCollision]._shapeIdA = contactM._lhsId;
    contactM._localMemory[contactM._numCollision]._shapeIdB = contactM._rhsId;
    auto globalPointA = cA + nA2B * contactM._lhs->_radius;
    auto globalPointB = cB - nA2B * contactM._rhs->_radius;
    contactM._localMemory[contactM._numCollision]._localPointA = contactM._lhs->_q.conjugate().toRotationMatrix()*(globalPointA-contactM._lhs->_x);
    contactM._localMemory[contactM._numCollision]._localPointB = contactM._rhs->_q.conjugate().toRotationMatrix()*(globalPointB-contactM._rhs->_x);
    contactM._localMemory[contactM._numCollision]._globalNormal = nA2B;
    contactM._localMemory[contactM._numCollision]._isValid = true;
    contactM._numCollision++;
    return 1;
  }
  static DEVICE_HOST void generateManifoldSphereShapeInternal(Constraint<T>& originCollision, T& originCollisionDepth, const Vec3T &cA,const Shape<T>* lhs,int lhs_idx, const Shape<T>* rhs,int rhs_idx) {
    Constraint<T> collision;
    auto cB1 = rhs->globalMinCorner().template cast<T>();
    auto cB2 = rhs->globalMaxCorner().template cast<T>();
    Vec3T n = cB2 - cB1;
    T nLenSqr = n.squaredNorm(), nLen = sqrt((double)nLenSqr);
    n /= nLen;
    T d = (cA - cB1).dot(n);
    T sumRad = lhs->_radius + rhs->_radius, sumRadSqr = sumRad * sumRad;
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
      globalPointA = cA + collision._globalNormal * lhs->_radius;
      globalPointB = cB1 - collision._globalNormal * rhs->_radius;
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
      globalPointA = cA + collision._globalNormal * lhs->_radius;
      globalPointB = cB2 - collision._globalNormal * rhs->_radius;
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
      globalPointA = cA + collision._globalNormal * lhs->_radius;
      globalPointB = cA - dir - collision._globalNormal * rhs->_radius;
    }
    collision._localPointA = lhs->_q.conjugate().toRotationMatrix()*(globalPointA-lhs->_x);
    collision._localPointB = rhs->_q.conjugate().toRotationMatrix()*(globalPointB-rhs->_x);
    auto collisionDepth = (globalPointA-globalPointB).dot(collision._globalNormal);
    if (!originCollision._isValid || originCollisionDepth < collisionDepth) {
      originCollisionDepth = collisionDepth;
      originCollision._shapeIdA = lhs_idx;
      originCollision._shapeIdB = rhs_idx;
      originCollision._localPointA = collision._localPointA;
      originCollision._localPointB = collision._localPointB;
      originCollision._globalNormal = collision._globalNormal;
      originCollision._isValid = true;
    }
  }
  static DEVICE_HOST int generateManifoldCapsuleBox(ContactManifold<T>& contactM, size_t maxCollisionsPerNode) {
    const Shape<T>* sA=contactM._lhs;
    const Shape<T>* sB=contactM._rhs;
    if(sA && sA->isCapsule() && sB) {
      // OMP_CRITICAL_ {
      //   //facets
      //   if(_facetCache.find(sA)==_facetCache.end())
      //     _facetCache[sA]=sA->facets();
      //   if(_facetCache.find(sB)==_facetCache.end())
      //     _facetCache[sB]=sB->facets();
      //   //edges
      //   if(_edgeCache.find(sA)==_edgeCache.end())
      //     _edgeCache[sA]=sA->edges();
      //   if(_edgeCache.find(sB)==_edgeCache.end())
      //     _edgeCache[sB]=sB->edges();
      // }
      // auto FA=_facetCache.find(sA)->second;
      FixedVector<Facet<T>,FACETSNUM> FB;
      sB->getFacets(FB);
      // auto EA=_edgeCache.find(sA)->second;
      // auto EB=_edgeCache.find(sB)->second;
      Vec3T pAL,pBL;
      bool intersect;
      const Trans<T> transA(sA->_x, sA->_q);
      const Trans<T> transB(sB->_x, sB->_q);
      T distSqr=GJK<T>::runGJK(sA,sB,transA,transB,pAL,pBL,&intersect);
      if(intersect || distSqr<epsDist*epsDist) {
        // SAT::generateManifold(sA,sB,FA,FB,EA,EB,m._tA,m._tB,m);
        // for(auto& p:m._points) {
        //   p._ptA+=p._nA2B*sA->radius();
        // }
      } else if(distSqr<sA->_radius*sA->_radius) {
        Facet<T> fA;
        Vec3T cA1=sA->globalMinCorner();
        Vec3T cA2=sA->globalMaxCorner();
        fA._boundary.push_back(ROT(transB).transpose()*(cA1-CTR(transB)));
        fA._boundary.push_back(ROT(transB).transpose()*(cA2-CTR(transB)));
        for(int i=0; i<FB.size(); i++) {
          const Facet<T>& fB = FB[i];
          if(abs(fB._n.dot(fA._boundary[0]-fA._boundary[1]))<epsDir && (fA._boundary[0]-fB._boundary[0]).dot(fB._n)>0) {
            //we can return multiple contacts
            SAT<T>::clip(fA,fB);
            int numFound = 0;
            for(int j=0; j<fA._boundary.size(); j++) {
              if(contactM._numCollision < maxCollisionsPerNode) {
                const Vec3T& pA = fA._boundary[j];
                Vec3T globalPointA=ROT(transB)*(pA-fB._n*sA->_radius)+CTR(transB);
                contactM._localMemory[contactM._numCollision]._shapeIdA = contactM._lhsId;
                contactM._localMemory[contactM._numCollision]._shapeIdB = contactM._rhsId;
                contactM._localMemory[contactM._numCollision]._localPointA = ROT(transA).transpose()*(globalPointA-CTR(transA));
                contactM._localMemory[contactM._numCollision]._localPointB = pA-fB._n*(pA-fB._boundary[0]).dot(fB._n);
                contactM._localMemory[contactM._numCollision]._globalNormal = -ROT(transB)*fB._n;
                contactM._localMemory[contactM._numCollision]._isValid = true;
                contactM._numCollision++;
                numFound++;
              } else break;
            }
            return numFound;
          }
        }
        //just return one closest point
        if(distSqr>epsDist*epsDist && contactM._numCollision < maxCollisionsPerNode) {
          Vec3T globalPointA = ROT(transA)*pAL+CTR(transA);
          Vec3T globalPointB = ROT(transB)*pBL+CTR(transB);
          Vec3T nA2B = (globalPointB-globalPointA)/sqrt((double)distSqr);
          globalPointA += nA2B*sA->_radius;
          contactM._localMemory[contactM._numCollision]._shapeIdA = contactM._lhsId;
          contactM._localMemory[contactM._numCollision]._shapeIdB = contactM._rhsId;
          contactM._localMemory[contactM._numCollision]._localPointA = ROT(transA).transpose()*(globalPointA-CTR(transA));
          contactM._localMemory[contactM._numCollision]._localPointB = pBL;
          contactM._localMemory[contactM._numCollision]._globalNormal = nA2B;
          contactM._localMemory[contactM._numCollision]._isValid = true;
          contactM._numCollision++;
          return 1;
        }
      }
    } else return 0;
    return 0;
  }
  static DEVICE_HOST int generateManifoldBoxBox(ContactManifold<T>& contactM, size_t maxCollisionsPerNode) {
    return 0;
  }
};
}

#endif
