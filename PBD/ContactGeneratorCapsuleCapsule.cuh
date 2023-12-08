#ifndef CONTACTGENERATORCAPSULECAPSULE_CUH
#define CONTACTGENERATORCAPSULECAPSULE_CUH
#include "Geometry.h"
#include "Capsule.h"

// 在胶囊体碰撞中使用的球体球体碰撞
template <typename T>
__device__ __host__ inline void generateManifoldSphereSphereInternal(
  GPUPBD::Collision<T> *localMemory, size_t &offset, int &numCollision,
  const Eigen::Matrix<T, 3, 1> &cA1, const Eigen::Matrix<T, 3, 1> &cA2, T cARadius, int lhs_idx,
  const Eigen::Matrix<T, 3, 1> &cB1, const Eigen::Matrix<T, 3, 1> &cB2, T cBRadius, int rhs_idx)
{

  typedef Eigen::Matrix<T, 3, 1> Vec3T;
  typedef Eigen::Matrix<T, 4, 1> Vec4T;
  Vec4T distSqrs((cA1 - cB1).squaredNorm(), (cA2 - cB1).squaredNorm(), (cA1 - cB2).squaredNorm(), (cA2 - cB2).squaredNorm());
  typename Vec4T::Index id;
  T distSqr = distSqrs.minCoeff(&id), dist = 0;
  T sumRad = cARadius + cBRadius, sumRadSqr = sumRad * sumRad;
  const Vec3T &cA = (id % 2 == 0) ? cA1 : cA2;
  const Vec3T &cB = (id < 2) ? cB1 : cB2;
  // not in contact
  if (distSqr > sumRadSqr)
    return;
  // in contact
  Vec3T nA2B;
  if (distSqr > epsDist * epsDist)
  {
    // a single contact point
    nA2B = cB - cA;
    nA2B /= dist = sqrt((double)distSqr);
  }
  else
  {
    // overlapping degenerate case
    distSqrs.maxCoeff(&id);
    const Vec3T &cAn = (id % 2 == 0) ? cA1 : cA2;
    const Vec3T &cBn = (id < 2) ? cB1 : cB2;
    nA2B = cBn - cAn;
    nA2B /= nA2B.template cast<double>().norm();
  }
  localMemory[offset]._capsuleIdA = lhs_idx;
  localMemory[offset]._capsuleIdB = rhs_idx;
  localMemory[offset]._localPointA = cA + nA2B * cARadius;
  localMemory[offset]._localPointB = cB - nA2B * cBRadius;
  localMemory[offset]._globalNormal = nA2B;
  localMemory[offset]._isValid = true;
  offset++;
  numCollision++;
}

template <typename T>
__device__ __host__ inline void generateManifoldSphereCapsuleInternal(
  GPUPBD::Collision<T> &originCollision,
  const Eigen::Matrix<T, 3, 1> &cA, T cARadius, int lhs_idx,
  const Eigen::Matrix<T, 3, 1> &cB1, const Eigen::Matrix<T, 3, 1> &cB2, T cBRadius, int rhs_idx)
{

  typedef Eigen::Matrix<T, 3, 1> Vec3T;

  GPUPBD::Collision<T> collision;
  Vec3T n = cB2 - cB1;
  T nLenSqr = n.squaredNorm(), nLen = sqrt((double)nLenSqr);
  n /= nLen;
  T d = (cA - cB1).dot(n);
  T sumRad = cARadius + cBRadius, sumRadSqr = sumRad * sumRad;
  // three cases
  if (d <= 0)
  {
    T distSqr = (cA - cB1).squaredNorm(), dist = 0;
    // not in contact
    if (distSqr > sumRadSqr)
      return;
    // in contact
    if (distSqr > epsDist * epsDist)
    {
      // a single contact point
      collision._globalNormal = cB1 - cA;
      collision._globalNormal /= dist = sqrt((double)distSqr);
    }
    else
    {
      collision._globalNormal = n;
    }
    collision._localPointA = cA + collision._globalNormal * cARadius;
    collision._localPointB = cB1 - collision._globalNormal * cBRadius;
  }
  else if (d >= nLen)
  {
    T distSqr = (cA - cB2).squaredNorm(), dist = 0;
    // not in contact
    if (distSqr > sumRadSqr)
      return;
    // in contact
    if (distSqr > epsDist * epsDist)
    {
      // a single contact point
      collision._globalNormal = cB2 - cA;
      collision._globalNormal /= dist = sqrt((double)distSqr);
    }
    else
    {
      collision._globalNormal = -n;
    }
    collision._localPointA = cA + collision._globalNormal * cARadius;
    collision._localPointB = cB2 - collision._globalNormal * cBRadius;
  }
  else if (d > 0 && d < nLen)
  {
    Vec3T dir = cA - cB1 - n * d;
    T distSqr = dir.squaredNorm(), dist = 0;
    // not in contact
    if (distSqr > sumRadSqr)
      return;
    // in contact
    if (distSqr > epsDist * epsDist)
    {
      collision._globalNormal = -dir;
      collision._globalNormal /= dist = sqrt((double)distSqr);
    }
    else
    {
      typename Vec3T::Index id;
      n.cwiseAbs().minCoeff(&id);
      collision._globalNormal = n.cross(Vec3T::Unit(id));
      collision._globalNormal /= collision._globalNormal.template cast<double>().norm();
    }
    collision._localPointA = cA + collision._globalNormal * cARadius;
    collision._localPointB = cA - dir - collision._globalNormal * cBRadius;
  }
  if (!originCollision._isValid || originCollision.depth() < collision.depth())
  {
    originCollision._capsuleIdA = lhs_idx;
    originCollision._capsuleIdB = rhs_idx;
    originCollision._localPointA = collision._localPointA;
    originCollision._localPointB = collision._localPointB;
    originCollision._globalNormal = collision._globalNormal;
    originCollision._isValid = true;
  }
}

// 胶囊体碰撞
// localMemory中存入的是lhs单个胶囊体的全部碰撞。
// offset是本次碰撞检测在localMemory中的插入位置
// maxCollisionsPerNode是单个胶囊体碰撞的最大数目
template <typename T>
__device__ __host__ inline int narrowPhaseCollision(
  const GPUPBD::Capsule<T> &lhs, int lhs_idx,
  const GPUPBD::Capsule<T> &rhs, int rhs_idx,
  GPUPBD::Collision<T> *localMemory, size_t &offset, size_t maxCollisionsPerNode) noexcept
{
  typedef Eigen::Matrix<T, 2, 1> Vec2T;
  typedef Eigen::Matrix<T, 3, 1> Vec3T;
  typedef Eigen::Matrix<T, 3, 2> Mat3X2T;

  int numCollision = 0;
  Vec3T cA1 = ROT(lhs._trans) * lhs.minCorner().template cast<T>() + CTR(lhs._trans);
  Vec3T cA2 = ROT(lhs._trans) * lhs.maxCorner().template cast<T>() + CTR(lhs._trans);
  Vec3T cB1 = ROT(rhs._trans) * rhs.minCorner().template cast<T>() + CTR(rhs._trans);
  Vec3T cB2 = ROT(rhs._trans) * rhs.maxCorner().template cast<T>() + CTR(rhs._trans);
  Vec3T nA = cA2 - cA1, nB = cB2 - cB1;
  T nLenASqr = nA.squaredNorm(), nLenA = sqrt((double)nLenASqr);
  T nLenBSqr = nB.squaredNorm(), nLenB = sqrt((double)nLenBSqr);
  nA /= nLenA;
  nB /= nLenB;

  if (abs(nA.dot(nB)) > 1 - epsDir)
  {
    // nearly parallel 产生1~2个碰撞点
    T dB1 = (cB1 - cA1).dot(nA);
    T dB2 = (cB2 - cA1).dot(nA);
    if (dB1 <= 0 && dB2 <= 0)
    {
      // sphere-sphere 产生0/1个碰撞点
      generateManifoldSphereSphereInternal(localMemory, offset, numCollision,
                                           cA1, cA2, lhs._radius, lhs_idx,
                                           cB1, cB2, rhs._radius, rhs_idx);
    }
    else if (dB1 >= nLenA && dB2 >= nLenA)
    {
      // sphere-sphere 产生0/1个碰撞点
      generateManifoldSphereSphereInternal(localMemory, offset, numCollision,
                                           cA1, cA2, lhs._radius, lhs_idx,
                                           cB1, cB2, rhs._radius, rhs_idx);
    }
    else if (maxCollisionsPerNode - offset > 2)
    { // 假如localMemory不夠就忽略了
      // range 产生0/2个碰撞点
      Vec3T dir = cB1 - cA1 - dB1 * nA;
      T distSqr = dir.squaredNorm(), dist = 0;
      T sumRad = lhs._radius + rhs._radius, sumRadSqr = sumRad * sumRad;
      // not in contact
      if (distSqr > sumRadSqr)
        return numCollision;
      // in contact
      Vec3T nA2B;
      if (distSqr > epsDist * epsDist)
      {
        nA2B = dir;
        nA2B /= dist = sqrt((double)distSqr);
      }
      else
      {
        typename Vec3T::Index id;
        nA.cwiseAbs().minCoeff(&id);
        nA2B = nA.cross(Vec3T::Unit(id));
        nA2B /= nA2B.template cast<double>().norm();
      }
      // two contacts
      Vec2T range(std::max<T>(0, std::min(dB1, dB2)), std::min(nLenA, std::max(dB1, dB2)));
      for (const T &r : range)
      {
        localMemory[offset]._capsuleIdA = lhs_idx;
        localMemory[offset]._capsuleIdB = rhs_idx;
        localMemory[offset]._localPointA = cA1 + nA * r + nA2B * lhs._radius;
        localMemory[offset]._localPointB = cA1 + nA * r + dir - nA2B * rhs._radius;
        localMemory[offset]._globalNormal = nA2B;
        localMemory[offset]._isValid = true;
        offset++;
        numCollision++;
      }
    }
  }
  else
  {
    // not parallel
    Mat3X2T LHS;
    LHS.col(0) = -(cA2 - cA1);
    LHS.col(1) = (cB2 - cB1);
    Vec2T bary = (LHS.transpose() * LHS).inverse() * (LHS.transpose() * (cA1 - cB1));
    if ((bary.array() >= 0).all() && (bary.array() <= 1).all())
    {
      Vec3T cA = cA1 * (1 - bary[0]) + cA2 * bary[0];
      Vec3T cB = cB1 * (1 - bary[1]) + cB2 * bary[1];
      T distSqr = (cA - cB).squaredNorm();
      T sumRad = lhs._radius + rhs._radius, sumRadSqr = sumRad * sumRad;
      if (distSqr > sumRadSqr)
        return numCollision;
      Vec3T nA2B = nA.cross(nB);
      nA2B /= nA2B.template cast<double>().norm();
      if ((cB - cA).dot(nA2B) < 0)
        nA2B *= -1;
      localMemory[offset]._capsuleIdA = lhs_idx;
      localMemory[offset]._capsuleIdB = rhs_idx;
      localMemory[offset]._localPointA = cA + nA2B * lhs._radius;
      localMemory[offset]._localPointB = cB - nA2B * rhs._radius;
      localMemory[offset]._globalNormal = nA2B;
      localMemory[offset]._isValid = true;
      offset++;
      numCollision++;
    }
    else
    {
      GPUPBD::Collision<T> collision;
      generateManifoldSphereCapsuleInternal(collision,
                                            cA1, lhs._radius, lhs_idx,
                                            cB1, cB2, rhs._radius, rhs_idx);
      generateManifoldSphereCapsuleInternal(collision,
                                            cA2, lhs._radius, lhs_idx,
                                            cB1, cB2, rhs._radius, rhs_idx);
      generateManifoldSphereCapsuleInternal(collision,
                                            cB1, rhs._radius, rhs_idx,
                                            cA1, cA2, lhs._radius, lhs_idx);
      generateManifoldSphereCapsuleInternal(collision,
                                            cB2, rhs._radius, rhs_idx,
                                            cA1, cA2, lhs._radius, lhs_idx);
      if (collision._isValid)
      {
        localMemory[offset]._capsuleIdA = collision._capsuleIdA;
        localMemory[offset]._capsuleIdB = collision._capsuleIdB;
        localMemory[offset]._localPointA = collision._localPointA;
        localMemory[offset]._localPointB = collision._localPointB;
        localMemory[offset]._globalNormal = collision._globalNormal;
        localMemory[offset]._isValid = true;
        offset++;
        numCollision++;
      }
    }
  }
  return numCollision;
}

#endif