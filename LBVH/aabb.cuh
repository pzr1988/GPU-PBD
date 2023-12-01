#ifndef LBVH_AABB_CUH
#define LBVH_AABB_CUH
#include "utility.cuh"
#include <thrust/swap.h>
#include <cmath>
#include "PBD/Geometry.cuh"

namespace lbvh
{

template<typename T>
struct aabb
{
    typename vector_of<T>::type upper;
    typename vector_of<T>::type lower;
};

template<typename T>
__device__ __host__
inline bool intersects(const aabb<T>& lhs, const aabb<T>& rhs) noexcept
{
    if(lhs.upper.x < rhs.lower.x || rhs.upper.x < lhs.lower.x) {return false;}
    if(lhs.upper.y < rhs.lower.y || rhs.upper.y < lhs.lower.y) {return false;}
    if(lhs.upper.z < rhs.lower.z || rhs.upper.z < lhs.lower.z) {return false;}
    return true;
}

// 在胶囊体碰撞中使用的球体球体碰撞
template<typename T>
__device__ __host__
inline void generateManifoldSphereSphereInternal(
    ::GPUPBD::Collision<T>* localMemory, size_t& offset, int& numCollision,
    const Eigen::Matrix<T,3,1>& cA1,const Eigen::Matrix<T,3,1>& cA2, T cARadius, int lhs_idx, 
    const Eigen::Matrix<T,3,1>& cB1,const Eigen::Matrix<T,3,1>& cB2, T cBRadius, int rhs_idx,
    T epsDist 
    ) {

    typedef Eigen::Matrix<T,3,1> Vec3T;
    typedef Eigen::Matrix<T,4,1> Vec4T;
    Vec4T distSqrs((cA1-cB1).squaredNorm(),(cA2-cB1).squaredNorm(),(cA1-cB2).squaredNorm(),(cA2-cB2).squaredNorm());
    typename Vec4T::Index id;
    T distSqr=distSqrs.minCoeff(&id),dist=0;
    T sumRad=cARadius+cBRadius,sumRadSqr=sumRad*sumRad;
    const Vec3T& cA=(id%2==0)?cA1:cA2;
    const Vec3T& cB=(id<2)?cB1:cB2;
    //not in contact
    if(distSqr>sumRadSqr)
        return;
    //in contact
    Vec3T nA2B;
    if(distSqr>epsDist*epsDist) {
        //a single contact point
        nA2B=cB-cA;
        nA2B/=dist=sqrt((double)distSqr);
    } else {
        //overlapping degenerate case
        distSqrs.maxCoeff(&id);
        const Vec3T& cAn=(id%2==0)?cA1:cA2;
        const Vec3T& cBn=(id<2)?cB1:cB2;
        nA2B=cBn-cAn;
        nA2B/=nA2B.template cast<double>().norm();
    }
    localMemory[offset]._capsuleIdA = lhs_idx;
    localMemory[offset]._capsuleIdB = rhs_idx;
    localMemory[offset]._localPointA = cA+nA2B*cARadius;
    localMemory[offset]._localPointB = cB-nA2B*cBRadius;
    localMemory[offset]._globalNormal = nA2B;
    localMemory[offset]._isValid = true;
    offset++;
    numCollision++;
}
template<typename T>
__device__ __host__
inline void generateManifoldSphereCapsuleInternal(
    ::GPUPBD::Collision<T>& originCollision,
    const Eigen::Matrix<T,3,1>& cA, T cARadius, int lhs_idx, 
    const Eigen::Matrix<T,3,1>& cB1,const Eigen::Matrix<T,3,1>& cB2, T cBRadius, int rhs_idx,
    T epsDist ) {

    typedef Eigen::Matrix<T,3,1> Vec3T;

    ::GPUPBD::Collision<T> collision;
    Vec3T n=cB2-cB1;
    T nLenSqr=n.squaredNorm(),nLen=sqrt((double)nLenSqr);
    n/=nLen;
    T d=(cA-cB1).dot(n);
    T sumRad=cARadius+cBRadius,sumRadSqr=sumRad*sumRad;
    //three cases
    if(d<=0) {
        T distSqr=(cA-cB1).squaredNorm(),dist=0;
        //not in contact
        if(distSqr>sumRadSqr)
        return;
        //in contact
        if(distSqr>epsDist*epsDist) {
        //a single contact point
        collision._globalNormal=cB1-cA;
        collision._globalNormal/=dist=sqrt((double)distSqr);
        } else {
        collision._globalNormal=n;
        }
        collision._localPointA=cA+collision._globalNormal*cARadius;
        collision._localPointB=cB1-collision._globalNormal*cBRadius;
    } else if(d>=nLen)
    {
        T distSqr=(cA-cB2).squaredNorm(),dist=0;
        //not in contact
        if(distSqr>sumRadSqr)
            return;
        //in contact
        if(distSqr>epsDist*epsDist) {
            //a single contact point
            collision._globalNormal=cB2-cA;
            collision._globalNormal/=dist=sqrt((double)distSqr);
        } else {
            collision._globalNormal=-n;
        }
        collision._localPointA=cA+collision._globalNormal*cARadius;
        collision._localPointB=cB2-collision._globalNormal*cBRadius;
    } else if(d>0 && d<nLen) {
        Vec3T dir=cA-cB1-n*d;
        T distSqr=dir.squaredNorm(),dist=0;
            //not in contact
        if(distSqr>sumRadSqr)
            return;
        //in contact
        if(distSqr>epsDist*epsDist) {
            collision._globalNormal=-dir;
            collision._globalNormal/=dist=sqrt((double)distSqr);
        } else {
            typename Vec3T::Index id;
            n.cwiseAbs().minCoeff(&id);
            collision._globalNormal=n.cross(Vec3T::Unit(id));
            collision._globalNormal/=collision._globalNormal.template cast<double>().norm();
        }
        collision._localPointA=cA+collision._globalNormal*cARadius;
        collision._localPointB=cA-dir-collision._globalNormal*cBRadius;
    }
    if(!originCollision._isValid || originCollision.depth() < collision.depth()){
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
template<typename T>
__device__ __host__
inline int narrowPhaseCollision(
        const ::GPUPBD::Capsule<T>& lhs, int lhs_idx,
        const ::GPUPBD::Capsule<T>& rhs, int rhs_idx,
        ::GPUPBD::Collision<T>* localMemory, size_t &offset, size_t maxCollisionsPerNode,
        T epsDir, T epsDist) noexcept
{
    typedef Eigen::Matrix<T,2,1> Vec2T;
    typedef Eigen::Matrix<T,3,1> Vec3T;
    typedef Eigen::Matrix<T,3,2> Mat3X2T;

    int numCollision = 0;
    Vec3T cA1=ROT(lhs._trans)*lhs.minCorner().template cast<T>()+CTR(lhs._trans);
    Vec3T cA2=ROT(lhs._trans)*lhs.maxCorner().template cast<T>()+CTR(lhs._trans);
    Vec3T cB1=ROT(rhs._trans)*rhs.minCorner().template cast<T>()+CTR(lhs._trans);
    Vec3T cB2=ROT(rhs._trans)*rhs.maxCorner().template cast<T>()+CTR(lhs._trans);
    Vec3T nA=cA2-cA1,nB=cB2-cB1;
    T nLenASqr=nA.squaredNorm(),nLenA=sqrt((double)nLenASqr);
    T nLenBSqr=nB.squaredNorm(),nLenB=sqrt((double)nLenBSqr);
    nA/=nLenA;
    nB/=nLenB;
    
    if(abs(nA.dot(nB))>1-epsDir)
    {
      //nearly parallel 产生1~2个碰撞点
      T dB1=(cB1-cA1).dot(nA);
      T dB2=(cB2-cA1).dot(nA);
      if(dB1<=0 && dB2<=0) {
        //sphere-sphere 产生0/1个碰撞点
        generateManifoldSphereSphereInternal(localMemory, offset, numCollision,
                                             cA1, cA2, lhs._radius, lhs_idx,
                                             cB1, cB2, rhs._radius, rhs_idx,
                                             epsDist);
      } else if(dB1>=nLenA && dB2>=nLenA)
      {   
        //sphere-sphere 产生0/1个碰撞点
        generateManifoldSphereSphereInternal(localMemory, offset, numCollision,
                                             cA1, cA2, lhs._radius, lhs_idx,
                                             cB1, cB2, rhs._radius, rhs_idx,
                                             epsDist);
      } else if (maxCollisionsPerNode-offset>2) // 假如localMemory不夠就忽略了
      {
        //range 产生0/2个碰撞点
        Vec3T dir=cB1-cA1-dB1*nA;
        T distSqr=dir.squaredNorm(),dist=0;
        T sumRad=lhs._radius+rhs._radius,sumRadSqr=sumRad*sumRad;
        //not in contact
        if(distSqr>sumRadSqr)
          return numCollision;
        //in contact
        Vec3T nA2B;
        if(distSqr>epsDist*epsDist) {
          nA2B=dir;
          nA2B/=dist=sqrt((double)distSqr);
        } else {
          typename Vec3T::Index id;
          nA.cwiseAbs().minCoeff(&id);
          nA2B=nA.cross(Vec3T::Unit(id));
          nA2B/=nA2B.template cast<double>().norm();
        }
        //two contacts
        Vec2T range(std::max<T>(0,std::min(dB1,dB2)),std::min(nLenA,std::max(dB1,dB2)));
        for(const T& r:range)
        {
            localMemory[offset]._capsuleIdA = lhs_idx;
            localMemory[offset]._capsuleIdB = rhs_idx;
            localMemory[offset]._localPointA = cA1+nA*r+nA2B*lhs._radius;
            localMemory[offset]._localPointB = cA1+nA*r+dir-nA2B*rhs._radius;
            localMemory[offset]._globalNormal = nA2B;
            localMemory[offset]._isValid = true;
            offset++;
            numCollision++;
        }
      }
    } else {
      //not parallel
      Mat3X2T LHS;
      LHS.col(0)=-(cA2-cA1);
      LHS.col(1)= (cB2-cB1);
      Vec2T bary=(LHS.transpose()*LHS).inverse()*(LHS.transpose()*(cA1-cB1));
      if((bary.array()>=0).all() && (bary.array()<=1).all()) {
        Vec3T cA=cA1*(1-bary[0])+cA2*bary[0];
        Vec3T cB=cB1*(1-bary[1])+cB2*bary[1];
        T distSqr=(cA-cB).squaredNorm();
        T sumRad=lhs._radius+rhs._radius,sumRadSqr=sumRad*sumRad;
        if(distSqr>sumRadSqr)
          return true;
        Vec3T nA2B=nA.cross(nB);
        nA2B/=nA2B.template cast<double>().norm();
        if((cB-cA).dot(nA2B)<0)
          nA2B*=-1;
        localMemory[offset]._capsuleIdA = lhs_idx;
        localMemory[offset]._capsuleIdB = rhs_idx;
        localMemory[offset]._localPointA = cA+nA2B*lhs._radius;
        localMemory[offset]._localPointB = cB-nA2B*rhs._radius;
        localMemory[offset]._globalNormal = nA2B;
        localMemory[offset]._isValid = true;
        offset++;
        numCollision++;
      } else {
        ::GPUPBD::Collision<T> collision;
        generateManifoldSphereCapsuleInternal(collision,
                                              cA1, lhs._radius, lhs_idx,
                                              cB1, cB2, rhs._radius, rhs_idx,
                                              epsDist);
        generateManifoldSphereCapsuleInternal(collision,
                                              cA2, lhs._radius, lhs_idx,
                                              cB1, cB2, rhs._radius, rhs_idx,
                                              epsDist);
        generateManifoldSphereCapsuleInternal(collision,
                                              cB1, rhs._radius, rhs_idx,
                                              cA1, cA2, lhs._radius, lhs_idx,
                                              epsDist);
        generateManifoldSphereCapsuleInternal(collision,
                                              cB2, rhs._radius, rhs_idx,
                                              cA1, cA2, lhs._radius, lhs_idx,
                                              epsDist);
        if(collision._isValid){
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

__device__ __host__
inline aabb<double> merge(const aabb<double>& lhs, const aabb<double>& rhs) noexcept
{
    aabb<double> merged;
    merged.upper.x = ::fmax(lhs.upper.x, rhs.upper.x);
    merged.upper.y = ::fmax(lhs.upper.y, rhs.upper.y);
    merged.upper.z = ::fmax(lhs.upper.z, rhs.upper.z);
    merged.lower.x = ::fmin(lhs.lower.x, rhs.lower.x);
    merged.lower.y = ::fmin(lhs.lower.y, rhs.lower.y);
    merged.lower.z = ::fmin(lhs.lower.z, rhs.lower.z);
    return merged;
}

__device__ __host__
inline aabb<float> merge(const aabb<float>& lhs, const aabb<float>& rhs) noexcept
{
    aabb<float> merged;
    merged.upper.x = ::fmaxf(lhs.upper.x, rhs.upper.x);
    merged.upper.y = ::fmaxf(lhs.upper.y, rhs.upper.y);
    merged.upper.z = ::fmaxf(lhs.upper.z, rhs.upper.z);
    merged.lower.x = ::fminf(lhs.lower.x, rhs.lower.x);
    merged.lower.y = ::fminf(lhs.lower.y, rhs.lower.y);
    merged.lower.z = ::fminf(lhs.lower.z, rhs.lower.z);
    return merged;
}

// metrics defined in
// Nearest Neighbor Queries (1995) ACS-SIGMOD
// - Nick Roussopoulos, Stephen Kelley FredericVincent

__device__ __host__
inline float mindist(const aabb<float>& lhs, const float4& rhs) noexcept
{
    const float dx = ::fminf(lhs.upper.x, ::fmaxf(lhs.lower.x, rhs.x)) - rhs.x;
    const float dy = ::fminf(lhs.upper.y, ::fmaxf(lhs.lower.y, rhs.y)) - rhs.y;
    const float dz = ::fminf(lhs.upper.z, ::fmaxf(lhs.lower.z, rhs.z)) - rhs.z;
    return dx * dx + dy * dy + dz * dz;
}

__device__ __host__
inline double mindist(const aabb<double>& lhs, const double4& rhs) noexcept
{
    const double dx = ::fmin(lhs.upper.x, ::fmax(lhs.lower.x, rhs.x)) - rhs.x;
    const double dy = ::fmin(lhs.upper.y, ::fmax(lhs.lower.y, rhs.y)) - rhs.y;
    const double dz = ::fmin(lhs.upper.z, ::fmax(lhs.lower.z, rhs.z)) - rhs.z;
    return dx * dx + dy * dy + dz * dz;
}

__device__ __host__
inline float minmaxdist(const aabb<float>& lhs, const float4& rhs) noexcept
{
    float3 rm_sq = make_float3((lhs.lower.x - rhs.x) * (lhs.lower.x - rhs.x),
                               (lhs.lower.y - rhs.y) * (lhs.lower.y - rhs.y),
                               (lhs.lower.z - rhs.z) * (lhs.lower.z - rhs.z));
    float3 rM_sq = make_float3((lhs.upper.x - rhs.x) * (lhs.upper.x - rhs.x),
                               (lhs.upper.y - rhs.y) * (lhs.upper.y - rhs.y),
                               (lhs.upper.z - rhs.z) * (lhs.upper.z - rhs.z));
    
    if((lhs.upper.x + lhs.lower.x) * 0.5f < rhs.x)
    {
        thrust::swap(rm_sq.x, rM_sq.x);
    }
    if((lhs.upper.y + lhs.lower.y) * 0.5f < rhs.y)
    {
        thrust::swap(rm_sq.y, rM_sq.y);
    }
    if((lhs.upper.z + lhs.lower.z) * 0.5f < rhs.z)
    {
        thrust::swap(rm_sq.z, rM_sq.z);
    }
    
    const float dx = rm_sq.x + rM_sq.y + rM_sq.z;
    const float dy = rM_sq.x + rm_sq.y + rM_sq.z;
    const float dz = rM_sq.x + rM_sq.y + rm_sq.z;
    return ::fminf(dx, ::fminf(dy, dz));
}

__device__ __host__
inline double minmaxdist(const aabb<double>& lhs, const double4& rhs) noexcept
{
    double3 rm_sq = make_double3((lhs.lower.x - rhs.x) * (lhs.lower.x - rhs.x),
                                 (lhs.lower.y - rhs.y) * (lhs.lower.y - rhs.y),
                                 (lhs.lower.z - rhs.z) * (lhs.lower.z - rhs.z));
    double3 rM_sq = make_double3((lhs.upper.x - rhs.x) * (lhs.upper.x - rhs.x),
                                 (lhs.upper.y - rhs.y) * (lhs.upper.y - rhs.y),
                                 (lhs.upper.z - rhs.z) * (lhs.upper.z - rhs.z));

    if((lhs.upper.x + lhs.lower.x) * 0.5 < rhs.x)
    {
        thrust::swap(rm_sq.x, rM_sq.x);
    }
    if((lhs.upper.y + lhs.lower.y) * 0.5 < rhs.y)
    {
        thrust::swap(rm_sq.y, rM_sq.y);
    }
    if((lhs.upper.z + lhs.lower.z) * 0.5 < rhs.z)
    {
        thrust::swap(rm_sq.z, rM_sq.z);
    }

    const double dx = rm_sq.x + rM_sq.y + rM_sq.z;
    const double dy = rM_sq.x + rm_sq.y + rM_sq.z;
    const double dz = rM_sq.x + rM_sq.y + rm_sq.z;
    return ::fmin(dx, ::fmin(dy, dz));
}

template<typename T>
__device__ __host__
inline typename vector_of<T>::type centroid(const aabb<T>& box) noexcept
{
    typename vector_of<T>::type c;
    c.x = (box.upper.x + box.lower.x) * 0.5;
    c.y = (box.upper.y + box.lower.y) * 0.5;
    c.z = (box.upper.z + box.lower.z) * 0.5;
    return c;
}

// 获得物体的bounding box
template<template<typename> class Geometry, typename T>
struct aabb_getter
{
    __device__
    lbvh::aabb<float> operator()(const Geometry<T>& c) const noexcept
    {
        lbvh::aabb<float> retval;
        return retval;
    }
};

// 获得胶囊体的bounding box
template<typename T>
struct aabb_getter<::GPUPBD::Capsule, T>
{
    __device__
    lbvh::aabb<float> operator()(const GPUPBD::Capsule<T> &c) const noexcept
    {
        lbvh::aabb<float> retval;
        Eigen::Matrix<float, 4, 1> end1(static_cast<float>(c._len)/2.0, 0, 0, 1); // 第一个端点
        Eigen::Matrix<float, 4, 1> end2(-static_cast<float>(c._len)/2.0, 0, 0, 1); // 第二个端点

        Eigen::Matrix<float, 3, 1> transformedEnd1 = c._trans * end1;
        Eigen::Matrix<float, 3, 1> transformedEnd2 = c._trans * end2;
        
        Eigen::Matrix<float, 3, 1> upper = transformedEnd1.head<3>().cwiseMax(transformedEnd2.head<3>());
        float radius = static_cast<float>(c._radius);
        retval.upper.x = upper.x() + radius;
        retval.upper.y = upper.y() + radius;
        retval.upper.z = upper.z() + radius;
        Eigen::Matrix<float, 3, 1> lower = transformedEnd1.head<3>().cwiseMin(transformedEnd2.head<3>()) ;
        retval.lower.x = lower.x() - radius;
        retval.lower.y = lower.y() - radius;
        retval.lower.z = lower.z() - radius;
        return retval;
    }
};

} // lbvh
#endif// LBVH_AABB_CUH
