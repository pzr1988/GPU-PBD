#ifndef PRAGMA_H
#define PRAGMA_H

//dll export
#if defined(_MSC_VER) and !defined(__CUDACC__)
#define _USE_MATH_DEFINES // for C++
#define DLL_EXPORT __declspec(dllexport)
#include <cmath>
#else
#define DLL_EXPORT
#endif

//cuda definitions
#ifdef __CUDACC__
#define DEVICE_HOST __device__ __host__
#define FINALIZE_GPU finalizeGPU();
#else
#define DEVICE_HOST
#define FINALIZE_GPU
#endif

//precision
#define LSCALAR float
#define epsDir 1e-3f
#define epsDist 1e-3f
#define MAXBOUNDARYSIZE 4
#define BOXFACENUM 6
#define BOXEDGENUM 12
#define CAPSULFACENUM 1
#define CAPSULEEDGENUM 1
const unsigned int maxCollisionPerObject=8;

//assert
#define ROT(A) A._q.toRotationMatrix()
#define CTR(A) A._x

#include <Eigen/Dense>

template <typename T>
class CopyableQuaternion : public Eigen::Quaternion<T> {
 public:
  DEVICE_HOST CopyableQuaternion():Eigen::Quaternion<T>(1,0,0,0) {}
  DEVICE_HOST CopyableQuaternion(const Eigen::AngleAxis<T>& aa):Eigen::Quaternion<T>(aa) {}
  DEVICE_HOST CopyableQuaternion(T w, T x, T y, T z):Eigen::Quaternion<T>(w,x,y,z) {}
  DEVICE_HOST CopyableQuaternion(const typename Eigen::Quaternion<T>::Coefficients& coeffs):Eigen::Quaternion<T>(coeffs) {}
  DEVICE_HOST CopyableQuaternion(const CopyableQuaternion<T>& other):Eigen::Quaternion<T>(other.w(),other.x(),other.y(),other.z()) {}
  DEVICE_HOST CopyableQuaternion(const Eigen::Quaternion<T>& other):Eigen::Quaternion<T>(other.w(),other.x(),other.y(),other.z()) {}
  DEVICE_HOST CopyableQuaternion<T>& operator=(const CopyableQuaternion<T>& other) {
    this->w() = other.w();
    this->x() = other.x();
    this->y() = other.y();
    this->z() = other.z();
    return *this;
  }
};

#define DECL_MAT_VEC_MAP_TYPES_T \
typedef Eigen::Matrix<T,2,1> Vec2T;\
typedef Eigen::Matrix<T,3,1> Vec3T;\
typedef Eigen::Matrix<T,4,1> Vec4T;\
typedef Eigen::Matrix<T,6,1> Vec6T;\
typedef Eigen::Matrix<T,3,3> Mat3T;\
typedef Eigen::Matrix<T,3,2> Mat3X2T;\
typedef CopyableQuaternion<T> QuatT;

#endif
