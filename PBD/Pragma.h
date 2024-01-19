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

#include <Eigen/Dense>

#define DECL_MAT_VEC_MAP_TYPES_T \
typedef Eigen::Matrix<T,2,1> Vec2T;\
typedef Eigen::Matrix<T,3,1> Vec3T;\
typedef Eigen::Matrix<T,4,1> Vec4T;\
typedef Eigen::Matrix<T,3,3> Mat3T;\
typedef Eigen::Matrix<T,3,2> Mat3X2T;\
typedef Eigen::Quaternion<T> QuatT;

#endif
