#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <experimental/filesystem>
#include "ParallelVector.h"
#include "Pragma.h"

//openmp macro prepare
#ifdef _MSC_VER
#define STRINGIFY_OMP(X) X
#define PRAGMA __pragma
#else
#define STRINGIFY_OMP(X) #X
#define PRAGMA _Pragma
#endif
#ifndef NO_OPENMP
//openmp convenient functions
#define OMP_PARALLEL_FOR_ PRAGMA(STRINGIFY_OMP(omp parallel for schedule(static)))
#define OMP_PARALLEL_FOR_I(...) PRAGMA(STRINGIFY_OMP(omp parallel for schedule(static) __VA_ARGS__))
#define OMP_PARALLEL_FOR_X(X) PRAGMA(STRINGIFY_OMP(omp parallel for num_threads(X) schedule(static)))
#define OMP_PARALLEL_FOR_XI(X,...) PRAGMA(STRINGIFY_OMP(omp parallel for num_threads(X) schedule(static) __VA_ARGS__))
#define OMP_ADD(...) reduction(+: __VA_ARGS__)
#define OMP_PRI(...) private(__VA_ARGS__)
#define OMP_FPRI(...) firstprivate(__VA_ARGS__)
#define OMP_ATOMIC_ PRAGMA(STRINGIFY_OMP(omp atomic))
#ifdef _MSC_VER
#define OMP_ATOMIC_CAPTURE_ PRAGMA(STRINGIFY_OMP(omp critical))	// VS doesn't support capture, use critical instead
#else
#define OMP_ATOMIC_CAPTURE_ PRAGMA(STRINGIFY_OMP(omp atomic capture))
#endif
#define OMP_CRITICAL_ PRAGMA(STRINGIFY_OMP(omp critical))
#define OMP_FLUSH_(X) PRAGMA(STRINGIFY_OMP(omp flush(X)))
#else
//openmp convenient functions
#define OMP_PARALLEL_FOR_
#define OMP_PARALLEL_FOR_I(...)
#define OMP_PARALLEL_FOR_X(X)
#define OMP_PARALLEL_FOR_XI(X,...)
#define OMP_ADD(...)
#define OMP_PRI(...)
#define OMP_FPRI(...)
#define OMP_ATOMIC_
#define OMP_ATOMIC_CAPTURE_
#define OMP_CRITICAL_
#define OMP_FLUSH_(X)
#endif

namespace GPUPBD {
//filesystem
bool exists(const std::experimental::filesystem::v1::path& path);
void removeDir(const std::experimental::filesystem::v1::path& path);
void create(const std::experimental::filesystem::v1::path& path);
void recreate(const std::experimental::filesystem::v1::path& path);
std::vector<std::experimental::filesystem::v1::path> files(const std::experimental::filesystem::v1::path& path);
std::vector<std::experimental::filesystem::v1::path> directories(const std::experimental::filesystem::v1::path& path);
void sortFilesByNumber(std::vector<std::experimental::filesystem::v1::path>& files);
bool isDir(const std::experimental::filesystem::v1::path& path);
size_t fileSize(const std::experimental::filesystem::v1::path& path);
//basic functions
template <typename T>
inline void sort2(T& a,T& b) {
  if(a>b)
    std::swap(a,b);
}
template <typename T>
inline void sort3(T& a,T& b,T& c) {
  if(a>b)
    std::swap(a,b);
  if(b>c)
    std::swap(b,c);
  if(a>b)
    std::swap(a,b);
}
template <typename T>
inline void parallelAdd(T& a,T b) {
  OMP_CRITICAL_
  a+=b;
}
inline void parallelAdd(double& a,double b) {
  OMP_ATOMIC_
  a+=b;
}
template <typename T>
inline T maxCoeff(const Eigen::SparseMatrix<T,0,int>& coef) {
  T ret=0;
  for(int k=0; k<coef.outerSize(); ++k)
    for(typename Eigen::SparseMatrix<T,0,int>::InnerIterator it(coef,k); it; ++it)
      ret=std::max<T>(ret,abs(it.value()));
  return ret;
}
template <typename Derived>
inline typename Derived::Scalar maxCoeff(const Eigen::MatrixBase<Derived>& coef) {
  typename Derived::Scalar ret=0;
  for(int r=0; r<coef.rows(); r++)
    for(int c=0; c<coef.cols(); c++)
      ret=std::max(ret,coef(r,c));
  return ret;
}
template <typename T>
inline void parallelAddMat3TStride(Eigen::Matrix<T,-1,-1>& a,int offr,int offc,int strider,int stridec,const Eigen::Matrix<T,3,3>& b) {
  for(int d=0; d<3; d++)
    for(int d2=0; d2<3; d2++)
      parallelAdd<T>(a(offr+d*strider,offc+d2*stridec),b(d,d2));
}
template <typename T>
inline void parallelAddMat3TStride(ParallelVector<Eigen::Triplet<T,int>>& a,int offr,int offc,int strider,int stridec,const Eigen::Matrix<T,3,3>& b) {
  for(int d=0; d<3; d++)
    for(int d2=0; d2<3; d2++)
      a.push_back(Eigen::Triplet<T,int>(offr+d*strider,offc+d2*stridec,b(d,d2)));
}
}

#endif
