#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <experimental/filesystem>
#include "ParallelVector.h"
#include "Pragma.h"

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
