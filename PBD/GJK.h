#ifndef GJK_H
#define GJK_H

#include "Pragma.h"
#include "Geometry.h"
namespace GPUPBD {

template <typename T>
DEVICE_HOST void sort2(T& a,T& b) {
  if(a>b) {
    T tmp = a;
    a = b;
    b = tmp;
  }
}
template <typename T>
DEVICE_HOST void sort3(T& a,T& b,T& c) {
  T tmp;
  if(a>b) {
    tmp = a;
    a = b;
    b = tmp;
  }
  if(a>c) {
    tmp = a;
    a = c;
    c = tmp;
  }
  if(b>c) {
    tmp = b;
    b = c;
    c = tmp;
  }
}
template <typename T>
struct Trans {
  DECL_MAT_VEC_MAP_TYPES_T
  Vec3T _x;
  QuatT _q;
};

template <typename T>
struct GJKPoint {
  DECL_MAT_VEC_MAP_TYPES_T
  DEVICE_HOST void calculate(const Trans<T>& transA,const Trans<T>& transB);
  Vec3T _ptAL,_ptBL,_ptAB;
  int _idA,_idB;
};
template <typename T>
struct GJK {
  DECL_MAT_VEC_MAP_TYPES_T
  static DEVICE_HOST Vec3T computeD(const GJKPoint<T> v[4],int nrP,T* bary,
                                    const Trans<T>& transA,
                                    const Trans<T>& transB,
                                    Vec3T& pAL,Vec3T& pBL);
  static DEVICE_HOST T runGJK(const Shape<T>* A,
                              const Shape<T>* B,
                              const Trans<T>& transA,
                              const Trans<T>& transB,
                              Vec3T& pAL,Vec3T& pBL,bool* intersect);
};
}
#endif
