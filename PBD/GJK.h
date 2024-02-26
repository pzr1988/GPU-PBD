#ifndef GJK_H
#define GJK_H

#include "Pragma.h"
#include "Geometry.h"
#include "DistanceFunction.h"
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
template <int DIM, typename T>
DEVICE_HOST void updateFeat(GJKPoint<T> v[4],int& nrP,int* feat) {
  if(feat[0]==-1)
    return;
  else if(DIM>=3 && feat[2]>=0) {
    sort3(feat[0],feat[1],feat[2]);
    v[0]=v[feat[0]];
    v[1]=v[feat[1]];
    v[2]=v[feat[2]];
    nrP=3;
    return;
  } else if(DIM>=2 && feat[1]>=0) {
    sort2(feat[0],feat[1]);
    v[0]=v[feat[0]];
    v[1]=v[feat[1]];
    nrP=2;
    return;
  } else if(DIM>=1 && feat[0]>=0) {
    v[0]=v[feat[0]];
    nrP=1;
    return;
  }
}
template <typename T>
DEVICE_HOST void GJKPoint<T>::calculate(const Trans<T>& transA,const Trans<T>& transB) {
  _ptAB = transA._q.toRotationMatrix()*_ptAL+transA._x;
  _ptAB-= transB._q.toRotationMatrix()*_ptBL+transB._x;
}
template <typename T>
DEVICE_HOST typename GJK<T>::Vec3T GJK<T>::computeD(const GJKPoint<T> v[4],int nrP,T* bary,
    const Trans<T>& transA,
    const Trans<T>& transB,
    Vec3T& pAL,Vec3T& pBL) {
  pAL=pBL=Vec3T::Zero();
  for(int d=0; d<nrP; d++) {
    pAL+=v[d]._ptAL*bary[d];
    pBL+=v[d]._ptBL*bary[d];
  }
  return (transA._q.toRotationMatrix()*pAL+transA._x)-(transB._q.toRotationMatrix()*pBL+transB._x);
}
template <typename T>
DEVICE_HOST T GJK<T>::runGJK(const Shape<T>* A,
                             const Shape<T>* B,
                             const Trans<T>& transA,
                             const Trans<T>& transB,
                             Vec3T& pAL,Vec3T& pBL,bool* intersect) {
  int nrP;
  Vec3T cp,D;
  GJKPoint<T> v[4];
  T dist,minDist;
  Eigen::Matrix<T,4,1> bary;
  //initialize
  nrP=1;
  v[0]._ptAL=A->support(Vec3T::Zero(),v[0]._idA);
  v[0]._ptBL=B->support(Vec3T::Zero(),v[0]._idB);
  v[0].calculate(transA,transB);
  dist=minDist=FLT_MAX;
  D=v[0]._ptAB;
  //main loop
  if(intersect)
    *intersect=false;
  while(true) {
    //insert new point
    v[nrP]._ptAL=A->support(-transA._q.inverse().toRotationMatrix()*D,v[nrP]._idA);
    v[nrP]._ptBL=B->support( transB._q.inverse().toRotationMatrix()*D,v[nrP]._idB);
    v[nrP++].calculate(transA,transB);
    //calculate distance
    if(nrP==4) {
      Eigen::Matrix<int,3,1> feat;
      Vec3T vAB[4]= {v[0]._ptAB,v[1]._ptAB,v[2]._ptAB,v[3]._ptAB};
      dist=distToSqrTetrahedron<T>(Vec3T::Zero(),vAB,Eigen::Map<Eigen::Matrix<T,4,1>>(bary.data()),cp,&feat);
      D=computeD(v,nrP,bary.data(),transA,transB,pAL,pBL);
      updateFeat<3,T>(v,nrP,feat.data());
      if(dist>=minDist)
        break;
      else {
        minDist=dist;
        if(feat[0]==-1) {
          if(intersect)
            *intersect=true;
          break;
        }
      }
    } else if(nrP==3) {
      Eigen::Matrix<int,2,1> feat;
      Vec3T vAB[3]= {v[0]._ptAB,v[1]._ptAB,v[2]._ptAB};
      dist=distToSqrTriangle<T>(Vec3T::Zero(),vAB,Eigen::Map<Eigen::Matrix<T,3,1>>(bary.data()),cp,&feat);
      D=computeD(v,nrP,bary.data(),transA,transB,pAL,pBL);
      updateFeat<2,T>(v,nrP,feat.data());
      if(dist>=minDist)
        break;
      else minDist=dist;
    } else if(nrP==2) {
      int feat;
      Vec3T vAB[2]= {v[0]._ptAB,v[1]._ptAB};
      dist=distToSqrLineSegment<T>(Vec3T::Zero(),vAB,Eigen::Map<Eigen::Matrix<T,2,1>>(bary.data()),cp,&feat);
      D=computeD(v,nrP,bary.data(),transA,transB,pAL,pBL);
      updateFeat<1,T>(v,nrP,&feat);
      if(dist>=minDist)
        break;
      else minDist=dist;
    }
  }
  return minDist;
}
}
#endif
