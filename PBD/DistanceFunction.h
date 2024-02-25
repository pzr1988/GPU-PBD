#ifndef DISTANCE_H
#define DISTANCE_H

#include "Pragma.h"

namespace GPUPBD {
template <typename T>
DEVICE_HOST T distToSqrLineSegment(const Eigen::Matrix<T,3,1>& pt,
                                   const Eigen::Matrix<T,3,1> v[2],
                                   Eigen::Map<Eigen::Matrix<T,2,1>> bary,
                                   Eigen::Matrix<T,3,1>& cp,
                                   int* feat) {
  DECL_MAT_VEC_MAP_TYPES_T
  Vec3T LHSE=v[1]-v[0],RHSE=pt-v[0];
  bool systemInvertible=true;
  T alpha=RHSE.dot(LHSE)/LHSE.dot(LHSE);
  if(!isfinite(alpha))
    systemInvertible=false;
  if(!systemInvertible) {
    alpha=0;
    if(feat)
      *feat=0;
  } else if(alpha>=0 && alpha<=1) {
    if(feat)
      *feat=-1;
  } else if(alpha<0) {
    alpha=0;
    if(feat)
      *feat=0;
  } else {
    alpha=1;
    if(feat)
      *feat=1;
  }
  bary[1]=alpha;
  bary[0]=1-alpha;
  cp=bary[0]*v[0]+bary[1]*v[1];
  return (pt-cp).squaredNorm();
}
template <typename T>
DEVICE_HOST T distToSqrTriangle(const Eigen::Matrix<T,3,1>& pt,
                                const Eigen::Matrix<T,3,1> v[3],
                                Eigen::Map<Eigen::Matrix<T,3,1>> bary,
                                Eigen::Matrix<T,3,1>& cp,
                                Eigen::Matrix<int,2,1>* feat) {
  DECL_MAT_VEC_MAP_TYPES_T
  Mat3X2T LHS;
  LHS.col(0)=v[1]-v[0];
  LHS.col(1)=v[2]-v[0];
  Vec3T RHS=pt-v[0];
  T alpha;
  //bary
  bool systemInvertible=true;
  bary.template segment<2>(1)=(LHS.transpose()*LHS).inverse()*(LHS.transpose()*RHS);
  bary[0]=1-bary.template segment<2>(1).sum();

  if(!bary.template cast<double>().array().isFinite().all())
    systemInvertible=false;

  if(!systemInvertible || bary.minCoeff()<0) {
    T dist,minDist=100000.f;
    //edge
    bool needTestV[3]= {true,true,true};
    for(int d=0; d<3; d++) {
      //|v[(d+1)%3+1]*alpha+v[d+1]*(1-alpha)-v[0]|^2
      Vec3T LHSE=v[(d+1)%3]-v[d],RHSE=pt-v[d];
      systemInvertible=true;
      alpha=RHSE.dot(LHSE)/LHSE.dot(LHSE);
      if(!isfinite(alpha))
        systemInvertible=false;
      if(systemInvertible && alpha>=0 && alpha<=1) {
        needTestV[(d+1)%3]=needTestV[d]=false;
        dist=(LHSE*alpha-RHSE).squaredNorm();
        if(dist<minDist) {
          if(feat)
            *feat=Eigen::Matrix<int,2,1>((d+1)%3,d);
          bary.setZero();
          bary[(d+1)%3]=alpha;
          bary[d]=1-alpha;
          minDist=dist;
        }
      }
    }
    //vertex
    for(int d=0; d<3; d++)
      if(needTestV[d]) {
        dist=(v[d]-pt).squaredNorm();
        if(dist<minDist) {
          if(feat)
            *feat=Eigen::Matrix<int,2,1>(d,-1);
          bary.setUnit(d);
          minDist=dist;
        }
      }
  } else if(feat)
    *feat=Eigen::Matrix<int,2,1>(-1,-1);
  cp=bary[0]*v[0]+bary[1]*v[1]+bary[2]*v[2];
  return (pt-cp).squaredNorm();
}

template <typename T>
DEVICE_HOST T distToSqrTetrahedron(const Eigen::Matrix<T,3,1>& pt,
                                   const Eigen::Matrix<T,3,1> v[4],
                                   Eigen::Map<Eigen::Matrix<T,4,1>> bary,
                                   Eigen::Matrix<T,3,1>& cp,
                                   Eigen::Matrix<int,3,1>* feat) {
  DECL_MAT_VEC_MAP_TYPES_T
  using namespace std;
  Mat3T LHS;
  LHS.col(0)=v[1]-v[0];
  LHS.col(1)=v[2]-v[0];
  LHS.col(2)=v[3]-v[0];
  Vec3T RHS=pt-v[0];
  bool inside=true;
  bary.template segment<3>(1)=LHS.inverse()*RHS;
  bary[0]=1-bary.template segment<3>(1).sum();
  if(!bary.template cast<double>().array().isFinite().all())
    inside=false;
  if(!inside || bary.minCoeff()<0) {
    T dist,minDist=std::numeric_limits<double>::max();
    //consider 4 triangles
    Eigen::Matrix<T,3,1> baryT;
    Eigen::Matrix<T,3,1> vt[3],cpT;
    Eigen::Matrix<int,2,1> featT;
    for(int d=0; d<4; d++) {
      vt[0]=v[(d+1)%4];
      vt[1]=v[(d+2)%4];
      vt[2]=v[(d+3)%4];
      dist=distToSqrTriangle(pt,vt,Eigen::Map<Vec3T>(baryT.data()),cpT,feat?&featT:NULL);
      if(dist<minDist) {
        minDist=dist;
        bary.setConstant(0);
        bary[(d+1)%4]=baryT[0];
        bary[(d+2)%4]=baryT[1];
        bary[(d+3)%4]=baryT[2];
        cp=cpT;
        if(feat) {
          if(featT==Eigen::Matrix<int,2,1>(-1,-1))
            *feat=Eigen::Matrix<int,3,1>((d+1)%4,(d+2)%4,(d+3)%4);
          else if(featT[1]==-1)
            *feat=Eigen::Matrix<int,3,1>((d+1+featT[0])%4,-1,-1);
          else
            *feat=Eigen::Matrix<int,3,1>((d+1+featT[0])%4,(d+1+featT[1])%4,-1);
        }
      }
    }
  } else if(feat)
    *feat=Eigen::Matrix<int,3,1>(-1,-1,-1);
  cp=bary[0]*v[0]+bary[1]*v[1]+bary[2]*v[2]+bary[3]*v[3];
  return (pt-cp).squaredNorm();
}
}
#endif