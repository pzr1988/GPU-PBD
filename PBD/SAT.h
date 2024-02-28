#ifndef SAT_H
#define SAT_H

#include "Pragma.h"
#include "Geometry.h"
#include "DistanceFunction.h"

namespace GPUPBD {
template <typename T,typename TF>
DEVICE_HOST T interp1D(const T& v0,const T& v1,
                       const TF& px) {
  return v0*(TF)(1.0f-px)+v1*(TF)px;
}
template <typename T>
struct SAT {
  DECL_MAT_VEC_MAP_TYPES_T
  struct ProjRange {
    int _fidA,_eidA;
    int _fidB,_eidB;
    Vec2T _rngA,_rngB;
    Vec3T _n;
    T _depth;
    //additional data for edge-edge check
    Vec3T _ptA,_ptB;
  };
  static DEVICE_HOST T depth(const ProjRange& rng);
  static DEVICE_HOST bool matchWitness(const ProjRange& rng,T witnessA,T* verboseB);
  static DEVICE_HOST bool matchWitness(const ProjRange& rng,T* verboseA,T witnessB);
  static DEVICE_HOST bool matchWitness(const ProjRange& rng,T witnessA,T witnessB);
  static Vec2T project(Shape<T>* S,const Trans<T>& trans,const Vec3T& n);
  static DEVICE_HOST void clip(Facet<T>& f,const Vec3T& pos,const Vec3T& inward);
  static DEVICE_HOST void clip(Facet<T>& f,const Facet<T>& ref);
};
template <typename T>
T SAT<T>::depth(const ProjRange& rng) {
  T a = rng._rngA[1]-rng._rngB[0];
  T b = rng._rngB[1]-rng._rngA[0];
  if(a<b) return a;
  return b;
}
template <typename T>
bool SAT<T>::matchWitness(const ProjRange& rng,T witnessA,T*) {
  if(rng._rngA[1]-rng._rngB[0]<rng._rngB[1]-rng._rngA[0]) {
    //check matching witness
    if(abs(witnessA-rng._rngA[1])<epsDist)
      return true;
  } else {
    //check matching witness
    if(abs(witnessA-rng._rngA[0])<epsDist)
      return true;
  }
  //witness does not match
  return false;
}
template <typename T>
bool SAT<T>::matchWitness(const ProjRange& rng,T*,T witnessB) {
  if(rng._rngA[1]-rng._rngB[0]<rng._rngB[1]-rng._rngA[0]) {
    //check matching witness
    if(abs(witnessB-rng._rngB[0])<epsDist)
      return true;
  } else {
    //check matching witness
    if(abs(witnessB-rng._rngB[1])<epsDist)
      return true;
  }
  //witness does not match
  return false;
}
template <typename T>
bool SAT<T>::matchWitness(const ProjRange& rng,T witnessA,T witnessB) {
  if(rng._rngA[1]-rng._rngB[0]<rng._rngB[1]-rng._rngA[0]) {
    //check matching witness
    if( abs(witnessA-rng._rngA[1])<epsDist &&
        abs(witnessB-rng._rngB[0])<epsDist)
      return true;
  } else {
    //check matching witness
    if( abs(witnessA-rng._rngA[0])<epsDist &&
        abs(witnessB-rng._rngB[1])<epsDist)
      return true;
  }
  //witness does not match
  return false;
}
template <typename T>
typename SAT<T>::Vec2T SAT<T>::project(Shape<T>* S,const Trans<T>& trans,const Vec3T& n) {
  Vec2T ret=S->project(ROT(trans).transpose()*n);
  ret=(ret.array()+n.dot(CTR(trans))).matrix();
  sort2(ret[0],ret[1]);
  return ret;
}
template <typename T>
void SAT<T>::clip(Facet<T>& f,const Vec3T& pos,const Vec3T& inward) {
  int nr=(int)f._boundary.size();
  bool hasIn=false,hasOut=false;
  FixedVector<T, MAXBOUNDARYSIZE> inside;
  for(int i=0; i<nr; i++) {
    inside[i]=(f._boundary[i]-pos).dot(inward);
    hasIn |=inside[i]>0;
    hasOut|=inside[i]<=0;
  }
  if(!hasIn) {
    f._boundary.clear();
    return;
  }
  if(!hasOut)
    return;
  //Sutherland–Hodgman algorithm
  FixedVector<Vec3T, MAXBOUNDARYSIZE> boundaryNew;
  for(int curr=0, next=1; curr<nr; curr++, next=(curr+1)%nr) {
    if(inside[curr]>0 && inside[next]>0)            //inside->inside
      boundaryNew.push_back(f._boundary[next]);
    else if(inside[curr]>0 && inside[next]<=0) {    //inside->outside
      if(abs(inside[curr]-inside[next])<epsDist)   //safety check
        continue;
      T alpha=inside[curr]/(inside[curr]-inside[next]);
      boundaryNew.push_back(interp1D<Vec3T,T>(f._boundary[curr],f._boundary[next],alpha));
    } else if(inside[curr]<=0 && inside[next]>0) {  //outside->inside
      if(abs(inside[curr]-inside[next])<epsDist)   //safety check
        continue;
      T alpha=inside[curr]/(inside[curr]-inside[next]);
      boundaryNew.push_back(interp1D<Vec3T,T>(f._boundary[curr],f._boundary[next],alpha));
      boundaryNew.push_back(f._boundary[next]);
    } //outside->outside
  }
  f._boundary.clear();
  for(int i=0; i<boundaryNew.size(); i++)
    f._boundary.push_back(boundaryNew[i]);
}

template <typename T>
void SAT<T>::clip(Facet<T>& f,const Facet<T>& ref) {
  int nr=(int)ref._boundary.size();
  for(int i=0; i<nr; i++) {
    Vec3T dir=ref._boundary[(i+1)%nr]-ref._boundary[i];
    clip(f,ref._boundary[i],ref._n.cross(dir));
    if(f._boundary.empty())
      return;
  }
}
}
#endif
