#ifndef SAT_H
#define SAT_H

#include "Pragma.h"
#include "Geometry.h"

namespace GPUPBD {
template <typename T,typename TF>
DEVICE_HOST T interp1D(const T& v0,const T& v1,
                       const TF& px) {
  return v0*(TF)(1.0f-px)+v1*(TF)px;
}

template <typename T>
struct SAT {
  DECL_MAT_VEC_MAP_TYPES_T
  static DEVICE_HOST void clip(Facet<T>& f,const Vec3T& pos,const Vec3T& inward);
  static DEVICE_HOST void clip(Facet<T>& f,const Facet<T>& ref);
};

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
