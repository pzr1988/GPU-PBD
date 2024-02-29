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
  static DEVICE_HOST Vec2T project(const Shape<T>* S,const Trans<T>& trans,const Vec3T& n);
  template<int F1, int F2, int E1, int E2>
  static DEVICE_HOST ProjRange runSAT(const Shape<T>* A,
                                      const Shape<T>* B,
                                      const FixedVector<Facet<T>, F1>& FA,
                                      const FixedVector<Facet<T>, F2>& FB,
                                      const FixedVector<Edge<T>, E1>& EA,
                                      const FixedVector<Edge<T>, E2>& EB,
                                      const Trans<T>& transA,
                                      const Trans<T>& transB,bool* intersect=NULL);
  template<int F1, int F2, int E1, int E2>
  static DEVICE_HOST ProjRange generateManifold(const Shape<T>* A,
      const Shape<T>* B,
      const FixedVector<Facet<T>, F1>& FA,
      const FixedVector<Facet<T>, F2>& FB,
      const FixedVector<Edge<T>, E1>& EA,
      const FixedVector<Edge<T>, E2>& EB,
      const Trans<T>& transA,
      const Trans<T>& transB,bool* Intersect=NULL);
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
typename SAT<T>::Vec2T SAT<T>::project(const Shape<T>* S,const Trans<T>& trans,const Vec3T& n) {
  Vec2T ret=S->project(ROT(trans).transpose()*n);
  ret=(ret.array()+n.dot(CTR(trans))).matrix();
  sort2(ret[0],ret[1]);
  return ret;
}
template <typename T>
template <int F1, int F2, int E1, int E2>
typename SAT<T>::ProjRange SAT<T>::runSAT(const Shape<T>* A,
    const Shape<T>* B,
    const FixedVector<Facet<T>, F1>& FA,
    const FixedVector<Facet<T>, F2>& FB,
    const FixedVector<Edge<T>, E1>& EA,
    const FixedVector<Edge<T>, E2>& EB,
    const Trans<T>& transA,
    const Trans<T>& transB,bool* intersect) {
  ProjRange minRng,rng;
  minRng._fidA=minRng._eidA=minRng._fidB=minRng._eidB=-1;
  minRng._depth=std::numeric_limits<double>::max();
  if(intersect)
    *intersect=false;
  //faceA
  rng._fidA=rng._eidA=rng._fidB=rng._eidB=-1;
  for(int i=0; i<(int)FA.size(); i++) {
    rng._fidA=i;
    rng._n=ROT(transA)*FA[i]._n;
    rng._n=rng._n.template cast<double>().normalized().template cast<T>();
    rng._rngA=project(A,transA,rng._n);
    rng._rngB=project(B,transB,rng._n);
    rng._depth=depth(rng);
    if(rng._depth<=0) {
      if(intersect)
        *intersect=false;
      minRng._fidA=minRng._eidA=minRng._fidB=minRng._eidB=-1;
      return minRng;
    } else if(rng._depth<minRng._depth) {
      Vec3T cA=ROT(transA)*FA[i]._boundary[0]+CTR(transA);
      if(!matchWitness(rng,cA.dot(rng._n),(T*)NULL))
        continue;
      if(intersect)
        *intersect=true;
      minRng=rng;
    }
  }
  //faceB
  rng._fidA=rng._eidA=rng._fidB=rng._eidB=-1;
  for(int i=0; i<(int)FB.size(); i++) {
    rng._fidB=i;
    rng._n=ROT(transB)*FB[i]._n;
    rng._n=rng._n.template cast<double>().normalized().template cast<T>();
    rng._rngA=project(A,transA,rng._n);
    rng._rngB=project(B,transB,rng._n);
    rng._depth=depth(rng);
    if(rng._depth<=0) {
      if(intersect)
        *intersect=false;
      minRng._fidA=minRng._eidA=minRng._fidB=minRng._eidB=-1;
      return minRng;
    } else if(rng._depth<minRng._depth) {
      Vec3T cB=ROT(transB)*FB[i]._boundary[0]+CTR(transB);
      if(!matchWitness(rng,(T*)NULL,cB.dot(rng._n)))
        continue;
      if(intersect)
        *intersect=true;
      minRng=rng;
    }
  }
  //edgeAB
  rng._fidA=rng._eidA=rng._fidB=rng._eidB=-1;
  for(rng._eidA=0; rng._eidA<(int)EA.size(); rng._eidA++) {
    Vec3T cA1=ROT(transA)*EA[rng._eidA]._a+CTR(transA);
    Vec3T cA2=ROT(transA)*EA[rng._eidA]._b+CTR(transA);
    Vec3T dA=ROT(transA)*(EA[rng._eidA]._b-EA[rng._eidA]._a);
    for(rng._eidB=0; rng._eidB<(int)EB.size(); rng._eidB++) {
      rng._n=dA.cross(ROT(transB)*(EB[rng._eidB]._b-EB[rng._eidB]._a));
      T nLenSqr=rng._n.squaredNorm();
      if(nLenSqr<epsDist)
        continue;
      //project range with witness, to ensure the closest point pair lies on the line segment
      rng._n/=sqrt((double)nLenSqr);
      rng._rngA=project(A,transA,rng._n);
      rng._rngB=project(B,transB,rng._n);
      rng._depth=depth(rng);
      if(rng._depth<=0) {
        if(intersect)
          *intersect=false;
        minRng._fidA=minRng._eidA=minRng._fidB=minRng._eidB=-1;
        return minRng;
      } else if(rng._depth<minRng._depth-epsDist) { //we prefer face intersection
        //edge-edge case, only a single contact is generated
        Vec3T cB1=ROT(transB)*EB[rng._eidB]._a+CTR(transB);
        Vec3T cB2=ROT(transB)*EB[rng._eidB]._b+CTR(transB);
        //find intersection coordinates
        Mat3X2T LHS;
        LHS.col(0)=-(cA2-cA1);
        LHS.col(1)= (cB2-cB1);
        Vec2T bary=(LHS.transpose()*LHS).inverse()*(LHS.transpose()*(cA1-cB1));
        if((bary.array()<0).any() || (bary.array()>1).any())
          continue;
        if(!matchWitness(rng,cA1.dot(rng._n),cB1.dot(rng._n)))
          continue;
        if(intersect)
          *intersect=true;
        //only (projected) intersecting line segments form potential closest point pair
        rng._ptA=cA1*(1-bary[0])+cA2*bary[0];
        rng._ptB=cB1*(1-bary[1])+cB2*bary[1];
        minRng=rng;
      }
    }
  }
  return minRng;
}
//not only run contact, but also generate manifold
template <typename T>
template <int F1, int F2, int E1, int E2>
typename SAT<T>::ProjRange SAT<T>::generateManifold(const Shape<T>* A,
    const Shape<T>* B,
    const FixedVector<Facet<T>, F1>& FA,
    const FixedVector<Facet<T>, F2>& FB,
    const FixedVector<Edge<T>, E1>& EA,
    const FixedVector<Edge<T>, E2>& EB,
    const Trans<T>& transA,
    const Trans<T>& transB,bool* Intersect) {
  bool intersect;
  ProjRange rng=runSAT<F1,F2,E1,E2>(A,B,FA,FB,EA,EB,transA,transB,&intersect);
  // if(Intersect)
  //   *Intersect=intersect;
  // if(!intersect)
  //   return rng;
  // if(rng._eidA>=0 && rng._eidB>=0) {
  //   //edge-edge case, only a single contact is generated
  //   ContactPoint p;
  //   p._ptA=rng._ptA;
  //   p._ptB=rng._ptB;
  //   p._nA2B=rng._n;
  //   if(p.depth()<0)
  //     p._nA2B*=-1;
  //   m._points.push_back(p);
  // } else if(rng._fidA>=0) {
  //   generateManifoldFace(A,B,FA,FB,EA,EB,transA,transB,rng._fidA,m._points);
  // } else {
  //   int nrP=(int)m._points.size();
  //   generateManifoldFace(B,A,FB,FA,EB,EA,transB,transA,rng._fidB,m._points);
  //   for(int i=nrP; i<(int)m._points.size(); i++)
  //     m._points[i].swap();
  // }
  return rng;
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
