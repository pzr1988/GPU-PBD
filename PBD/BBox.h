#ifndef BBOX_H
#define BBOX_H

#include "Utils.h"

namespace GPUPBD {
struct DLL_EXPORT BBox {
  typedef LSCALAR T;
  DECL_MAT_VEC_MAP_TYPES_T
  BBox();
  BBox(const Vec3T& p);
  BBox(const Vec3T& minC,const Vec3T& maxC);
  BBox(const BBox& other);
  virtual ~BBox();
  template <typename T2>
  BBox& operator=(const BBox& other) {
    copy(other);
    return *this;
  }
  static BBox createMM(const Vec3T& minC,const Vec3T& maxC);
  static BBox createME(const Vec3T& minC,const Vec3T& extent);
  static BBox createCE(const Vec3T& center,const Vec3T& extent);
  BBox getIntersect(const BBox& other) const;
  BBox getUnion(const BBox& other) const;
  BBox getUnion(const Vec3T& point) const;
  BBox getUnion(const Vec3T& ctr,const T& rad) const;
  void setIntersect(const BBox& other);
  void setUnion(const BBox& other);
  void setUnion(const Vec3T& point);
  void setUnion(const Vec3T& ctr,const T& rad);
  void setPoints(const Vec3T& a,const Vec3T& b,const Vec3T& c);
  void setPoints(const Vec3T& a,const Vec3T& b,const Vec3T& c,const Vec3T& d);
  const Vec3T& minCorner() const;
  const Vec3T& maxCorner() const;
  void enlargedEps(T eps);
  BBox enlargeEps(T eps) const;
  void enlarged(T len);
  BBox enlarge(T len) const;
  Vec3T lerp(const Vec3T& frac) const;
  bool empty() const;
  template <int DIM2>
  bool containDim(const Vec3T& point) const {
    for(int i=0; i<DIM2; i++)
      if(_minC[i] > point[i] || _maxC[i] < point[i])
        return false;
    return true;
  }
  bool contain(const BBox& other) const;
  bool contain(const Vec3T& point) const;
  bool contain(const Vec3T& point,const T& rad) const;
  void reset();
  Vec3T getExtent() const;
  T distTo(const BBox& other) const;
  T distTo(const Vec3T& pt) const;
  T distToSqr(const Vec3T& pt) const;
  Vec3T closestTo(const Vec3T& pt) const;
  bool intersect(const Vec3T& p,const Vec3T& q) const;
  bool intersect(const Vec3T& p,const Vec3T& q,T& s,T& t) const;
  bool intersect(const BBox& other) const;
  Vec2T project(const Vec3T& a) const;
  BBox& copy(const BBox& other);
  T perimeter() const;
  Vec3T _minC,_maxC;
};
}

#endif
