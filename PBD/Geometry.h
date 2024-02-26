#ifndef GEOMETRY_H
#define GEOMETRY_H
#include "Pragma.h"
#include "LBVH/aabb.cuh"
#include <thrust/device_vector.h>

namespace GPUPBD {
//This struct represents a shape-shaped object
//This struct represents a shaped object
//The capsule shape's centerline extends from (-_len/2,0,0) to (_len/2,0,0)
//The box shape's centerline extends from (-_len/2,-_width/2,-_height/2) to (_len/2,_width/2,_height/2)
//The radius of the capsule is stored in _radius
//The height and width of the box is stored in _height, _width
//The local to global transformation is stored in _x and _q
//The 3x3 rotational matrix is: _q.toRotationMatrix()
//The 3x1 translational vector is: _x
enum class ShapeType {
  Capsule,
  Box,
  Unknown,
};
template <typename T>
struct Shape {
  DECL_MAT_VEC_MAP_TYPES_T
  /* Constant quantities */
  ShapeType _type=ShapeType::Unknown;
  T _len;
  T _radius; //capsule
  T _width, _height; //box
  T _mass;
  Mat3T _Ibody, _Ibodyinv; //body space inertia tensor
  bool _isDynamic;
  /* State variables */
  Vec3T _x;
  QuatT _q;
  Vec3T _xPrev; //Auxilliary variable for XPBD
  QuatT _qPrev;
  Vec3T _v; //linear velocity
  Vec3T _w; //angular velocity
  /* Derived quantities (auxiliary variables) */
  Mat3T _Iinv; // inverse of inertia tensor
  Mat3T _R; //rotation matrix
  /* Computed quantities */
  Vec3T _force;
  Vec3T _torque;
  // Collision exclusion
  int _parent = -1;
  void initInertiaTensor(T rho=1);
  DEVICE_HOST Vec3T minCorner() const {
    if(ShapeType::Capsule==_type)
      return Vec3T(-_len/2, 0, 0);
    else if(ShapeType::Box==_type)
      return Vec3T(-_len/2,-_width/2,-_height/2);
    return Vec3T(0, 0, 0);
  }
  DEVICE_HOST Vec3T maxCorner() const {
    if(ShapeType::Capsule==_type)
      return Vec3T(_len/2, 0, 0);
    else if(ShapeType::Box==_type)
      return Vec3T(_len/2,_width/2,_height/2);
    return Vec3T(0, 0, 0);
  }
  DEVICE_HOST Vec3T globalMinCorner() const {
    return _q.toRotationMatrix()*minCorner()+_x;
  }
  DEVICE_HOST Vec3T globalMaxCorner() const {
    return _q.toRotationMatrix()*maxCorner()+_x;
  }
  DEVICE_HOST Mat3T getInertiaTensorInv() const {
    auto R = _q.toRotationMatrix();
    return R*_Ibodyinv*R.transpose();
  }
  DEVICE_HOST Vec3T support(const Vec3T& D,int& id) const {
    id=-1;
    Vec3T v,maxV;
    T dist,maxDist=-FLT_MAX;
    Vec3T _minC = minCorner();
    Vec3T _maxC = maxCorner();
    for(int i=0; i<8; i++) {
      v[0]=(i&1)?_maxC[0]:_minC[0];
      v[1]=(i&2)?_maxC[1]:_minC[1];
      v[2]=(i&4)?_maxC[2]:_minC[2];
      dist=v.dot(D);
      if(dist>maxDist) {
        maxDist=dist;
        maxV=v;
        id=i;
      }
    }
    return maxV;
  }
  DEVICE_HOST bool isCapsule() const {
    return ShapeType::Capsule==_type; 
  }
};

//The geometry stores a vector of shapes
//The vector is over-sized and pre-allocated
//The number of shapes in use is stored in _nrShape
template <typename T>
struct Geometry {
  DECL_MAT_VEC_MAP_TYPES_T
  //Get/Set the actual number of shapes used
  size_t size() const;
  void resize(size_t nrShape);
  //Set the pre-allocated shape list
  void reserve(size_t nrShape);
  //CPU->GPU transfer: copying the list of shapes to GPU
  void setShape(const std::vector<Shape<T>>& c);
  //Set the id-th shape
  void setShape(size_t id, const Shape<T>& c);
  //Get the id-th shape
  Shape<T> operator[](size_t id) const;
  //Get All Shapes
  typename thrust::device_vector<Shape<T>>::iterator begin();
  typename thrust::device_vector<Shape<T>>::iterator end();
  typename thrust::device_vector<Shape<T>>::const_iterator begin() const;
  typename thrust::device_vector<Shape<T>>::const_iterator end() const;
  typename thrust::device_ptr<const Shape<T>> getShapes() const;
  typename thrust::device_ptr<Shape<T>> getShapes();
 protected:
  thrust::device_vector<Shape<T>> _shapes;
  size_t _nrShape=0;
};

// General function to get AABB
template <template<typename> class Geometry, typename T>
struct AABBGetter;
// AABBGetter for Shape<T>
template <typename T>
struct AABBGetter<Shape, T> {
  DEVICE_HOST lbvh::aabb<T> operator()(const Shape<T> &c) const noexcept {
    lbvh::aabb<T> retval;
    auto transformedEnd1 = c.globalMaxCorner();
    auto transformedEnd2 = c.globalMinCorner();
    auto upper = transformedEnd1.cwiseMax(transformedEnd2);
    T radius = static_cast<T>(c._radius);
    retval.upper.x = upper.x() + radius;
    retval.upper.y = upper.y() + radius;
    retval.upper.z = upper.z() + radius;
    auto lower = transformedEnd1.cwiseMin(transformedEnd2);
    retval.lower.x = lower.x() - radius;
    retval.lower.y = lower.y() - radius;
    retval.lower.z = lower.z() - radius;
    return retval;
  }
};
}
#endif
