#include "Geometry.h"

namespace GPUPBD {
template <typename T>
void Shape<T>::initInertiaTensor(T rho) {
  if(ShapeType::Capsule==_type) {
    //https://gamedev.net/tutorials/programming/math-and-physics/shape-inertia-tensor-r3856/
    _mass = 3.14f * (_radius*_radius*_len) + 3.14f * (4.f/3.f)*_radius*_radius*_radius;
    _mass *= rho;
    T mCy = _mass * (_radius*_radius*_len)/(_radius*_radius*_len + (4.f/3.f)*_radius*_radius*_radius);
    T mHSph = _mass * ((2.f/3.f)*_radius*_radius*_radius)/(_radius*_radius*_len + (4.f/3.f)*_radius*_radius*_radius);
    T Ixx = mCy*(_radius*_radius/2.f)
            + 2.f*mHSph*(2.f*_radius*_radius/5.f);
    T Iyy = mCy*(_radius*_radius/4.f + _len*_len/12.f)
            + 2.f*mHSph*(2.f*_radius*_radius/5.f+_len*_len/2.f+3.f*_len*_radius/8);
    _Ibody.setIdentity();
    _Ibody.diagonal()=Vec3T(Ixx,Iyy,Iyy);
  } else if(ShapeType::Box==_type) {
    _mass = _len*_width*_height;
    _mass *= rho;
    T Ixx = _mass*(_width*_width+_height*_height)/12.f;
    T Iyy = _mass*(_len*_len+_height*_height)/12.f;
    T Izz = _mass*(_len*_len+_width*_width)/12.f;
    _Ibody.setIdentity();
    _Ibody.diagonal()=Vec3T(Ixx,Iyy,Izz);
  } else if(ShapeType::Sphere==_type) {
    _mass = (4.f/3.f)*3.14f*_radius*_radius*_radius;
    _mass *= rho;
    T Ixx = (2.f/5.f)*_mass*_radius*_radius;
    _Ibody.diagonal()=Vec3T(Ixx,Ixx,Ixx);
    _Ibody.setIdentity();
  } else {
    _mass = 1.f;
    _Ibody.setIdentity();
  }
  _Ibodyinv=_Ibody.inverse();
}

template <typename T>
size_t Geometry<T>::size() const {
  return _nrShape;
}
template <typename T>
void Geometry<T>::resize(size_t nrShape) {
  _nrShape=nrShape;
  if(_nrShape>_shapes.size())
    _shapes.resize(_nrShape);
}
template <typename T>
void Geometry<T>::reserve(size_t nrShape) {
  _shapes.resize(std::max(nrShape,_nrShape));
}
template <typename T>
void Geometry<T>::setShape(const std::vector<Shape<T>>& c) {
  if(c.size()!=_nrShape)
    throw std::runtime_error("#Shape mismatch!");
  thrust::copy(c.begin(),c.end(),_shapes.begin());
}
template <typename T>
void Geometry<T>::setShape(size_t id, const Shape<T>& c) {
  _shapes[id]=c;
}
template <typename T>
Shape<T> Geometry<T>::operator[](size_t id) const {
  return _shapes[id];
}
template <typename T>
typename thrust::device_vector<Shape<T>>::iterator Geometry<T>::begin() {
  return _shapes.begin();
}
template <typename T>
typename thrust::device_vector<Shape<T>>::iterator Geometry<T>::end() {
  return _shapes.begin()+_nrShape;
}
template <typename T>
typename thrust::device_vector<Shape<T>>::const_iterator Geometry<T>::begin() const {
  return _shapes.begin();
}
template <typename T>
typename thrust::device_vector<Shape<T>>::const_iterator Geometry<T>::end() const {
  return _shapes.begin()+_nrShape;
}
template <typename T>
typename thrust::device_ptr<const Shape<T>> Geometry<T>::getShapes() const {
  return _shapes.data();
}
template <typename T>
typename thrust::device_ptr<Shape<T>> Geometry<T>::getShapes() {
  return _shapes.data();
}

//declare instance
template struct Shape<LSCALAR>;
template struct Geometry<LSCALAR>;
}
