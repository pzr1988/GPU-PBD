#include "XPBD.h"

namespace GPUPBD {
template <typename T>
XPBD<T>::XPBD(std::shared_ptr<Geometry<T>> geometry,T dt,int nRelax)
  :_geometry(geometry),_detector(new CollisionDetector<T>(geometry)),_dt(dt),_nRelax(nRelax) {}
template <typename T>
void XPBD<T>::step() {
  integrate();
  _detector->detectCollisions();
  for(int i=0; i<_nRelax; i++)
    relaxConstraint();
  updateVelocity();
}
template <typename T>
void XPBD<T>::integrate() {
  T dt = _dt;
  auto bIter = _geometry->getMutableCapsules().begin();
  auto eIter = _geometry->getMutableCapsules().end();
  thrust::for_each(thrust::device, bIter, eIter, [=] __host__ __device__ (Capsule<T>& capsule) {
    capsule._v += capsule._force*dt/capsule._mass;
    capsule._xNext = capsule._x+capsule._v*dt;
    capsule._R = capsule._q.toRotationMatrix();
    capsule._Iinv = capsule._R*capsule._Ibodyinv*capsule._R.transpose();
    capsule._w += capsule._Iinv*(capsule._torque
                                 - capsule._w.cross(capsule._R*capsule._Ibody*capsule._R.transpose()*capsule._w))*dt;
    Eigen::Quaternion<T> wQuat(0,capsule._w.x(),capsule._w.y(),capsule._w.z());
    Eigen::Quaternion<T> updatedQuat = Eigen::Quaternion<T>(0.5*dt,0,0,0)*wQuat*capsule._q;
    capsule._qNext = Eigen::Quaternion<T>(capsule._q.coeffs() + updatedQuat.coeffs());
    capsule._qNext.normalize();
  });
  cudaDeviceSynchronize();
}
template <typename T>
void XPBD<T>::relaxConstraint() {
  //to be implemeneted
}
template <typename T>
void XPBD<T>::updateVelocity() {
  T dt = _dt;
  auto bIter = _geometry->getMutableCapsules().begin();
  auto eIter = _geometry->getMutableCapsules().end();
  thrust::for_each(thrust::device, bIter, eIter, [=] __host__ __device__ (Capsule<T>& capsule) {
    capsule._v = (capsule._xNext-capsule._x)/dt;
    capsule._x = capsule._xNext;
    auto deltaQ = capsule._qNext*capsule._q.inverse();
    capsule._w = 2*deltaQ.vec()/dt;
    capsule._w = deltaQ.w() >=0 ? capsule._w : -capsule._w;
    capsule._q = capsule._qNext;
  });
  cudaDeviceSynchronize();
}
template <typename T>
const CollisionDetector<T>& XPBD<T>::getDetector() const {
  if (!_detector) {
    throw std::runtime_error("Detector is not initialized");
  }
  return *_detector;
}
//declare instance
template struct XPBD<LSCALAR>;
}
