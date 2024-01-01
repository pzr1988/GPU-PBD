#include "XPBD.h"

namespace GPUPBD {
template <typename T>
XPBD<T>::XPBD(std::shared_ptr<Geometry<T>> geometry,T dt,int nRelax)
  :_geometry(geometry),_detector(new CollisionDetector<T>(geometry)),_dt(dt),_nRelax(nRelax) {}
template <typename T>
void XPBD<T>::step() {
  integrate();
  _detector->detectCollisions();
  initRelaxConstraint();
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
    capsule._xPrev = capsule._x;
    capsule._v += capsule._force*dt/capsule._mass;
    capsule._x += capsule._v*dt;

    capsule._qPrev = capsule._q;
    capsule._R = capsule._q.toRotationMatrix();
    capsule._Iinv = capsule._R*capsule._Ibodyinv*capsule._R.transpose();
    capsule._w += capsule._Iinv*(capsule._torque
                                 - capsule._w.cross(capsule._R*capsule._Ibody*capsule._R.transpose()*capsule._w))*dt;
    Eigen::Quaternion<T> wQuat(0,capsule._w.x(),capsule._w.y(),capsule._w.z());
    Eigen::Quaternion<T> updatedQuat = Eigen::Quaternion<T>(0.5*dt,0,0,0)*wQuat*capsule._q;
    capsule._q = Eigen::Quaternion<T>(capsule._q.coeffs() + updatedQuat.coeffs());
    capsule._q.normalize();
  });
  cudaDeviceSynchronize();
}
template <typename T>
void XPBD<T>::initRelaxConstraint() {
  if(_detector->size() == 0) {
    return;
  }
  _lambda.clear();
  _lambda.resize(_detector->size());
}
template <typename T>
void XPBD<T>::relaxConstraint() {
  if(_detector->size() == 0) {
    return;
  }
  const auto& collisions = _detector->getCollisions();
  auto& capsules = _geometry->getMutableCapsules();
  Capsule<T>* d_capsules = thrust::raw_pointer_cast(capsules.data());
  auto bIter = collisions.begin();
  auto eIter = collisions.end();

  thrust::for_each(thrust::device, bIter, eIter, [d_capsules] __host__ __device__ (const Collision<T>& collision) {
    // TODO
    int idA = collision._capsuleIdA;
    int idB = collision._capsuleIdB;

  });

}
template <typename T>
void XPBD<T>::updateVelocity() {
  T dt = _dt;
  auto bIter = _geometry->getMutableCapsules().begin();
  auto eIter = _geometry->getMutableCapsules().end();
  thrust::for_each(thrust::device, bIter, eIter, [=] __host__ __device__ (Capsule<T>& capsule) {
    capsule._v = (capsule._x-capsule._xPrev)/dt;
    auto deltaQ = capsule._q*capsule._qPrev.inverse();
    capsule._w = 2*deltaQ.vec()/dt;
    capsule._w = deltaQ.w() >=0 ? capsule._w : -capsule._w;
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
