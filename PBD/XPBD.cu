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
  const Collision<T>* d_collisions = thrust::raw_pointer_cast(collisions.data());
  auto& capsules = _geometry->getMutableCapsules();
  Capsule<T>* d_capsules = thrust::raw_pointer_cast(capsules.data());
  T* d_lambda = thrust::raw_pointer_cast(_lambda.data());
  T dt = _dt;
  thrust::for_each(thrust::device,
                   thrust::make_counting_iterator(0),
                   thrust::make_counting_iterator(static_cast<int>(collisions.size())),
  [=] __host__ __device__ (int idx) {
    auto& collision = d_collisions[idx];
    auto& cA = d_capsules[collision._capsuleIdA];
    auto& cB = d_capsules[collision._capsuleIdB];
    auto globalPointA = cA._q.toRotationMatrix()*collision._localPointA+cA._x;
    auto globalPointB = cB._q.toRotationMatrix()*collision._localPointB+cB._x;
    auto cARCrossN = globalPointA.cross(collision._globalNormal); //r_1 x n
    auto cAR = cA._q.toRotationMatrix();
    auto cAIinv = cAR*cA._Ibodyinv*cAR.transpose(); //inverse of interia tensor
    auto cBRCrossN = globalPointB.cross(collision._globalNormal); //r_2 x n
    auto cBR = cB._q.toRotationMatrix();
    auto cBIinv = cBR*cB._Ibodyinv*cBR.transpose(); //inverse of interia tensor
    T wA = 1.0/cA._mass+cARCrossN.transpose()*cAIinv*cARCrossN;
    T wB = 1.0/cB._mass+cBRCrossN.transpose()*cBIinv*cBRCrossN;
    T c = (globalPointA-globalPointB).dot(collision._globalNormal);
    T alpha= collision._alpha / (dt*dt);
    T deltaLambda = (-c-d_lambda[idx]*alpha)/(wA + wB + alpha);
    d_lambda[idx] += deltaLambda;
    Vec3T p = deltaLambda*collision._globalNormal;
    // TODO 原子
    cA._x += p/cA._mass;
    cB._x -= p/cB._mass;
    Vec3T cAIinvRCrossP = cAIinv * (globalPointA.cross(p)); // I^{-1}(r x p)
    Eigen::Quaternion<T> cAIinvRCrossPQuat(0,cAIinvRCrossP.x(),cAIinvRCrossP.y(),cAIinvRCrossP.z());
    Eigen::Quaternion<T> cAQUpdated = Eigen::Quaternion<T>(0.5,0,0,0)*cAIinvRCrossPQuat*cA._q;
    cA._q = Eigen::Quaternion<T>(cA._q.coeffs() + cAQUpdated.coeffs());
    Vec3T cBIinvRCrossP = cBIinv * (globalPointB.cross(p)); // I^{-1}(r x p)
    Eigen::Quaternion<T> cBIinvRCrossPQuat(0,cBIinvRCrossP.x(),cBIinvRCrossP.y(),cBIinvRCrossP.z());
    Eigen::Quaternion<T> cBQUpdated = Eigen::Quaternion<T>(0.5,0,0,0)*cBIinvRCrossPQuat*cB._q;
    cB._q = Eigen::Quaternion<T>(cB._q.coeffs() - cBQUpdated.coeffs());
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
