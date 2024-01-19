#include "XPBD.h"
#include <thrust/sort.h>

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
    capsule._qPrev = capsule._q;
    if(capsule._isDynamic) {
      capsule._v += capsule._force*dt/capsule._mass;
      capsule._x += capsule._v*dt;
      capsule._R = capsule._q.toRotationMatrix();
      capsule._Iinv = capsule._R*capsule._Ibodyinv*capsule._R.transpose();
      capsule._w += capsule._Iinv*(capsule._torque
                                   - capsule._w.cross(capsule._R*capsule._Ibody*capsule._R.transpose()*capsule._w))*dt;
      QuatT wQuat(0,capsule._w.x(),capsule._w.y(),capsule._w.z());
      QuatT updatedQuat = QuatT(0.5*dt,0,0,0)*wQuat*capsule._q;
      capsule._q = QuatT(capsule._q.coeffs() + updatedQuat.coeffs());
      capsule._q.normalize();
    }
  });
  cudaDeviceSynchronize();
}
template <typename T>
void XPBD<T>::initRelaxConstraint() {
  int numCollision = _detector->size();
  if(numCollision == 0) {
    return;
  }
  _lambda.clear();
  _lambda.resize(numCollision);
  _collisionCapsuleId.resize(numCollision*2); //each collision contains 2 capsules
  _update.resize(numCollision*2);
  _reduceCapsuleId.resize(numCollision*2);
  _reduceUpdate.resize(numCollision*2);
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
  int* d_collisionCapsuleId = thrust::raw_pointer_cast(_collisionCapsuleId.data());
  auto d_update = thrust::raw_pointer_cast(_update.data());
  T dt = _dt;
  thrust::for_each(thrust::device,
                   thrust::make_counting_iterator(0),
                   thrust::make_counting_iterator(static_cast<int>(collisions.size())),
  [=] __host__ __device__ (int idx) {
    auto& collision = d_collisions[idx];
    auto& cA = d_capsules[collision._capsuleIdA];
    auto& cB = d_capsules[collision._capsuleIdB];
    auto placementPointA = cA._q.toRotationMatrix()*collision._localPointA;
    auto placementPointB = cB._q.toRotationMatrix()*collision._localPointB;
    auto globalPointA = placementPointA+cA._x;
    auto globalPointB = placementPointB+cB._x;
    auto wA = computeGeneralizedInversMass(cA, placementPointA,
                                           collision._globalNormal);
    auto wB = computeGeneralizedInversMass(cB, placementPointB,
                                           collision._globalNormal);
    auto collisionDepth = (globalPointA-globalPointB).dot(collision._globalNormal);
    auto alpha = collision._alpha/(dt*dt);
    auto deltaLambda = (-collisionDepth-d_lambda[idx]*alpha)/(wA+wB+alpha);
    d_lambda[idx] += deltaLambda;
    Vec3T pulse = deltaLambda*collision._globalNormal;
    if(collisionDepth<=0) {
      pulse = Vec3T(0,0,0);
    }
    // To avoid multi write problem, first cache update
    d_collisionCapsuleId[2*idx] = collision._capsuleIdA;
    d_collisionCapsuleId[2*idx+1] = collision._capsuleIdB;
    d_update[2*idx]._x = pulse/cA._mass;
    d_update[2*idx+1]._x = -pulse/cB._mass;
    d_update[2*idx]._q = getDeltaRot(cA, placementPointA, pulse).coeffs();
    d_update[2*idx+1]._q = -getDeltaRot(cB, placementPointB, pulse).coeffs();
  });
  updateCapsuleState();
}
template <typename T>
void XPBD<T>::updateVelocity() {
  T dt = _dt;
  auto bIter = _geometry->getMutableCapsules().begin();
  auto eIter = _geometry->getMutableCapsules().end();
  thrust::for_each(thrust::device, bIter, eIter, [=] __host__ __device__ (Capsule<T>& capsule) {
    if(capsule._isDynamic) {
      capsule._v = (capsule._x-capsule._xPrev)/dt;
      auto deltaQ = capsule._q*capsule._qPrev.inverse();
      capsule._w = 2*deltaQ.vec()/dt;
      capsule._w = deltaQ.w() >=0 ? capsule._w : -capsule._w;
    }
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
template <typename T>
DEVICE_HOST T XPBD<T>::computeGeneralizedInversMass(const Capsule<T>& c, const Vec3T& n, const Vec3T& r) {
  if(!c._isDynamic) {
    return .0;
  }
  auto Iinv = c.getInertiaTensorInv();
  auto rCrossN = r.cross(n);
  auto w = 1.0/c._mass+rCrossN.transpose()*Iinv*rCrossN;
  return w;
}
template <typename T>
DEVICE_HOST typename XPBD<T>::QuatT XPBD<T>::getDeltaRot(const Capsule<T>& c, const Vec3T& r, const Vec3T& pulse) {
  auto cIinv = c.getInertiaTensorInv();
  auto cIinvRCrossP = cIinv * (r.cross(pulse)); // I^{-1}(r x p)
  QuatT cIinvRCrossPQuat(0,cIinvRCrossP.x(),cIinvRCrossP.y(),cIinvRCrossP.z());
  auto qUpdated = QuatT(0.5,0,0,0)*cIinvRCrossPQuat*c._q;
  return qUpdated;
}
template <typename T>
void XPBD<T>::updateCapsuleState() {
  auto& capsules = _geometry->getMutableCapsules();
  Capsule<T>* d_capsules = thrust::raw_pointer_cast(capsules.data());
  //Reduce multi collisions of one capsule, then write
  thrust::sort_by_key(_collisionCapsuleId.begin(), _collisionCapsuleId.end(),
                      _update.begin(), thrust::greater<int>());
  auto end = thrust::reduce_by_key(_collisionCapsuleId.begin(), _collisionCapsuleId.end(),
                                   _update.begin(), _reduceCapsuleId.begin(), _reduceUpdate.begin(),
                                   thrust::equal_to<int>(), UpdateAdd());
  _reduceCapsuleId.erase(end.first, _reduceCapsuleId.end());
  int * d_reduceCapsuleId = thrust::raw_pointer_cast(_reduceCapsuleId.data());
  auto d_reduceUpdate = thrust::raw_pointer_cast(_reduceUpdate.data());
  thrust::for_each(thrust::device,
                   thrust::make_counting_iterator(0),
                   thrust::make_counting_iterator(static_cast<int>(_reduceCapsuleId.size())),
  [=] __host__ __device__ (int idx) {
    if(d_capsules[d_reduceCapsuleId[idx]]._isDynamic) {
      d_capsules[d_reduceCapsuleId[idx]]._x = d_capsules[d_reduceCapsuleId[idx]]._x + d_reduceUpdate[idx]._x;
      d_capsules[d_reduceCapsuleId[idx]]._q = QuatT(d_capsules[d_reduceCapsuleId[idx]]._q.coeffs() + d_reduceUpdate[idx]._q);
      d_capsules[d_reduceCapsuleId[idx]]._q.normalize();
    }
  });
}

//declare instance
template struct XPBD<LSCALAR>;
}
