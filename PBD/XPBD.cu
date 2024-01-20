#include "XPBD.h"
#include <thrust/sort.h>

namespace GPUPBD {
template <typename T>
XPBD<T>::XPBD(std::shared_ptr<Geometry<T>> geometry,T dt,int nRelax)
  :_geometry(geometry),_detector(new CollisionDetector<T>(geometry)),_dt(dt),_nRelax(nRelax),_collisionGroupAssigned(false) {}
template <typename T>
void XPBD<T>::step() {
  assignCollisionGroup();
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
  thrust::for_each(thrust::device, _geometry->begin(), _geometry->end(), [=] __host__ __device__ (Capsule<T>& capsule) {
    capsule._xPrev = capsule._x;
    capsule._qPrev = capsule._q;
    if(capsule._isDynamic) {
      capsule._v += capsule._force*dt/capsule._mass;
      capsule._x += capsule._v*dt;
      capsule._R = capsule._q.toRotationMatrix();
      capsule._Iinv = capsule._R*capsule._Ibodyinv*capsule._R.transpose();
      capsule._w += capsule._Iinv*(capsule._torque - capsule._w.cross(capsule._R*capsule._Ibody*capsule._R.transpose()*capsule._w))*dt;
      QuatT wQuat(0,capsule._w.x(),capsule._w.y(),capsule._w.z());
      QuatT updatedQuat = QuatT(0.5f*dt,0,0,0)*wQuat*capsule._q;
      capsule._q = QuatT(capsule._q.coeffs() + updatedQuat.coeffs());
      capsule._q.normalize();
    }
  });
  cudaDeviceSynchronize();
}
template <typename T>
void XPBD<T>::initRelaxConstraint() {
  if(numConstraints() == 0)
    return;
  _constraintCapsuleId.resize(numConstraints()*2); //each collision contains 2 capsules
  _update.resize(numConstraints()*2);
  _reduceCapsuleId.resize(numConstraints()*2);
  _reduceUpdate.resize(numConstraints()*2);
}
template <typename T>
void XPBD<T>::relaxConstraint() {
  if(numConstraints() == 0)
    return;
  size_t numJoints = _joints.size();
  Constraint<T>* d_joints = thrust::raw_pointer_cast(_joints.data());
  Constraint<T>* d_collisions = thrust::raw_pointer_cast(_detector->getCollisions());
  Capsule<T>* d_capsules = thrust::raw_pointer_cast(_geometry->getCapsules());
  T* d_lambda = thrust::raw_pointer_cast(_lambda.data());
  int* d_constraintCapsuleId = thrust::raw_pointer_cast(_constraintCapsuleId.data());
  auto d_update = thrust::raw_pointer_cast(_update.data());
  T dt = _dt;
  thrust::for_each(thrust::device,
                   thrust::make_counting_iterator(0),
                   thrust::make_counting_iterator((int)numConstraints()),
  [=] __host__ __device__ (int idx) {
    Vec3T pulse;
    auto& constraint = idx < numJoints ? d_joints[idx] : d_collisions[idx];
    auto& cA = d_capsules[constraint._capsuleIdA];
    auto& cB = d_capsules[constraint._capsuleIdB];
    auto placementPointA = cA._q.toRotationMatrix()*constraint._localPointA;
    auto placementPointB = cB._q.toRotationMatrix()*constraint._localPointB;
    auto globalPointA = placementPointA+cA._x;
    auto globalPointB = placementPointB+cB._x;
    if(constraint._type == Joint) {
      auto distSqr = (globalPointA - globalPointB).squaredNorm();
      if(distSqr > epsDist * epsDist)
        constraint._globalNormal = (globalPointA - globalPointB) / sqrt(distSqr);
      else constraint._globalNormal.setZero();
    }
    auto wA = computeGeneralizedInversMass(cA, placementPointA, constraint._globalNormal);
    auto wB = computeGeneralizedInversMass(cB, placementPointB, constraint._globalNormal);
    auto constraintViolation = (globalPointA-globalPointB).dot(constraint._globalNormal);
    auto alpha = constraint._alpha/(dt*dt);
    auto deltaLambda = (-constraintViolation-d_lambda[idx]*alpha)/(wA+wB+alpha);
    d_lambda[idx] += deltaLambda;
    pulse = deltaLambda * constraint._globalNormal;
    if(constraint._type == Collision && constraintViolation <= 0)
      pulse.setZero();
    // To avoid multi write problem, first cache update
    d_constraintCapsuleId[2*idx] = constraint._capsuleIdA;
    d_constraintCapsuleId[2*idx+1] = constraint._capsuleIdB;
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
  thrust::for_each(thrust::device, _geometry->begin(), _geometry->end(), [=] __host__ __device__ (Capsule<T>& capsule) {
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
  if (!_detector)
    throw std::runtime_error("Detector is not initialized");
  return *_detector;
}
template <typename T>
size_t XPBD<T>::numConstraints() const {
  return _joints.size() + _detector->size();
}
template <typename T>
void XPBD<T>::addJoint(size_t idA, size_t idB, const Vec3T& localA, const Vec3T& localB) {
  if(idA==idB)
    throw std::runtime_error("Cannot add joint to the same Capsule<T>!");
  Constraint<T> c;
  c._isValid=true;
  c._type=Joint;
  c._capsuleIdA=(int)idA;
  c._capsuleIdB=(int)idB;
  c._localPointA=localA;
  c._localPointB=localB;
  _joints.push_back(c);
  _collisionGroupAssigned=false;    //need to re-assign
}
template <typename T>
void XPBD<T>::assignCollisionGroup() {
  if(_collisionGroupAssigned)
    return;
  _collisionGroupAssigned=true;
}
template <typename T>
DEVICE_HOST T XPBD<T>::computeGeneralizedInversMass(const Capsule<T>& c, const Vec3T& n, const Vec3T& r) {
  if(!c._isDynamic)
    return 0;
  auto Iinv = c.getInertiaTensorInv();
  auto rCrossN = r.cross(n);
  auto w = 1.0f/c._mass+rCrossN.transpose()*Iinv*rCrossN;
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
  Capsule<T>* d_capsules = thrust::raw_pointer_cast(_geometry->getCapsules());
  //Reduce multi collisions of one capsule, then write
  thrust::sort_by_key(_constraintCapsuleId.begin(), _constraintCapsuleId.end(), _update.begin(), thrust::greater<int>());
  auto end = thrust::reduce_by_key(_constraintCapsuleId.begin(), _constraintCapsuleId.end(),
                                   _update.begin(), _reduceCapsuleId.begin(), _reduceUpdate.begin(),
                                   thrust::equal_to<int>(), UpdateAdd());
  _reduceCapsuleId.erase(end.first, _reduceCapsuleId.end());
  int * d_reduceCapsuleId = thrust::raw_pointer_cast(_reduceCapsuleId.data());
  auto d_reduceUpdate = thrust::raw_pointer_cast(_reduceUpdate.data());
  thrust::for_each(thrust::device,
                   thrust::make_counting_iterator(0),
                   thrust::make_counting_iterator((int)_reduceCapsuleId.size()),
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
