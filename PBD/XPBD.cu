#include "Pragma.h"
#include "XPBD.h"
#include <cmath>
#include <thrust/sort.h>
#include <unordered_map>

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
  for(int i=0; i<_nRelax; i++) {
    updateJointConstraint();
    relaxConstraint();
  }
  updateVelocity();
}
template <typename T>
void XPBD<T>::integrate() {
  T dt = _dt;
  thrust::for_each(thrust::device, _geometry->begin(), _geometry->end(), [=] __host__ __device__ (Shape<T>& shape) {
    shape._xPrev = shape._x;
    shape._qPrev = shape._q;
    if(shape._isDynamic) {
      shape._v += shape._force*dt/max(epsDir,shape._mass);
      shape._x += shape._v*dt;
      shape._R = shape._q.toRotationMatrix();
      shape._Iinv = shape._R*shape._Ibodyinv*shape._R.transpose();
      shape._w += shape._Iinv*(shape._torque - shape._w.cross(shape._R*shape._Ibody*shape._R.transpose()*shape._w))*dt;
      QuatT wQuat(0,shape._w.x(),shape._w.y(),shape._w.z());
      QuatT updatedQuat = QuatT(0.5f*dt,0,0,0)*wQuat*shape._q;
      shape._q = QuatT(shape._q.coeffs() + updatedQuat.coeffs());
      shape._q.normalize();
    }
  });
  cudaDeviceSynchronize();
}
template <typename T>
void XPBD<T>::initRelaxConstraint() {
  if(numConstraints() == 0)
    return;
  _lambda.clear();
  _lambda.resize(numConstraints());
  _constraintShapeId.resize(numConstraints()*2); //each collision contains 2 shapes
  _update.resize(numConstraints()*2);
  _reduceShapeId.resize(numConstraints()*2);
  _reduceUpdate.resize(numConstraints()*2);
}
template <typename T>
void XPBD<T>::relaxConstraint() {
  if(numConstraints() == 0)
    return;
  size_t numJointPositions = _jointPositions.size();
  size_t numJointAngulars = _jointAngulars.size();
  Constraint<T>* d_jointPositions = thrust::raw_pointer_cast(_jointPositions.data());
  Constraint<T>* d_jointAngulars = thrust::raw_pointer_cast(_jointAngulars.data());
  Constraint<T>* d_collisions = thrust::raw_pointer_cast(_detector->getCollisions());
  Shape<T>* d_shapes = thrust::raw_pointer_cast(_geometry->getShapes());
  T* d_lambda = thrust::raw_pointer_cast(_lambda.data());
  int* d_constraintShapeId = thrust::raw_pointer_cast(_constraintShapeId.data());
  auto d_update = thrust::raw_pointer_cast(_update.data());
  T dt = _dt;
  thrust::for_each(thrust::device,
                   thrust::make_counting_iterator(0),
                   thrust::make_counting_iterator((int)numConstraints()),
  [=] __host__ __device__ (int idx) {
    Vec3T pulse;
    auto& constraint = idx < numJointPositions ? d_jointPositions[idx] :
                       (idx < numJointPositions + numJointAngulars ?
                        d_jointAngulars[idx - numJointPositions] :
                        d_collisions[idx - numJointPositions - numJointAngulars]);
    auto& cA = d_shapes[constraint._shapeIdA];
    auto& cB = d_shapes[constraint._shapeIdB];
    auto placementPointA = cA._q.toRotationMatrix()*constraint._localPointA;
    auto placementPointB = cB._q.toRotationMatrix()*constraint._localPointB;
    auto globalPointA = placementPointA+cA._x;
    auto globalPointB = placementPointB+cB._x;
    auto wA = constraint._type == ConstraintType::JointAngular ? computeGeneralizedInversMass(cA, constraint._axis)
              : computeGeneralizedInversMass(cA, placementPointA, constraint._globalNormal);
    auto wB = constraint._type == ConstraintType::JointAngular ? computeGeneralizedInversMass(cB, constraint._axis)
              : computeGeneralizedInversMass(cB, placementPointB, constraint._globalNormal);
    auto constraintViolation = constraint._type == ConstraintType::JointAngular ? constraint._theta
                               : (globalPointA-globalPointB).dot(constraint._globalNormal);
    auto alpha = constraint._alpha / (dt*dt);
    auto deltaLambda = (-constraintViolation-d_lambda[idx]*alpha)/max(epsDir,(wA+wB+alpha));
    d_lambda[idx] += deltaLambda;
    pulse = constraint._type == ConstraintType::JointAngular ? deltaLambda * constraint._axis
            : deltaLambda * constraint._globalNormal;
    if(constraint._type == ConstraintType::Collision && constraintViolation <= 0)
      pulse.setZero();
    // To avoid multi write problem, first cache update
    d_constraintShapeId[2*idx] = constraint._shapeIdA;
    d_constraintShapeId[2*idx+1] = constraint._shapeIdB;
    if(constraint._type == ConstraintType::JointAngular) {
      d_update[2*idx]._x.setZero();
      d_update[2*idx+1]._x.setZero();
      d_update[2*idx]._q = getDeltaRot(cA, pulse).coeffs();
      d_update[2*idx+1]._q = -getDeltaRot(cB, pulse).coeffs();
    } else {
      d_update[2*idx]._x = pulse/max(epsDir,cA._mass);
      d_update[2*idx+1]._x = -pulse/max(epsDir,cB._mass);
      d_update[2*idx]._q = getDeltaRot(cA, placementPointA, pulse).coeffs();
      d_update[2*idx+1]._q = -getDeltaRot(cB, placementPointB, pulse).coeffs();
    }
  });
  updateShapeState();
}
template <typename T>
void XPBD<T>::updateJointConstraint() {
  Constraint<T>* d_jointPositions = thrust::raw_pointer_cast(_jointPositions.data());
  Shape<T>* d_shapes = thrust::raw_pointer_cast(_geometry->getShapes());
  thrust::for_each(thrust::device,
                   _jointPositions.begin(),
                   _jointPositions.end(),
  [=]  __device__ (Constraint<T>& constraint) {
    auto& cA = d_shapes[constraint._shapeIdA];
    auto& cB = d_shapes[constraint._shapeIdB];
    auto globalPointA = cA._q.toRotationMatrix()*constraint._localPointA+cA._x;
    auto globalPointB = cB._q.toRotationMatrix()*constraint._localPointB+cB._x;
    auto distSqr = (globalPointA - globalPointB).squaredNorm();
    if(distSqr > epsDist * epsDist)
      constraint._globalNormal = (globalPointA-globalPointB)/max(epsDir,sqrt(distSqr));
    else
      constraint._globalNormal.setZero();
  });
  thrust::for_each(thrust::device,
                   _jointAngulars.begin(),
                   _jointAngulars.end(),
  [=]  __device__ (Constraint<T>& constraint) {
    auto& cA = d_shapes[constraint._shapeIdA];
    auto& cB = d_shapes[constraint._shapeIdB];
    auto deltaQ = (constraint._targetQ.conjugate()*cA._q*cB._q.conjugate()).vec();
    // auto deltaQ = (constraint._targetQ.conjugate()*cB._q*cA._q.conjugate()).vec();
    auto len = sqrt(deltaQ.squaredNorm());
    constraint._theta = 2*asin(len);
    if(constraint._theta > epsDir)
      constraint._axis = deltaQ/max(epsDir,len);
    else {
      constraint._axis=Vec3T(0,0,1);
      constraint._theta=0;
    }
  });
}
template <typename T>
void XPBD<T>::updateVelocity() {
  T dt = _dt;
  thrust::for_each(thrust::device, _geometry->begin(), _geometry->end(), [=] __host__ __device__ (Shape<T>& shape) {
    if(shape._isDynamic) {
      shape._v = (shape._x-shape._xPrev)/dt;
      auto deltaQ = shape._q*shape._qPrev.conjugate();
      shape._w = 2*deltaQ.vec()/dt;
      shape._w = deltaQ.w() >=0 ? shape._w : -shape._w;
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
  return _detector->size() + _jointPositions.size() + _jointAngulars.size();
}
template <typename T>
void XPBD<T>::addJoint(size_t idA, size_t idB, const Vec3T& localA, const Vec3T& localB) {
  if(idA==idB)
    throw std::runtime_error("Cannot add joint to the same Shape<T>!");
  Constraint<T> c;
  c._isValid=true;
  c._type=ConstraintType::JointPosition;
  c._shapeIdA=(int)idA;
  c._shapeIdB=(int)idB;
  c._localPointA=localA;
  c._localPointB=localB;
  _jointPositions.push_back(c);
  _collisionGroupAssigned=false;    //need to re-assign
}
template <typename T>
void XPBD<T>::addJointAngular(size_t idA, size_t idB, const XPBD<T>::QuatT& targetQ) {
  if(idA==idB)
    throw std::runtime_error("Cannot add joint to the same Shape<T>!");
  Constraint<T> c;
  c._isValid=true;
  c._type=ConstraintType::JointAngular;
  c._shapeIdA=(int)idA;
  c._shapeIdB=(int)idB;
  c._targetQ=targetQ;
  _jointAngulars.push_back(c);
}
template <typename T>
void XPBD<T>::assignCollisionGroup() {
  if(_collisionGroupAssigned)
    return;
  //Step 1: do label propagation in shapes.
  thrust::device_vector<int> _parents(_geometry->size());
  thrust::device_vector<int> _changes(_geometry->size());
  thrust::sequence(thrust::device, _parents.begin(), _parents.end());
  int* d_parents = thrust::raw_pointer_cast(_parents.data());
  bool changed;
  do {
    thrust::transform(thrust::device, _jointPositions.begin(), _jointPositions.end(), _changes.begin(),
    [d_parents] __device__ (const Constraint<T>& c) {
      if(c._type!=ConstraintType::JointPosition || !c._isValid)
        return false;
      int parentA = d_parents[c._shapeIdA];
      int parentB = c._shapeIdB;
      if(parentA != parentB) {
        atomicExch(&d_parents[c._shapeIdA], parentB);
        return true;
      }
      return false;
    });
    int change_count = thrust::count_if(thrust::device, _changes.begin(), _changes.end(), thrust::identity<bool>());
    changed = change_count > 0;
  } while (changed);
  //Step 2: assign collision group id
  Shape<T>* d_shapes = thrust::raw_pointer_cast(_geometry->getShapes());
  thrust::for_each(thrust::device,
                   thrust::make_counting_iterator(0),
                   thrust::make_counting_iterator((int)_geometry->size()),
  [d_shapes, d_parents] __device__ (int idx) {
    d_shapes[idx]._parent = d_parents[idx];
  });
  _collisionGroupAssigned=true;
}
//Position Constraint
template <typename T>
DEVICE_HOST T XPBD<T>::computeGeneralizedInversMass(const Shape<T>& c, const Vec3T& n, const Vec3T& r) {
  if(!c._isDynamic)
    return 0;
  auto Iinv = c.getInertiaTensorInv();
  auto rCrossN = r.cross(n);
  auto w = 1.0f/max(epsDir,c._mass)+rCrossN.transpose()*Iinv*rCrossN;
  return w;
}
template <typename T>
DEVICE_HOST typename XPBD<T>::QuatT XPBD<T>::getDeltaRot(const Shape<T>& c, const Vec3T& r, const Vec3T& pulse) {
  auto cIinv = c.getInertiaTensorInv();
  auto cIinvRCrossP = cIinv * (r.cross(pulse)); // I^{-1}(r x p)
  QuatT cIinvRCrossPQuat(0,cIinvRCrossP.x(),cIinvRCrossP.y(),cIinvRCrossP.z());
  auto qUpdated = QuatT(0.5,0,0,0)*cIinvRCrossPQuat*c._q;
  return qUpdated;
}
//Angular Constraint
template <typename T>
DEVICE_HOST T XPBD<T>::computeGeneralizedInversMass(const Shape<T>& c, const Vec3T& n) {
  if(!c._isDynamic)
    return 0;
  auto Iinv = c.getInertiaTensorInv();
  auto w = n.transpose()*Iinv*n;
  return w;
}
template <typename T>
DEVICE_HOST typename XPBD<T>::QuatT XPBD<T>::getDeltaRot(const Shape<T>& c, const Vec3T& pulse) {
  auto cIinv = c.getInertiaTensorInv();
  auto cIinvP = cIinv * pulse;
  QuatT cIinvPQuat(0,cIinvP.x(),cIinvP.y(),cIinvP.z());
  auto qUpdated = QuatT(0.5,0,0,0)*cIinvPQuat*c._q;
  return qUpdated;
}
template <typename T>
void XPBD<T>::updateShapeState() {
  Shape<T>* d_shapes = thrust::raw_pointer_cast(_geometry->getShapes());
  //Reduce multi collisions of one shape, then write
  thrust::sort_by_key(_constraintShapeId.begin(), _constraintShapeId.end(), _update.begin(), thrust::greater<int>());
  auto end = thrust::reduce_by_key(_constraintShapeId.begin(), _constraintShapeId.end(),
                                   _update.begin(), _reduceShapeId.begin(), _reduceUpdate.begin(),
                                   thrust::equal_to<int>(), UpdateAdd());
  _reduceShapeId.erase(end.first, _reduceShapeId.end());
  int * d_reduceShapeId = thrust::raw_pointer_cast(_reduceShapeId.data());
  auto d_reduceUpdate = thrust::raw_pointer_cast(_reduceUpdate.data());
  thrust::for_each(thrust::device,
                   thrust::make_counting_iterator(0),
                   thrust::make_counting_iterator((int)_reduceShapeId.size()),
  [=] __host__ __device__ (int idx) {
    if(d_shapes[d_reduceShapeId[idx]]._isDynamic) {
      d_shapes[d_reduceShapeId[idx]]._x = d_shapes[d_reduceShapeId[idx]]._x + d_reduceUpdate[idx]._x;
      d_shapes[d_reduceShapeId[idx]]._q = QuatT(d_shapes[d_reduceShapeId[idx]]._q.coeffs() + d_reduceUpdate[idx]._q);
      d_shapes[d_reduceShapeId[idx]]._q.normalize();
    }
  });
}
template <typename T>
const thrust::device_vector<Constraint<T>>& XPBD<T>::getJointPositions() const {
  return _jointPositions;
}

//declare instance
template struct XPBD<LSCALAR>;
}
