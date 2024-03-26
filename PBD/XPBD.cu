#include "Pragma.h"
#include "XPBD.h"
#include <cmath>
#include <thrust/sort.h>
#include <unordered_map>

namespace GPUPBD {
template <typename T>
XPBD<T>::XPBD(std::shared_ptr<Geometry<T>> geometry,T dt,int nRelax)
  :_geometry(geometry),_detector(new CollisionDetector<T>(geometry)),_dt(dt),_nRelax(nRelax),_isPlay(false),_collisionGroupAssigned(false) {}
template <typename T>
void XPBD<T>::step() {
  assignCollisionGroup();
  integrate();
  playAnimation();
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
    // 推导方式1
    // deltaQ * cB._q = cB_target._q // deltaQ是对parent的修正量。
    // cB_target._q = cB world target q * sB
    // = cA world q * psQ.inverse() * sB
    // = cA._q * sA.inverse() * psQ.inverse() * sB
    // =>
    // deltaQ = cA._q * sA.inverse() * psQ.inverse() * sB * cB._q.inverse()
    // 推导方式2
    // cB._q = cB's world.q * sB  // sB是胶囊体在局部坐标系的旋转，cB's world.q 是局部坐标系的旋转， cB._q是胶囊体在全局坐标系的旋转
    // =》 cB_world.q = cB._q * sB.inverse()
    // deltaQ.inverse() * cA._q = cA_target._q // cA_target._q是胶囊体的目标旋转角度（全局坐标系下），deltaQ是rigid bodyXPBD公式14中的角度修正。
    // = cA_world._q * sA
    // = cB_world._q * psQ * sA
    // = cB._q * sB.inverse() * psQ * sA
    // =》
    // deltaQ.inverse() = cB._q * sB.inverse() * psQ * sA * cA._q.inverse()
    // =》
    // deltaQ = cA._q * sA.inverse() * psQ.inverse() * sB * cB._q.inverse()
    auto deltaQ = (cA._q*constraint._aQ.conjugate()*constraint._targetQ.conjugate()*constraint._bQ*cB._q.conjugate()).vec();
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
void XPBD<T>::addJoint(size_t idA, size_t idB, const Vec3T& localA, const Vec3T& localB, T alpha) {
  if(idA==idB)
    throw std::runtime_error("Cannot add joint to the same Shape<T>!");
  Constraint<T> c;
  c._isValid=true;
  c._alpha=alpha;
  c._type=ConstraintType::JointPosition;
  c._shapeIdA=(int)idA;
  c._shapeIdB=(int)idB;
  c._localPointA=localA;
  c._localPointB=localB;
  _jointPositions.push_back(c);
  _collisionGroupAssigned=false;    //need to re-assign
}
template <typename T>
void XPBD<T>::addJointAngular(size_t idA, size_t idB, const XPBD<T>::QuatT& targetQ, T alpha, const XPBD<T>::QuatT& aQ, const XPBD<T>::QuatT& bQ) {
  if(idA==idB)
    throw std::runtime_error("Cannot add joint to the same Shape<T>!");
  Constraint<T> c;
  c._isValid=true;
  c._alpha=alpha;
  c._type=ConstraintType::JointAngular;
  c._shapeIdA=(int)idA;
  c._shapeIdB=(int)idB;
  c._targetQ=targetQ;
  c._aQ=aQ;
  c._bQ=bQ;
  _jointAngulars.push_back(c);
}
template <typename T>
void XPBD<T>::addAnimation(int frameNum, QCIter angularB, QCIter angularE, QCIter rootQB, QCIter rootQE, XCIter rootXB, XCIter rootXE) {
  _animationFrameId=0;
  _animationFrameNum=frameNum;
  _isPlay=true;
  thrust::host_vector<QuatT> tempHostAngular(angularB, angularE);
  _animationData.resize(angularE-angularB);
  _animationData=tempHostAngular;
  thrust::host_vector<QuatT> tempHostRootQ(rootQB, rootQE);
  _rootAnimationQ.resize(frameNum);
  _rootAnimationQ=tempHostRootQ;
  thrust::host_vector<Vec3T> tempHostRootX(rootXB, rootXE);
  _rootAnimationX.resize(frameNum);
  _rootAnimationX=tempHostRootX;
}
template <typename T>
void XPBD<T>::playAnimation() {
  if(!_isPlay) return;
  int numJoints = _animationData.size()/_animationFrameNum;
  const auto b = _animationData.begin() + numJoints*_animationFrameId;
  updateJointAngular(b+1, b+numJoints);
  Vec3T* d_rootX = thrust::raw_pointer_cast(_rootAnimationX.data());
  QuatT* d_rootQ = thrust::raw_pointer_cast(_rootAnimationQ.data());
  Shape<T>* d_shapes = thrust::raw_pointer_cast(_geometry->getShapes());
  // TODO  Each SK has its own root, so there will be a list of roots and rootQ.
  int animationFrameId=_animationFrameId;
  thrust::for_each(thrust::device,
                   thrust::make_counting_iterator(0),
                   thrust::make_counting_iterator(1),
  [d_shapes, d_rootX, d_rootQ, animationFrameId] __device__ (int idx) {
    Shape<T>& root = d_shapes[idx];
    root._q = d_rootQ[animationFrameId];
    root._x = d_rootX[animationFrameId];
  });
  _animationFrameId=(_animationFrameId+1)%_animationFrameNum;
}
template <typename T>
void XPBD<T>::updateJointAngular(typename thrust::device_vector<QuatT>::const_iterator b, typename thrust::device_vector<QuatT>::const_iterator e) {
  assert(e-b == _jointAngulars.size());
  thrust::for_each(
    thrust::make_zip_iterator(thrust::make_tuple(_jointAngulars.begin(), b)),
    thrust::make_zip_iterator(thrust::make_tuple(_jointAngulars.end(), e)),
  [] __device__ (thrust::tuple<Constraint<T>&, const QuatT&> t) {
    auto& aJ = t.template get<0>();
    auto& q = t.template get<1>();
    aJ._targetQ=q;
  }
  );
}
template <typename T>
void XPBD<T>::addGroupLink(int a, int b) {
  _groupLinks.push_back(groupLink(a,b));
}
template <typename T>
void XPBD<T>::assignCollisionGroup() {
  if(_collisionGroupAssigned)
    return;
  //generate parent id
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
  //Step 2: assign collision parent id
  Shape<T>* d_shapes = thrust::raw_pointer_cast(_geometry->getShapes());
  thrust::for_each(thrust::device,
                   thrust::make_counting_iterator(0),
                   thrust::make_counting_iterator((int)_geometry->size()),
  [d_shapes, d_parents] __device__ (int idx) {
    d_shapes[idx]._parent = d_parents[idx];
  });
  //generate group id
  //Step 1: do union find in shapes.
  int* d_groupid = d_parents;
  auto& _groupid = _parents;
  thrust::sequence(thrust::device, _groupid.begin(), _groupid.end());
  do {
    thrust::transform(thrust::device, _groupLinks.begin(), _groupLinks.end(), _changes.begin(),
    [d_groupid] __device__ (const groupLink& c) {
      int groupidA = d_groupid[c._shapeIdA];
      int groupidB = d_groupid[c._shapeIdB];
      if (groupidA != groupidB) {
        int groupid = min(groupidA, groupidB);
        atomicMin(&d_groupid[c._shapeIdA], groupid);
        atomicMin(&d_groupid[c._shapeIdB], groupid);
        return true;
      }
      return false;
    });
    int change_count = thrust::count_if(thrust::device, _changes.begin(), _changes.end(), thrust::identity<bool>());
    changed = change_count > 0;
  } while (changed);
  //Step 2: assign collision parent id
  thrust::for_each(thrust::device,
                   thrust::make_counting_iterator(0),
                   thrust::make_counting_iterator((int)_geometry->size()),
  [d_shapes, d_groupid] __device__ (int idx) {
    d_shapes[idx]._groupid = d_groupid[idx];
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
