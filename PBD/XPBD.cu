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
    capsule._qPrev = capsule._q;
    if(capsule._isDynamic) {
      capsule._v += capsule._force*dt/capsule._mass;
      capsule._x += capsule._v*dt;
      capsule._R = capsule._q.toRotationMatrix();
      capsule._Iinv = capsule._R*capsule._Ibodyinv*capsule._R.transpose();
      capsule._w += capsule._Iinv*(capsule._torque
                                   - capsule._w.cross(capsule._R*capsule._Ibody*capsule._R.transpose()*capsule._w))*dt;
      Eigen::Quaternion<T> wQuat(0,capsule._w.x(),capsule._w.y(),capsule._w.z());
      Eigen::Quaternion<T> updatedQuat = Eigen::Quaternion<T>(0.5*dt,0,0,0)*wQuat*capsule._q;
      capsule._q = Eigen::Quaternion<T>(capsule._q.coeffs() + updatedQuat.coeffs());
      capsule._q.normalize();
    }
  });
  cudaDeviceSynchronize();
}
__device__ int binary_search(const int* data, int left, int right, int value) {
  while (left <= right) {
    int mid = left + (right - left) / 2;
    if (data[mid] == value) {
      return mid;
    } else if (data[mid] < value) {
      left = mid + 1;
    } else {
      right = mid - 1;
    }
  }
  return -1;
}
template <typename T>
struct UpdateLabels {
  int n;
  int* idMap;
  int* labels;
  UpdateLabels(int _n, int* _idMap, int* _labels) : n(_n), idMap(_idMap), labels(_labels) {}
  __device__ bool operator()(const Collision<T>& c) {
    int u = binary_search(idMap, 0, n-1, c._capsuleIdA);
    int v = binary_search(idMap, 0, n-1, c._capsuleIdB);
    int label_u = labels[u];
    int label_v = labels[v];

    if (label_u != label_v) {
      int min_label = min(label_u, label_v);
      atomicMin(&labels[u], min_label);
      atomicMin(&labels[v], min_label);
      return true;
    }
    return false;
  }
};
template <typename T>
void XPBD<T>::initRelaxConstraint() {
  if(_detector->size() == 0) {
    return;
  }
  _lambda.clear();
  int collisionNum = _detector->size();
  _lambda.resize(collisionNum);
  _changes.resize(collisionNum);
  _capsuleIds.resize(collisionNum*2);
  int* d_capsuleIds = thrust::raw_pointer_cast(_capsuleIds.data());
  const auto& collisions = _detector->getCollisions();
  const Collision<T>* d_collisions = thrust::raw_pointer_cast(collisions.data());
  thrust::for_each(thrust::device,
                   thrust::make_counting_iterator(0),
                   thrust::make_counting_iterator(collisionNum),
  [=] __host__ __device__ (int idx) {
    auto& collision = d_collisions[idx];
    d_capsuleIds[2*idx] = collision._capsuleIdA;
    d_capsuleIds[2*idx+1] = collision._capsuleIdB;
  });
  auto endUnique = thrust::unique(_capsuleIds.begin(), _capsuleIds.end());
  _capsuleIds.resize(thrust::distance(_capsuleIds.begin(), endUnique));
  int numCapsuleIds = _capsuleIds.size();
  thrust::sort(_capsuleIds.begin(), _capsuleIds.end());
  _labels.resize(numCapsuleIds);
  thrust::copy(_capsuleIds.begin(), _capsuleIds.end(), _labels.begin());
  int* d_labels = thrust::raw_pointer_cast(_labels.data());
  bool changed;
  do {
    thrust::transform(thrust::device, collisions.begin(), collisions.end(),
                      _changes.begin(), UpdateLabels<T>(numCapsuleIds, d_capsuleIds, d_labels));
    changed = thrust::reduce(thrust::device, _changes.begin(), _changes.end(), false, thrust::logical_or<bool>());
  } while (changed);
  _uniqueLabels.resize(_labels.size());
  auto endUniqueLabels = thrust::unique_copy(_labels.begin(), _labels.end(), _uniqueLabels.begin());
  _uniqueLabels.resize(thrust::distance(_uniqueLabels.begin(), endUniqueLabels));

}
template <typename T>
void XPBD<T>::relaxConstraint() {
  if(_detector->size() == 0) {
    return;
  }
  const auto& collisions = _detector->getCollisions();
  int numCollision = collisions.size();
  const Collision<T>* d_collisions = thrust::raw_pointer_cast(collisions.data());
  auto& capsules = _geometry->getMutableCapsules();
  Capsule<T>* d_capsules = thrust::raw_pointer_cast(capsules.data());
  T* d_lambda = thrust::raw_pointer_cast(_lambda.data());
  T dt = _dt;
  int numCapsuleIds = _capsuleIds.size();
  int* d_capsuleIds = thrust::raw_pointer_cast(_capsuleIds.data());
  int* d_labels = thrust::raw_pointer_cast(_labels.data());
  _labelOfCollision.resize(numCollision);
  int* d_labelOfCollision = thrust::raw_pointer_cast(_labelOfCollision.data());
  thrust::for_each(thrust::device,
                   thrust::make_counting_iterator<std::size_t>(0),
                   thrust::make_counting_iterator<std::size_t>(numCollision),
  [=] __device__ (int idx) {
    auto& collision = d_collisions[idx];
    d_labelOfCollision[idx] = d_labels[binary_search(d_capsuleIds, 0, numCapsuleIds, collision._capsuleIdA)];
  });
  thrust::for_each(thrust::device,
                   _uniqueLabels.begin(),
                   _uniqueLabels.end(),
  [=] __device__ (int targetLabel) {
    for(int idx=0; idx < numCollision; idx++) {
      int label = d_labels[binary_search(d_capsuleIds, 0, numCapsuleIds, d_collisions[idx]._capsuleIdA)];
      if(label == targetLabel) {
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
        if(collisionDepth<=0)
          return;
        auto alpha = collision._alpha/(dt*dt);
        auto deltaLambda = (-collisionDepth-d_lambda[idx]*alpha)/(wA+wB+alpha);
        d_lambda[idx] += deltaLambda;
        auto pulse = deltaLambda*collision._globalNormal;
        if(cA._isDynamic) {
          cA._x += pulse/cA._mass;
          auto cAQUpdated = getDeltaRot(cA, placementPointA, pulse);
          cA._q = Eigen::Quaternion<T>(cA._q.coeffs()+cAQUpdated.coeffs());
          cA._q.normalize();
        }
        if(cB._isDynamic) {
          cB._x -= pulse/cB._mass;
          auto cBQUpdated = getDeltaRot(cB, placementPointB, pulse);
          cB._q = Eigen::Quaternion<T>(cB._q.coeffs()-cBQUpdated.coeffs());
          cB._q.normalize();
        }
      }
    }
  });
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
DEVICE_HOST Eigen::Quaternion<T> XPBD<T>::getDeltaRot(const Capsule<T>& c, const Vec3T& r, const Vec3T& pulse) {
  auto cIinv = c.getInertiaTensorInv();
  auto cIinvRCrossP = cIinv * (r.cross(pulse)); // I^{-1}(r x p)
  Eigen::Quaternion<T> cIinvRCrossPQuat(0,cIinvRCrossP.x(),cIinvRCrossP.y(),cIinvRCrossP.z());
  auto qUpdated = Eigen::Quaternion<T>(0.5,0,0,0)*cIinvRCrossPQuat*c._q;
  return qUpdated;
}
//declare instance
template struct XPBD<LSCALAR>;
}
