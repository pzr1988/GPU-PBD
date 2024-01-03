#ifndef XPBD_H
#define XPBD_H
#include "Collision.h"
#include "PBD/Geometry.h"

namespace GPUPBD {
//The collisionDetector has the capability of detecting all pairs of collisions between all pairs of capsules
//The main function (detectCollisions()) is supposed to fill up the vector _collisions with active collision constraints
template <typename T>
class XPBD {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  //The contructor takes the set of capsules
  XPBD(std::shared_ptr<Geometry<T>> geometry,T dt,int nRelax=8);
  //Timestepping the physical system
  void step();
  //Initial time integration for updating Capsule._transNext
  void integrate();
  //Relax all the constraints
  void relaxConstraint();
  //We first compute Capsule._v and Capsule._w from Capsule._trans(Next)
  //We then set Capsule._trans=Capsule._transNext;
  void updateVelocity();
  //Use detector to visulize collisions
  const CollisionDetector<T>& getDetector() const;
 protected:
  std::shared_ptr<Geometry<T>> _geometry;
  std::shared_ptr<CollisionDetector<T>> _detector;
  T _dt;
  int _nRelax;
 private:
  //Init for relax all the constraints process
  void initRelaxConstraint();
  DEVICE_HOST static T computeGeneralizedInversMass(const Capsule<T>& capsule, const Vec3T& normal, const Vec3T& placementPoint);
  DEVICE_HOST static Eigen::Quaternion<T> getDeltaRot(const Capsule<T>& capsule, const Vec3T& placementPoint, const Vec3T& pulse);
 private:
  thrust::device_vector<T> _lambda;
  //Cache deltaX and deltaQ during relaxConstraint to avoid multi write problem
  thrust::device_vector<int> _collisionCapsuleId;
  thrust::device_vector<Vec3T> _deltaX;
  thrust::device_vector<Vec4T> _deltaQ;
  thrust::device_vector<int> _reduceCapsuleId;
  thrust::device_vector<Vec3T> _reduceDeltaX;
  thrust::device_vector<Vec4T> _reduceDeltaQ;
};

}
#endif
