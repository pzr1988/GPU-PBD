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
  //Init for relax all the constraints process
  void initRelaxConstraint();
  //Relax all the constraints
  void relaxConstraint();
  //During relaxConstraints, the state's update will be cached.
  //Then we update state in this function.
  void updateVelocity();
  //Use detector to visulize collisions
  const CollisionDetector<T>& getDetector() const;
  DEVICE_HOST static T computeGeneralizedInversMass(const Capsule<T>& capsule, const Vec3T& normal, const Vec3T& placementPoint);
  DEVICE_HOST static Eigen::Quaternion<T> getDeltaRot(const Capsule<T>& capsule, const Vec3T& placementPoint, const Vec3T& pulse);
 protected:
  std::shared_ptr<Geometry<T>> _geometry;
  std::shared_ptr<CollisionDetector<T>> _detector;
  T _dt;
  int _nRelax;
 private:
 private:
  thrust::device_vector<T> _lambda;
  thrust::device_vector<int> _capsuleIds;
  thrust::device_vector<int> _labels;
  thrust::device_vector<int> _uniqueLabels;
  thrust::device_vector<bool> _changes;
  thrust::device_vector<int> _labelOfCollision;
};

}
#endif
