#ifndef XPBD_H
#define XPBD_H

#include "Geometry.h"
#include "Collision.h"

namespace GPUPBD {
//The collisionDetector has the capability of detecting all pairs of collisions between all pairs of capsules
//The main function (detectCollisions()) is supposed to fill up the vector _collisions with active collision constraints
template <typename T>
class XPBD {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  struct Update {
    Vec3T _x;
    Vec4T _q;
  };
  struct UpdateAdd {
    DEVICE_HOST Update operator()(const Update& a, const Update& b) const {
      Update result;
      result._x = a._x + b._x;
      result._q = a._q + b._q;
      return result;
    }
  };
  //The contructor takes the set of capsules
  XPBD(std::shared_ptr<Geometry<T>> geometry,T dt,int nRelax=8);
  //Timestepping the physical system
  void step();
  //Initial time integration for updating Capsule._transNext
  void integrate();
  //Relax all the constraints
  void relaxConstraint();
  //During relaxConstraints, the state's update will be cached.
  //Then we update state in this function.
  void updateCapsuleState();
  //We first compute Capsule._v and Capsule._w from Capsule._trans(Next)
  //We then set Capsule._trans=Capsule._transNext;
  void updateVelocity();
  //Use detector to visulize collisions
  const CollisionDetector<T>& getDetector() const;
  //Get number of constraints
  size_t numConstraints() const;
  //Add joint between a pair of capsules
  void addJoint(size_t idA, size_t idB, const Vec3T& localA, const Vec3T& localB);
  //Automatically avoid collision between a Capsule<T> and its children
  void assignCollisionGroup();
 private:
  //Init for relax all the constraints process
  void initRelaxConstraint();
  DEVICE_HOST static T computeGeneralizedInversMass(const Capsule<T>& capsule, const Vec3T& normal, const Vec3T& placementPoint);
  DEVICE_HOST static QuatT getDeltaRot(const Capsule<T>& capsule, const Vec3T& placementPoint, const Vec3T& pulse);
  //data
  std::shared_ptr<Geometry<T>> _geometry;
  std::shared_ptr<CollisionDetector<T>> _detector;
  thrust::device_vector<Constraint<T>> _joints;
  T _dt;
  int _nRelax;
  thrust::device_vector<T> _lambda;
  bool _collisionGroupAssigned;
  //Cache deltaX and deltaQ during relaxConstraint to avoid multi write problem
  //Jacobian aggregation
  thrust::device_vector<int> _constraintCapsuleId;
  thrust::device_vector<Update> _update;
  thrust::device_vector<int> _reduceCapsuleId;
  thrust::device_vector<Update> _reduceUpdate;
};

}
#endif
