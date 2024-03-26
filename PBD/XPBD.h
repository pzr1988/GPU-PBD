#ifndef XPBD_H
#define XPBD_H

#include "Geometry.h"
#include "Collision.h"

namespace GPUPBD {
//The collisionDetector has the capability of detecting all pairs of collisions between all pairs of shapes
//The main function (detectCollisions()) is supposed to fill up the vector _collisions with active collision constraints
template <typename T>
class XPBD {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  typedef typename std::vector<QuatT>::const_iterator QCIter;
  typedef typename std::vector<Vec3T>::const_iterator XCIter;
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
  struct groupLink {
    int _shapeIdA;
    int _shapeIdB;
    groupLink(int a, int b):_shapeIdA(a),_shapeIdB(b) {}
  };
  //The contructor takes the set of shapes
  XPBD(std::shared_ptr<Geometry<T>> geometry,T dt,int nRelax=8);
  //Timestepping the physical system
  void step();
  //Initial time integration for updating Shape._transNext
  void integrate();
  //Relax all the constraints
  void relaxConstraint();
  //During relaxConstraints, the state's update will be cached.
  //Then we update state in this function.
  void updateShapeState();
  //We first compute Shape._v and Shape._w from Shape._trans(Next)
  //We then set Shape._trans=Shape._transNext;
  void updateVelocity();
  //Use detector to visulize collisions
  const CollisionDetector<T>& getDetector() const;
  //Get jointPositions for visulization
  const thrust::device_vector<Constraint<T>>& getJointPositions() const;
  //Get number of constraints
  size_t numConstraints() const;
  //Add joint between a pair of shapes
  void addJoint(size_t idA, size_t idB, const Vec3T& localA, const Vec3T& localB, T alpha=.0001f);
  //Fix angle between a pair of shapes, from B to A. aQ is the rotation of shape in its space. bQ is the rotation of shape in its space.
  void addJointAngular(size_t idA, size_t idB, const QuatT& targetQ, T alpha=.0001f, const QuatT& aQ=QuatT::Identity(), const QuatT& bQ=QuatT::Identity());
  //Add animation data
  void addAnimation(int frameNum, QCIter angularB, QCIter angularE, QCIter rootQB, QCIter rootQE);
  //Play Animation
  void playAnimation();
  //Update joint angular to play animation
  void updateJointAngular(typename thrust::device_vector<QuatT>::const_iterator b, typename thrust::device_vector<QuatT>::const_iterator e);
  //Assign a and b to the same group. Shapes in the same group will not collide.
  void addGroupLink(int a, int b);
  //Automatically avoid collision between a Shape<T> and its children
  void assignCollisionGroup();
  //Calculate globalnormal of joint position constraint, and axis/theta of joint angular constraint.
  void updateJointConstraint();
 private:
  //Init for relax all the constraints process
  void initRelaxConstraint();
  //Position constraint
  DEVICE_HOST static T computeGeneralizedInversMass(const Shape<T>& shape, const Vec3T& normal, const Vec3T& placementPoint);
  DEVICE_HOST static QuatT getDeltaRot(const Shape<T>& shape, const Vec3T& placementPoint, const Vec3T& pulse);
  //Angular constrant
  DEVICE_HOST static T computeGeneralizedInversMass(const Shape<T>& c, const Vec3T& n);
  DEVICE_HOST static QuatT getDeltaRot(const Shape<T>& c, const Vec3T& pulse);
  //data
  std::shared_ptr<Geometry<T>> _geometry;
  std::shared_ptr<CollisionDetector<T>> _detector;
  thrust::device_vector<Constraint<T>> _jointPositions;
  thrust::device_vector<Constraint<T>> _jointAngulars;
  thrust::device_vector<groupLink> _groupLinks;
  T _dt;
  int _nRelax;
  int _animationFrameId;
  int _animationFrameNum;
  bool _isPlay;
  thrust::device_vector<QuatT> _animationData;
  thrust::device_vector<QuatT> _rootAnimationQ;
  thrust::device_vector<Vec3T> _rootAnimationX;
  thrust::device_vector<T> _lambda;
  bool _collisionGroupAssigned;
  //Cache deltaX and deltaQ during relaxConstraint to avoid multi write problem
  //Jacobian aggregation
  thrust::device_vector<int> _constraintShapeId;
  thrust::device_vector<Update> _update;
  thrust::device_vector<int> _reduceShapeId;
  thrust::device_vector<Update> _reduceUpdate;
};

}
#endif
