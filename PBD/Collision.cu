#include "Collision.h"

namespace GPUPBD {
template <typename T>
CollisionDetector<T>::CollisionDetector(std::shared_ptr<Geometry<T>> geometry)
  :_geometry(geometry) {}
template <typename T>
void CollisionDetector<T>::detectCollisions() {
  //To see how to implement this function, refer to the tutorial:
  //http://media.steampowered.com/apps/valve/2015/DirkGregorius_Contacts.pdf
  //To see a sample code of how this can be implemented on CPU, see:
  //https://bitbucket.org/runningblade/libdifferentiable/src/ConvexHullPBAD/Environment/ContactGenerator.h
  FUNCTION_NOT_IMPLEMENTED
}
template <typename T>
Collision<T> CollisionDetector<T>::operator[](int id) {
  return _collisions[id];
}
template <typename T>
int CollisionDetector<T>::size() const {
  return _collisions.size();
}
}
