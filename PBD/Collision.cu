#include "Collision.h"

namespace GPUPBD {
template <typename T>
CollisionDetector<T>::CollisionDetector(std::shared_ptr<Geometry<T>> geometry)
  :_geometry(geometry) {}
template <typename T>
void CollisionDetector<T>::detectCollisions() {
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
