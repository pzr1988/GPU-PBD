#include "XPBD.h"

namespace GPUPBD {
template <typename T>
XPBD<T>::XPBD(std::shared_ptr<Geometry<T>> geometry,T dt,int nRelax)
  :_geometry(geometry),_detector(new CollisionDetector<T>(geometry)),_dt(dt),_nRelax(nRelax) {}
template <typename T>
void XPBD<T>::step() {
  _detector->detectCollisions();
  integrate();
  for(int i=0; i<_nRelax; i++)
    relaxConstraint();
  updateVelocity();
}
template <typename T>
void XPBD<T>::integrate() {
  //to be implemeneted
}
template <typename T>
void XPBD<T>::relaxConstraint() {
  //to be implemeneted
}
template <typename T>
void XPBD<T>::updateVelocity() {
  //to be implemeneted
}
template <typename T>
const CollisionDetector<T>& XPBD<T>::getDetector() const {
  if (!_detector) {
    throw std::runtime_error("Detector is not initialized");
  }
  return *_detector;
}
//declare instance
template struct XPBD<LSCALAR>;
}
