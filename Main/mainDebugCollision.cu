#include <PBD/Utils.h>
#include <PBD/Geometry.h>
#include <PBD/Collision.h>
#include <random>

using namespace GPUPBD;

int main() {
  typedef LSCALAR T;
  constexpr std::size_t N=10;
  std::vector<Capsule<T>> ps(N);

  std::mt19937 mt(123456789);
  std::uniform_real_distribution<T> uni(0.0, 1.0);

  for(auto& p : ps) {
    p._len = uni(mt);
    p._radius = uni(mt);
    Eigen::Quaternion<T> q;
    q.x() = uni(mt);
    q.y() = uni(mt);
    q.z() = uni(mt);
    q.w() = uni(mt);
    Eigen::Matrix<T, 3, 3> rot = q.toRotationMatrix();
    Eigen::Matrix<T, 3, 1> trans;
    trans.setRandom();
    p._trans << rot, trans;
  }

  std::shared_ptr<GPUPBD::Geometry<T>> geometry = std::make_shared<GPUPBD::Geometry<T>>();
  geometry->reserve(ps.size());
  geometry->resize(ps.size());
  geometry->setCapsule(ps);
  GPUPBD::CollisionDetector<T> detector(geometry);
  detector.detectCollisions();
  return 0;
}
