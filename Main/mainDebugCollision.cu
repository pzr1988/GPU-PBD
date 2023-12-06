#include <PBD/Utils.h>
#include <PBD/Geometry.h>
#include <PBD/Collision.h>
#include <random>
#include <fstream>
#include <iostream>


using namespace GPUPBD;

int main() {
  typedef LSCALAR T;
  constexpr std::size_t N=10;
  std::vector<Capsule<T>> ps(N);

  std::mt19937 mt(123456789);
  std::uniform_real_distribution<T> uni(0.0, 1.0);

  for(auto& p : ps) {
    p._len = uni(mt);
    p._radius = uni(mt)/3.0;
    Eigen::Quaternion<T> q;
    q.x() = uni(mt);
    q.y() = uni(mt);
    q.z() = uni(mt);
    q.w() = uni(mt);
    q.normalize();
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

  // write file for TinyVisualizer
  std::ofstream file("collision.txt");
  // 格式：Capsule _radius _len _trans
  if (!file.is_open()) {
    std::cerr << "Unable to open file" << std::endl;
    return; // 或者使用其他错误处理方法
  }
  for(auto& p : ps) {
    file << "Capsule" << std::endl;
    file << p << std::endl;
  }

  size_t num_collision = detector.size();
  for(size_t i = 0; i < num_collision; i++) {
    auto c = detector[i];
    file << "Collision" << std::endl;
    file << c << std::endl;
  }

  file.close();


  return 0;
}
