#include <PBD/Utils.h>
#include <PBD/Capsule.h>
#include <PBD/Geometry.h>
#include <PBD/Collision.h>
#include <random>
#include <fstream>
#include <iostream>


using namespace GPUPBD;

template <typename T>
std::ostream &operator<<(std::ostream &out, const Collision<T> &collision)
{
  out << collision._capsuleIdA << " "
      << collision._capsuleIdB << " ";
  for (int i = 0; i < collision._localPointA.rows(); ++i)
  {
    out << collision._localPointA(i) << " ";
  }
  for (int i = 0; i < collision._localPointB.rows(); ++i)
  {
    out << collision._localPointB(i) << " ";
  }
  for (int i = 0; i < collision._globalNormal.rows(); ++i)
  {
    out << collision._globalNormal(i) << " ";
  }
  out << collision._isValid;
  return out;
}

template <typename T>
std::istream &operator>>(std::istream &in, Collision<T> &collision)
{
  in >> collision._capsuleIdA >> collision._capsuleIdB;
  for (int i = 0; i < collision._localPointA.rows(); ++i)
  {
    in >> collision._localPointA(i);
  }
  for (int i = 0; i < collision._localPointB.rows(); ++i)
  {
    in >> collision._localPointB(i);
  }
  for (int i = 0; i < collision._globalNormal.rows(); ++i)
  {
    in >> collision._globalNormal(i);
  }
  in >> collision._isValid;
  return in;
}

template <typename T>
std::ostream &operator<<(std::ostream &out, const Capsule<T> &capsule)
{
  out << capsule._len << " " << capsule._radius << " ";
  for (int i = 0; i < capsule._trans.rows(); ++i)
  {
    for (int j = 0; j < capsule._trans.cols(); ++j)
    {
      out << capsule._trans(i, j) << " ";
    }
  }
  return out;
}

template <typename T>
std::istream &operator>>(std::istream &in, Capsule<T> &capsule)
{
  in >> capsule._len >> capsule._radius;
  for (int i = 0; i < capsule._trans.rows(); ++i)
  {
    for (int j = 0; j < capsule._trans.cols(); ++j)
    {
      in >> capsule._trans(i, j);
    }
  }
  return in;
}


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
