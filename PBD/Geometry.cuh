#ifndef GEOMETRY_CUH
#define GEOMETRY_CUH
#include "Utils.h"
#include <thrust/device_vector.h>
#include <LBVH/lbvh.cuh>

namespace GPUPBD {
//This struct represents a capsule-shaped object
//The capsule's centerline extends from (-_len/2,0,0) to (_len/2,0,0)
//The radius of the capsule is stored in _radius
//The local to global transformation is stored in _trans
//The 3x3 rotational matrix is: _trans.template block<3,3>(0,0) (you can also use macro: ROT(_trans))
//The 3x1 translational vector is: _trans.template block<3,1>(0,3) (you can also use macro: CTR(_trans))
template <typename T>
struct Capsule {
  DECL_MAT_VEC_MAP_TYPES_T
  T _len,_radius;
  Mat3X4T _trans;
};

// 获得胶囊体的bounding box
template<typename T>
struct aabb_getter
{
    __device__
    lbvh::aabb<float> operator()(const Capsule<T> c) const noexcept
    {
        lbvh::aabb<float> retval;
        // 假设胶囊体长轴与z轴对齐
        Eigen::Matrix<float, 4, 1> end1(0, 0, 0, 1); // 第一个端点
        Eigen::Matrix<float, 4, 1> end2(0, 0, static_cast<float>(c._len), 1); // 第二个端点

        Eigen::Matrix<float, 3, 1> transformedEnd1 = c._trans * end1;
        Eigen::Matrix<float, 3, 1> transformedEnd2 = c._trans * end2;
        
        Eigen::Matrix<float, 3, 1> upper = transformedEnd1.head<3>().cwiseMax(transformedEnd2.head<3>());
        float radius = static_cast<float>(c._radius);
        retval.upper.x = upper.x() + radius;
        retval.upper.y = upper.y() + radius;
        retval.upper.z = upper.z() + radius;
        Eigen::Matrix<float, 3, 1> lower = transformedEnd1.head<3>().cwiseMin(transformedEnd2.head<3>()) ;
        retval.lower.x = lower.x() - radius;
        retval.lower.y = lower.y() - radius;
        retval.lower.z = lower.z() - radius;
        return retval;
    }
};

//The geometry stores a vector of capsules
//The vector is over-sized and pre-allocated
//The number of capsules in use is stored in _nrCapsule
template <typename T>
class Geometry {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  //Set the actual number of capsules used
  void resize(int nrCapsule);
  //Set the pre-allocated capsule list
  void reserve(int nrCapsule);
  //CPU->GPU transfer: copying the list of capsules to GPU
  void setCapsule(const std::vector<Capsule<T>>& c);
  //Set the id-th capsule
  void setCapsule(int id,const Capsule<T>& c);
  //Get the id-th capsule
  Capsule<T> operator[](int id) const;
  //Get All Capsules
  const thrust::device_vector<Capsule<T>>& getCapsules() const;
 protected:
  thrust::device_vector<Capsule<T>> _capsules;
  int _nrCapsule=0;
};


template <typename T>
void Geometry<T>::resize(int nrCapsule) {
  _nrCapsule=nrCapsule;
  if(_nrCapsule<(int)_capsules.size())
    _capsules.resize(_nrCapsule);
}
template <typename T>
void Geometry<T>::reserve(int nrCapsule) {
  _capsules.resize(std::max<int>(nrCapsule,_nrCapsule));
}
template <typename T>
void Geometry<T>::setCapsule(const std::vector<Capsule<T>>& c) {
  // ASSERT_MSG((int)c.size()==_nrCapsule,"Incorrect number of capsules on host!")
  thrust::copy(c.begin(),c.end(),_capsules.begin());
}
template <typename T>
void Geometry<T>::setCapsule(int id,const Capsule<T>& c) {
  _capsules[id]=c;
}
template <typename T>
Capsule<T> Geometry<T>::operator[](int id) const {
  return _capsules[id];
}
template <typename T>
const thrust::device_vector<Capsule<T>>& Geometry<T>::getCapsules() const {
    return _capsules;
}
}

#endif