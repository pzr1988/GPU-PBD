#ifndef ANIMATION_PARSER_H
#define ANIMATION_PARSER_H

#include <PBD/Pragma.h>
#include <PBD/Geometry.h>
#include <fstream>

namespace PHYSICSMOTION {

template <typename T>
struct AnimationData {
  DECL_MAT_VEC_MAP_TYPES_T
  AnimationData(const std::string& animationFile, const std::string& rootTranslationFile);
  std::vector<QuatT> _animation;
  std::vector<QuatT> _rootQ;
  std::vector<Vec3T> _rootX;
  int _frameNum;
};

template <typename T>
AnimationData<T>::AnimationData(const std::string& animationFileName, const std::string& rootTranslationFileName) {
  std::ifstream animationFile(animationFileName);
  if (!animationFile) {
    std::cerr << "Unable to open file" << std::endl;
    return;
  }
  std::vector<std::vector<double>> rawData;
  double value;
  while (!animationFile.eof()) {
    std::vector<double> row;
    for (int i = 0; i < 4; ++i) {
      if (animationFile >> value) {
        row.push_back(value);
      }
    }
    if (!row.empty()) {
      rawData.push_back(row);
    }
  }
  animationFile.close();
  int numJoints = 20; // TODO fix hardcode
  _frameNum = rawData.size() / numJoints;
  _animation.resize(rawData.size());
  _rootQ.resize(_frameNum);
  // TODO fix hardcode of root local rotation.
  QuatT rootLocalQ = QuatT::FromTwoVectors(Vec3T::UnitX(),Vec3T(0.0237-0.0335, -0.0861-0.0849, -0.0278-(-0.0278)));
  for(int i=0; i<rawData.size(); i++) {
    const auto& v = rawData.at(i);
    QuatT q(v[0],v[1],v[2],v[3]);
    _animation[i]=q;
    if(i%numJoints==0) {
      _rootQ[i/numJoints]=q*rootLocalQ;
    }
  }
  std::ifstream rootTranslationFile(rootTranslationFileName);
  if (!rootTranslationFile) {
    std::cerr << "Unable to open file" << std::endl;
    return;
  }
  while (!rootTranslationFile.eof()) {
    std::vector<double> row;
    for (int i = 0; i < 3; ++i) {
      if (rootTranslationFile >> value) {
        row.push_back(value);
      }
    }
    if (!row.empty()) {
      // TODO fix hardcode of root local translation.
      _rootX.push_back(Vec3T(row[0],row[1],row[2])-Vec3T(3.1954-0.0208, -0.0175, 0.9931-1.0399));
    }
  }
  rootTranslationFile.close();
}

//declare instance
template struct AnimationData<LSCALAR>;
}
#endif
