#ifndef ANIMATION_PARSER_H
#define ANIMATION_PARSER_H

#include <PBD/Pragma.h>
#include <PBD/Geometry.h>
#include <fstream>

namespace PHYSICSMOTION {

template <typename T>
struct AnimationData {
  DECL_MAT_VEC_MAP_TYPES_T
  AnimationData(const std::string& animationFile, const std::string& rootTranslationFile, const std::string& parentIndicesFile);
  void moveToMJCFRootPos(const QuatT& rootLocalQ, const Vec3T& rootLocalX);
  std::vector<int> _parentIndices;
  std::vector<QuatT> _animation;
  std::vector<QuatT> _rootQ;
  std::vector<Vec3T> _rootX;
  int _frameNum;
  int _numJoints;
};

template <typename T>
AnimationData<T>::AnimationData(const std::string& animationFileName, const std::string& rootTranslationFileName, const std::string& parentIndicesFileName) {
  std::ifstream parentIndicesFile(parentIndicesFileName);
  if (!parentIndicesFile) {
    std::cerr << "Unable to open file" << std::endl;
    return;
  }
  int indices;
  while (!parentIndicesFile.eof()) {
    parentIndicesFile >> indices;
    _parentIndices.push_back(indices);
  }
  parentIndicesFile.close();
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
  _numJoints = _parentIndices.size();
  _frameNum = rawData.size() / _numJoints;
  _animation.resize(rawData.size());
  _rootQ.resize(_frameNum);
  for(int i=0; i<rawData.size(); i++) {
    const auto& v = rawData.at(i);
    _animation[i]=QuatT(v[0],v[1],v[2],v[3]);
    if(i%_numJoints==0) {
      _rootQ[i/_numJoints]=_animation[i];
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
      _rootX.push_back(Vec3T(row[0],row[1],row[2]));
    }
  }
  rootTranslationFile.close();
}

template <typename T>
void AnimationData<T>::moveToMJCFRootPos(const QuatT& rootLocalQ, const Vec3T& rootLocalX) {
  for(int i=0; i<_rootQ.size(); i++)
    _rootQ[i]=_rootQ[i]*rootLocalQ;
  Vec3T shift = _rootX[0]-rootLocalX;
  for(int i=0; i<_rootX.size(); i++) {
    _rootX[i] = _rootX[i]-shift;
  }
}
//declare instance
template struct AnimationData<LSCALAR>;
}
#endif
