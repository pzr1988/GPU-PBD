#ifndef MJCFPARSER_H
#define MJCFPARSER_H

#include <PBD/Pragma.h>
#include <PBD/Geometry.h>
#include <SKParser/tinyxml2.h>
#include <SKParser/Utils.h>
#include <SKParser/AnimationParser.h>
#include <fstream>

namespace PHYSICSMOTION {

template <typename T>
struct PositionConstraint {
  DECL_MAT_VEC_MAP_TYPES_T
  bool _isValid=false;
  int _cA=-1;//son
  int _cB=-1;//parent
  Vec3T _cAPos;
  Vec3T _cBPos;
};
template <typename T>
struct AngularConstraint {
  DECL_MAT_VEC_MAP_TYPES_T
  bool _isValid=false;
  int _cA=-1;//son
  int _cB=-1;//parent
  QuatT _pQ; //the rotation of parent shape in it's space.
  QuatT _sQ; //the rotation of son shape in it's space.
  QuatT _psQ; //the rotation from parent's space to son's space.
};

template <typename T>
class MJCFParser {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  MJCFParser(const std::string& file);
  void modifyInitGestureByAnimation(const AnimationData<T>& animationData);
  void getRootPos(QuatT& rootLocalQ, Vec3T& rootLocalX);
  void getShape(std::vector<GPUPBD::Shape<T>>& ps);
  void getPositionConstraint(std::vector<PositionConstraint<T>>& pc);
  void getAngularConstraint(std::vector<AngularConstraint<T>>& ac);
 private:
  struct Body {
    bool _isValid=false;
    int _parent;
    int _depth;
    T _radius;
    Vec3T _localX;
    QuatT _localQ;
    Vec3T _globalX;
    QuatT _globalQ;
    Vec6T _ft;
    GPUPBD::ShapeType _type;
    Vec3T _boxPos;
    QuatT _boxQuat;
    Vec3T _boxSize;
    Vec3T _spherePos;
    std::string _name;
    GPUPBD::Shape<T> _c;
    PositionConstraint<T> _pc;
    AngularConstraint<T> _ac;
  };
  void readMJCF(const std::string& file);
  void readBodies(int parentId, const tinyxml2::XMLElement* g);
  void updateShape();
  std::vector<Body> _bodies;
};
}
#endif
