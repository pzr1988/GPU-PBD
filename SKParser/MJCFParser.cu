#include <SKParser/MJCFParser.h>

using namespace GPUPBD;
namespace PHYSICSMOTION {
template <typename T>
MJCFParser<T>::MJCFParser(const std::string& file) {
  readMJCF(file);
  std::cout<<"==========================original info==========================" << std::endl;
  for(int i=0; i<_bodies.size(); i++) {
    Body& b = _bodies[i];
    std::string indent(b._depth * 2, ' ');
    std::cout << indent  << b._name <<": x:" << b._localX[0] << " " << b._localX[1] << " " << b._localX[2]
              << ", q:" << b._localQ.w() << " " << b._localQ.x() << " " << b._localQ.y() << " " << b._localQ.z()
              << ", radius:" << b._radius
              << ", fromto:" << b._ft[0] << " " << b._ft[1] << " "<< b._ft[2] << " "<< b._ft[3] << " "<< b._ft[4] << " "<< b._ft[5]
              << ", parent Name: " << _bodies[b._parent>-1?b._parent:0]._name
              << ", joint._isValid: " << b._pc._isValid
              << ", joint._cAName: " << _bodies[b._pc._cA>-1?b._pc._cA:0]._name
              << ", joint._cBName: " << _bodies[b._pc._cB>-1?b._pc._cB:0]._name
              << std::endl;
  }
  updateShape();
  std::cout<<"==========================updated info==========================" << std::endl;
  for(int i=0; i<_bodies.size(); i++) {
    Body& b = _bodies[i];
    std::string indent(b._depth * 2, ' ');
    auto tmpx = b._c._q.toRotationMatrix() * Vec3T(b._c._len, 0, 0);
    std::cout << indent  << b._name <<": x:" << b._globalX[0] << " " << b._globalX[1] << " " << b._globalX[2]
              << ", q:" << b._globalQ.w() << " " << b._globalQ.x() << " " << b._globalQ.y() << " " << b._globalQ.z()
              << ", shape radius: " << b._c._radius << ", len:" << b._c._len
              <<", x:" << b._c._x[0] << " " << b._c._x[1] << " " << b._c._x[2]
              << ", q:" << b._c._q.w() << " " << b._c._q.x() << " " << b._c._q.y() << " " << b._c._q.z()
              << ", parent Name: " << _bodies[b._parent>-1?b._parent:0]._name
              << std::endl;
  }
}
template <typename T>
void MJCFParser<T>::modifyInitPosByAnimation(AnimationData<T>& animationData) {

}
template <typename T>
void MJCFParser<T>::getShape(std::vector<Shape<T>>& ps) {
  ps.clear();
  for(auto& b : _bodies) {
    Shape<T> c = b._c;
    c._v.setZero();
    c._w.setZero();
    c._torque.setZero();
    c.initInertiaTensor();
    c._force.setZero();
    c._isDynamic = true;
    ps.push_back(c);
  }
}
template <typename T>
void MJCFParser<T>::getPositionConstraint(std::vector<PositionConstraint<T>>& pc) {
  std::cout<<"===================position constraint info===================" << std::endl;
  pc.clear();
  for(auto& b : _bodies) {
    auto& j = b._pc;
    if(j._isValid) {
      pc.push_back(j);
      // xpbd.addJoint(j._cA,j._cB,j._cAPos,j._cBPos);
      std::cout<<"Parent: " << _bodies[j._cB]._name << ", Self: " << _bodies[j._cA]._name
               << ", Parent Pos: " << j._cBPos[0] << ", "  << j._cBPos[2] << ", " << j._cBPos[2]
               << ", Self Pos: " << j._cAPos[0] << ", "  << j._cAPos[2] << ", " << j._cAPos[2]
               << std::endl;
    }
  }
}
template <typename T>
void MJCFParser<T>::getAngularConstraint(std::vector<AngularConstraint<T>>& ac) {
  std::cout<<"===================angular constraint info===================" << std::endl;
  for(auto& b :_bodies) {
    auto& j = b._ac;
    if(j._isValid) {
      ac.push_back(j);
      std::cout<<"Parent: " << _bodies[j._cB]._name << ", Self: " << _bodies[j._cA]._name
               << ", _pq: " << j._pQ.w() << ", "  << j._pQ.x() << ", " << j._pQ.y() << ", " << j._pQ.z()
               << ", _sq: " << j._sQ.w() << ", "  << j._sQ.x() << ", " << j._sQ.y() << ", " << j._sQ.z()
               << ", _psq: " << j._psQ.w() << ", "  << j._psQ.x() << ", " << j._psQ.y() << ", " << j._psQ.z()
               << std::endl;
    }
  }
}
template <typename T>
void MJCFParser<T>::readBodies(int parentId, const tinyxml2::XMLElement* g) {
  Body body;
  body._parent=parentId;
  body._name=g->Attribute("name");
  //compute depth
  body._depth=parentId==-1?0:_bodies[parentId]._depth+1;
  int _currId=(int)_bodies.size();
  //read trans
  {
    body._localX=parsePtreeDef<Vec3T>(*g,"<xmlattr>.pos","0 0 0");
    Vec4T tmpQ=parsePtreeDef<Vec4T>(*g,"<xmlattr>.quat","1 0 0 0");
    body._localQ=QuatT(tmpQ[0],tmpQ[1],tmpQ[2],tmpQ[3]);
  }
  //read geometry
  if(g->FirstChildElement("geom")->FindAttribute("type") == NULL) {
    //capsule
    body._isValid=true;
    body._type = ShapeType::Capsule;
    const tinyxml2::XMLElement* gg=g->FirstChildElement("geom");
    body._ft=parsePtreeDef<Vec6T>(*gg,"<xmlattr>.fromto","0 0 0 0 0 0");
    body._radius=get<T>(*gg,"<xmlattr>.size");
  } else if (std::string(g->FirstChildElement("geom")->FindAttribute("type")->Value()) == "box") {
    //box
    body._isValid=true;
    body._type = ShapeType::Box;
    const tinyxml2::XMLElement* gg=g->FirstChildElement("geom");
    body._boxPos=parsePtreeDef<Vec3T>(*gg,"<xmlattr>.pos","0 0 0");
    Vec4T tmpQ=parsePtreeDef<Vec4T>(*gg,"<xmlattr>.quat","1 0 0 0");
    body._boxQuat=QuatT(tmpQ[0],tmpQ[1],tmpQ[2],tmpQ[3]);
    body._boxSize=2*parsePtreeDef<Vec3T>(*gg,"<xmlattr>.size","0 0 0");
  } else if (std::string(g->FirstChildElement("geom")->FindAttribute("type")->Value()) == "sphere") {
    //sphere
    body._isValid=true;
    body._type = ShapeType::Sphere;
    const tinyxml2::XMLElement* gg=g->FirstChildElement("geom");
    body._spherePos=parsePtreeDef<Vec3T>(*gg,"<xmlattr>.pos","0 0 0");
    body._radius=get<T>(*gg,"<xmlattr>.size");
  } else {
    body._type = ShapeType::Unknown;
    body._isValid=false;
  }
  //read joints, temporarily ignore all angular constraints.
  if(g->FirstChildElement("joint")) {
    if(body._isValid) {
      body._pc._isValid=true;
      body._ac._isValid=true;
      body._ac._psQ=body._localQ;
      // parent:
      Body& p = _bodies[body._parent];
      body._pc._cB=body._ac._cB=body._parent;
      if(p._type == ShapeType::Capsule) {
        Vec3T pC1 = Vec3T(p._ft[0], p._ft[1], p._ft[2]);
        Vec3T pC2 = Vec3T(p._ft[3], p._ft[4], p._ft[5]);
        Vec3T pX = (pC1+pC2)/2;
        QuatT pQ = QuatT::FromTwoVectors(Vec3T::UnitX(),(pC2-pC1).normalized());
        body._pc._cBPos=pQ.inverse().toRotationMatrix()*(body._localX-pX);
        body._ac._pQ=pQ;
      } else if(p._type == ShapeType::Sphere) {
        body._pc._cBPos=p._localQ.inverse().toRotationMatrix()*(body._localX-p._spherePos);
        body._ac._pQ=QuatT::Identity();
      } else {
        body._pc._isValid=false;
        body._ac._isValid=false;
      }
      // son:
      body._pc._cA=body._ac._cA=_currId;
      if(body._type == ShapeType::Capsule) {
        Vec3T sC1 = Vec3T(body._ft[0], body._ft[1], body._ft[2]);
        Vec3T sC2 = Vec3T(body._ft[3], body._ft[4], body._ft[5]);
        Vec3T sX = (sC1+sC2)/2;
        QuatT sQ = QuatT::FromTwoVectors(Vec3T::UnitX(),(sC2-sC1).normalized());
        body._pc._cAPos=-sQ.inverse().toRotationMatrix()*(sX);
        body._ac._sQ=sQ;
      } else if (ShapeType::Box == body._type) {
        body._pc._cAPos=-body._boxQuat.inverse().toRotationMatrix()*(body._boxPos);
        body._ac._sQ=body._boxQuat;
      } else if(ShapeType::Sphere == body._type) {
        body._pc._cAPos=-body._spherePos;
        body._ac._sQ=QuatT::Identity();
      } else {
        body._pc._isValid=false;
        body._ac._isValid=false;
      }
    } else {
      body._pc._isValid=false;
      body._ac._isValid=false;
    }
  } else {
    body._pc._isValid=false;
    body._ac._isValid=false;
  }
  if(!body._isValid) return;
  _bodies.push_back(body);
  for(const tinyxml2::XMLElement* gc=g->FirstChildElement(); gc; gc=gc->NextSiblingElement())
    if(std::string(gc->Name()) == "body")
      readBodies(_currId,gc);
}
template <typename T>
void MJCFParser<T>::readMJCF(const std::string& file) {
  tinyxml2::XMLDocument pt;
  pt.LoadFile(file.c_str());
  tinyxml2::XMLElement* link=pt.RootElement();
  for(const tinyxml2::XMLElement* g=link->FirstChildElement(); g; g=g->NextSiblingElement())
    if(std::string(g->Name()) == "worldbody")
      readBodies(-1,g->FirstChildElement("body"));
}
template <typename T>
void MJCFParser<T>::updateShape() {
  for(auto& body : _bodies) {
    if(body._parent>=0) {
      const auto& p = _bodies[body._parent];
      body._globalX = p._globalX + p._globalQ.toRotationMatrix()*body._localX;
      body._globalQ = (p._globalQ * body._localQ).normalized();
    } else {
      body._globalX=body._localX;
      body._globalQ=body._localQ;
    }
    if(ShapeType::Capsule == body._type) {
      body._c._type = ShapeType::Capsule;
      Vec3T c1 = Vec3T(body._ft[0], body._ft[1], body._ft[2]);
      Vec3T c2 = Vec3T(body._ft[3], body._ft[4], body._ft[5]);
      body._c._radius=body._radius;
      body._c._len=(c2-c1).norm();
      body._c._x=body._globalQ.toRotationMatrix()*(c1+c2)/2+body._globalX;
      body._c._q=body._globalQ*QuatT::FromTwoVectors(Vec3T::UnitX(),(c2-c1).normalized());
      body._c._q.normalize();
    } else if (ShapeType::Box == body._type) {
      body._c._type = ShapeType::Box;
      body._c._x = body._globalX + body._globalQ.toRotationMatrix() * body._boxPos;
      body._c._q = (body._globalQ * body._boxQuat).normalized();
      body._c._radius = 0;
      body._c._len = body._boxSize[0];
      body._c._width = body._boxSize[1];
      body._c._height = body._boxSize[2];
    } else if (ShapeType::Sphere == body._type) {
      body._c._type = ShapeType::Sphere;
      body._c._x = body._globalX + body._globalQ.toRotationMatrix() * body._spherePos;
      body._c._q = body._globalQ;
      body._c._radius = body._radius;
    }
  }
}

//declare instance
template class PHYSICSMOTION::MJCFParser<LSCALAR>;
}