#include <PBD/Geometry.h>
#include <PBD/Collision.h>
#include <PBD/XPBD.h>
#include <PBD/Visualizer.h>
#include <TinyVisualizer/FirstPersonCameraManipulator.h>
#include <TinyVisualizer/CameraExportPlugin.h>
#include <TinyVisualizer/CaptureGIFPlugin.h>
#include <TinyVisualizer/ImGuiPlugin.h>
#include <TinyVisualizer/Camera3D.h>
#include <SKParser/tinyxml2.h>
#include <SKParser/Utils.h>
#include <fstream>

using namespace GPUPBD;

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
struct Body {
  DECL_MAT_VEC_MAP_TYPES_T
  bool _isValid=false;
  int _parent;
  int _depth;
  T _radius;
  Vec3T _x;
  QuatT _q;
  Vec6T _ft;
  ShapeType _type;
  Vec3T _boxPos;
  QuatT _boxQuat;
  Vec3T _boxSize;
  Vec3T _spherePos;
  std::string _name;
  Shape<T> _c;
  PositionConstraint<T> _pc;
  AngularConstraint<T> _ac;
};
template <typename T>
void readBodies(std::vector<Body<T>>& bodies, int parentId, const tinyxml2::XMLElement* g) {
  DECL_MAT_VEC_MAP_TYPES_T
  Body<T> body;
  body._parent=parentId;
  body._name=g->Attribute("name");
  //compute depth
  body._depth=parentId==-1?0:bodies[parentId]._depth+1;
  int _currId=(int)bodies.size();
  //read trans
  {
    body._x=PHYSICSMOTION::parsePtreeDef<Vec3T>(*g,"<xmlattr>.pos","0 0 0");
    Vec4T tmpQ=PHYSICSMOTION::parsePtreeDef<Vec4T>(*g,"<xmlattr>.quat","1 0 0 0");
    body._q=QuatT(tmpQ[0],tmpQ[1],tmpQ[2],tmpQ[3]);
  }
  //read geometry
  if(g->FirstChildElement("geom")->FindAttribute("type") == NULL) {
    //capsule
    body._isValid=true;
    body._type = ShapeType::Capsule;
    const tinyxml2::XMLElement* gg=g->FirstChildElement("geom");
    body._ft=PHYSICSMOTION::parsePtreeDef<Vec6T>(*gg,"<xmlattr>.fromto","0 0 0 0 0 0");
    body._radius=PHYSICSMOTION::get<T>(*gg,"<xmlattr>.size");
  } else if (std::string(g->FirstChildElement("geom")->FindAttribute("type")->Value()) == "box") {
    //box
    body._isValid=true;
    body._type = ShapeType::Box;
    const tinyxml2::XMLElement* gg=g->FirstChildElement("geom");
    body._boxPos=PHYSICSMOTION::parsePtreeDef<Vec3T>(*gg,"<xmlattr>.pos","0 0 0");
    Vec4T tmpQ=PHYSICSMOTION::parsePtreeDef<Vec4T>(*gg,"<xmlattr>.quat","1 0 0 0");
    body._boxQuat=QuatT(tmpQ[0],tmpQ[1],tmpQ[2],tmpQ[3]);
    body._boxSize=2*PHYSICSMOTION::parsePtreeDef<Vec3T>(*gg,"<xmlattr>.size","0 0 0");
  } else if (std::string(g->FirstChildElement("geom")->FindAttribute("type")->Value()) == "sphere") {
    //sphere
    body._isValid=true;
    body._type = ShapeType::Sphere;
    const tinyxml2::XMLElement* gg=g->FirstChildElement("geom");
    body._spherePos=PHYSICSMOTION::parsePtreeDef<Vec3T>(*gg,"<xmlattr>.pos","0 0 0");
    body._radius=PHYSICSMOTION::get<T>(*gg,"<xmlattr>.size");
  } else {
    body._type = ShapeType::Unknown;
    body._isValid=false;
  }

  //read joints, temporarily ignore all angular constraints.
  if(g->FirstChildElement("joint")) {
    if(body._isValid) {
      body._pc._isValid=true;
      body._ac._isValid=true;
      body._ac._psQ=body._q;
      // parent:
      Body<T>& p = bodies[body._parent];
      body._pc._cB=body._ac._cB=body._parent;
      if(p._type == ShapeType::Capsule) {
        Vec3T pC1 = Vec3T(p._ft[0], p._ft[1], p._ft[2]);
        Vec3T pC2 = Vec3T(p._ft[3], p._ft[4], p._ft[5]);
        Vec3T pX = (pC1+pC2)/2;
        QuatT pQ = QuatT::FromTwoVectors(Vec3T::UnitX(),(pC2-pC1).normalized());
        body._pc._cBPos=pQ.inverse().toRotationMatrix()*(body._x-pX);
        body._ac._pQ=pQ;
      } else if(p._type == ShapeType::Sphere) {
        body._pc._cBPos=p._q.inverse().toRotationMatrix()*(body._x-p._spherePos);
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

  bodies.push_back(body);
  for(const tinyxml2::XMLElement* gc=g->FirstChildElement(); gc; gc=gc->NextSiblingElement())
    if(std::string(gc->Name()) == "body")
      readBodies(bodies,_currId,gc);
}

template <typename T>
void readMJCF(std::vector<Body<T>>& bodies, const std::string& file) {
  tinyxml2::XMLDocument pt;
  pt.LoadFile(file.c_str());
  tinyxml2::XMLElement* link=pt.RootElement();
  for(const tinyxml2::XMLElement* g=link->FirstChildElement(); g; g=g->NextSiblingElement())
    if(std::string(g->Name()) == "worldbody")
      readBodies(bodies,-1,g->FirstChildElement("body"));
}

template <typename T>
void updateShape(std::vector<Body<T>>& bodies) {
  DECL_MAT_VEC_MAP_TYPES_T
  for(auto& body : bodies) {
    if(ShapeType::Capsule == body._type) {
      body._c._type = ShapeType::Capsule;
      Vec3T c1 = Vec3T(body._ft[0], body._ft[1], body._ft[2]);
      Vec3T c2 = Vec3T(body._ft[3], body._ft[4], body._ft[5]);
      body._c._q = QuatT::FromTwoVectors(Vec3T::UnitX(),(c2-c1).normalized());
      c1 = body._x + body._q.toRotationMatrix() * c1;
      c2 = body._x + body._q.toRotationMatrix() * c2;
      body._c._q = body._q * body._c._q;
      body._c._q.normalize();
      if(body._parent>=0) {
        const auto& p = bodies[body._parent];
        c1 = p._x + p._q.toRotationMatrix()*c1;
        c2 = p._x + p._q.toRotationMatrix()*c2;
        body._c._q = p._q * body._c._q;
        body._c._q.normalize();
        body._x = p._x + p._q.toRotationMatrix()*body._x;
        body._q = (p._q * body._q).normalized();
      }
      body._c._radius=body._radius;
      body._c._len=(c2-c1).norm();
      body._c._x=(c1+c2)/2;
    } else if (ShapeType::Box == body._type) {
      body._c._type = ShapeType::Box;
      Vec3T x = body._x + body._q.toRotationMatrix() * body._boxPos;
      QuatT q = (body._q * body._boxQuat).normalized();
      if(body._parent>=0) {
        const auto& p = bodies[body._parent];
        x = p._x + p._q.toRotationMatrix()*x;
        q = (p._q * q).normalized();
        body._x = p._x + p._q.toRotationMatrix()*body._x;
        body._q = (p._q * body._q).normalized();
      }
      body._c._x = x;
      body._c._q = q;
      body._c._radius = 0;
      body._c._len = body._boxSize[0];
      body._c._width = body._boxSize[1];
      body._c._height = body._boxSize[2];
    } else if (ShapeType::Sphere == body._type) {
      body._c._type = ShapeType::Sphere;
      Vec3T x = body._x + body._q.toRotationMatrix() * body._spherePos;
      QuatT q = body._q;
      if(body._parent>=0) {
        const auto& p = bodies[body._parent];
        x = p._x + p._q.toRotationMatrix()*x;
        q = (p._q * q).normalized();
        body._x = p._x + p._q.toRotationMatrix()*body._x;
        body._q = (p._q * body._q).normalized();
      }
      body._c._x = x;
      body._c._q = q;
      body._c._radius = body._radius;
    }
  }
}


int main(int argc,char** argv) {
  typedef LSCALAR T;
  DECL_MAT_VEC_MAP_TYPES_T
  std::vector<Body<T>> bodies;
  readMJCF(bodies, "/data/GPU-PBD/SKParser/MarathonCharacter_PhysicsAsset2.xml");
  std::cout<<"==========================original info==========================" << std::endl;
  for(int i=0; i<bodies.size(); i++) {
    Body<T>& b = bodies[i];
    std::string indent(b._depth * 2, ' ');
    std::cout << indent  << b._name <<": x:" << b._x[0] << " " << b._x[1] << " " << b._x[2]
              << ", q:" << b._q.w() << " " << b._q.x() << " " << b._q.y() << " " << b._q.z()
              << ", radius:" << b._radius
              << ", fromto:" << b._ft[0] << " " << b._ft[1] << " "<< b._ft[2] << " "<< b._ft[3] << " "<< b._ft[4] << " "<< b._ft[5]
              << ", parent Name: " << bodies[b._parent>-1?b._parent:0]._name
              << ", joint._isValid: " << b._pc._isValid
              << ", joint._cAName: " << bodies[b._pc._cA>-1?b._pc._cA:0]._name
              << ", joint._cBName: " << bodies[b._pc._cB>-1?b._pc._cB:0]._name
              << std::endl;
  }

  updateShape(bodies);
  std::cout<<"==========================updated info==========================" << std::endl;
  for(int i=0; i<bodies.size(); i++) {
    Body<T>& b = bodies[i];
    std::string indent(b._depth * 2, ' ');
    auto tmpx = b._c._q.toRotationMatrix() * Vec3T(b._c._len, 0, 0);
    std::cout << indent  << b._name <<": x:" << b._x[0] << " " << b._x[1] << " " << b._x[2]
              << ", q:" << b._q.w() << " " << b._q.x() << " " << b._q.y() << " " << b._q.z()
              << ", shape radius: " << b._c._radius << ", len:" << b._c._len
              <<", x:" << b._c._x[0] << " " << b._c._x[1] << " " << b._c._x[2]
              << ", q:" << b._c._q.w() << " " << b._c._q.x() << " " << b._c._q.y() << " " << b._c._q.z()
              << ", parent Name: " << bodies[b._parent>-1?b._parent:0]._name
              << std::endl;
  }

  std::vector<Shape<T>> ps;
  for(auto& b : bodies) {
    Shape<T> c = b._c;
    c._v.setZero();
    c._w.setZero();
    c._torque.setZero();
    c.initInertiaTensor();
    // c._force = Vec3T(0, -9.8f*c._mass,0);
    c._force.setZero();
    c._isDynamic = true;
    ps.push_back(c);
  }

  // floor
  Shape<T> b_1;
  b_1._type = ShapeType::Box;
  b_1._len=10;
  b_1._width=10;
  b_1._height=1;
  b_1._x = Vec3T(0,0,-0.5);
  b_1._q = QuatT(1, 0, 0, 0);
  b_1.initInertiaTensor();
  b_1._isDynamic = false;
  ps.push_back(b_1);

  std::shared_ptr<Geometry<T>> geometry(new Geometry<T>);
  geometry->resize(ps.size());
  geometry->setShape(ps);
  XPBD<T> xpbd(geometry, 1.0f/60);
  // addPositionConstraint
  std::cout<<"===================position constraint info===================" << std::endl;
  for(auto& b : bodies) {
    auto& j = b._pc;
    if(j._isValid) {
      xpbd.addJoint(j._cA,j._cB,j._cAPos,j._cBPos);
      std::cout<<"Parent: " << bodies[j._cB]._name << ", Self: " << bodies[j._cA]._name
               << ", Parent Pos: " << j._cBPos[0] << ", "  << j._cBPos[2] << ", " << j._cBPos[2]
               << ", Self Pos: " << j._cAPos[0] << ", "  << j._cAPos[2] << ", " << j._cAPos[2]
               << std::endl;
    }
  }
  std::cout<<"===================angular constraint info===================" << std::endl;
  for(auto& b : bodies) {
    auto& j = b._ac;
    if(j._isValid) {
      xpbd.addJointAngular(j._cA,j._cB,j._psQ, 0.0001f, j._sQ, j._pQ);
      std::cout<<"Parent: " << bodies[j._cB]._name << ", Self: " << bodies[j._cA]._name
               << ", _pq: " << j._pQ.w() << ", "  << j._pQ.x() << ", " << j._pQ.y() << ", " << j._pQ.z()
               << ", _sq: " << j._sQ.w() << ", "  << j._sQ.x() << ", " << j._sQ.y() << ", " << j._sQ.z()
               << ", _psq: " << j._psQ.w() << ", "  << j._psQ.x() << ", " << j._psQ.y() << ", " << j._psQ.z()
               << std::endl;
    }
  }
  std::cout<<"===================add group info===================" << std::endl;
  std::vector<std::pair<int, int>> groupLinks= {{0,7},{7,8},{8,9},{9,10},{10,11},{9,14},{14,15},{9,16},{16,17},{1,4}};
  for(const auto& g : groupLinks) {
    xpbd.addGroupLink(g.first, g.second);
    std::cout << "Add " << bodies[g.first]._name << " and " << bodies[g.second]._name << " into the same group." << std::endl;
  }

  std::cout<<"==========================animation info==========================" << std::endl;
  std::ifstream inFile("/data/GPU-PBD/SKParser/animation.data");
  if (!inFile) {
    std::cerr << "Unable to open file" << std::endl;
    return 1;
  }
  std::vector<std::vector<double>> rawData;
  std::vector<QuatT> animation;
  double value;
  while (!inFile.eof()) {
    std::vector<double> row;
    for (int i = 0; i < 4; ++i) {
      if (inFile >> value) {
        row.push_back(value);
      }
    }
    if (!row.empty()) {
      rawData.push_back(row);
    }
  }
  inFile.close();
  int numJoints = 20;
  int framNum = rawData.size() / numJoints;
  animation.resize(rawData.size());
  for(int i=0; i<rawData.size(); i++) {
    const auto& v = rawData.at(i);
    QuatT q(v[0],v[1],v[2],v[3]);
    animation[i]=q;
  }
  thrust::device_vector<QuatT> d_animation(animation.begin(), animation.end());

  DRAWER::Drawer drawer(argc,argv);
  drawer.addPlugin(std::shared_ptr<DRAWER::Plugin>(new DRAWER::CameraExportPlugin(GLFW_KEY_2,GLFW_KEY_3,"camera.dat")));
  drawer.addPlugin(std::shared_ptr<DRAWER::Plugin>(new DRAWER::CaptureGIFPlugin(GLFW_KEY_1,"record.gif",drawer.FPS())));
  auto shapeGeometry=visualizeOrUpdateGeometry(*geometry);
  auto shapeCollision=visualizeOrUpdateCollision(*geometry,xpbd.getDetector(),xpbd.getJointPositions());
  drawer.addShape(shapeGeometry);
  drawer.addShape(shapeCollision);
  // drawer.addCamera3D(90,Eigen::Matrix<GLfloat,3,1>(0,1,0),Eigen::Matrix<GLfloat,3,1>(0,0,5),Eigen::Matrix<GLfloat,3,1>(0,0,-1));
  // drawer.getCamera3D()->setManipulator(std::shared_ptr<DRAWER::CameraManipulator>(new DRAWER::FirstPersonCameraManipulator(drawer.getCamera3D())));
  // drawer.addPlugin(std::shared_ptr<DRAWER::Plugin>(new DRAWER::ImGuiPlugin([&]() {
  //   drawer.getCamera3D()->getManipulator()->imGuiCallback();
  // })));
  bool sim=false;
  int frameId = 0;
  drawer.setFrameFunc([&](std::shared_ptr<DRAWER::SceneNode>& root) {
    if(sim) {
      xpbd.step();
      visualizeOrUpdateGeometry(*geometry,shapeGeometry);
      visualizeOrUpdateCollision(*geometry,xpbd.getDetector(),xpbd.getJointPositions(),shapeCollision);
      const auto b = d_animation.begin() + numJoints*frameId;
      xpbd.updateJointAngular(b+1, b+numJoints);
      ++frameId;
      frameId=frameId%framNum;
    }
  });
  //press R to run simulation
  drawer.setKeyFunc([&](GLFWwindow* wnd,int key,int scan,int action,int mods,bool captured) {
    if(captured)
      return;
    else if(key==GLFW_KEY_R && action==GLFW_PRESS)
      sim=!sim;
  });
  drawer.mainLoop();

  return 0;
}
