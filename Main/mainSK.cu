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

using namespace GPUPBD;

template <typename T>
struct Joint {
  DECL_MAT_VEC_MAP_TYPES_T
  bool _isValid=false;
  int _cA=-1;//son
  int _cB=-1;//parent
  Vec3T _cAPos;
  Vec3T _cBPos;
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
  std::string _name;
  Shape<T> _c;
  Joint<T> _j;
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
    body._boxQuat=PHYSICSMOTION::parsePtreeDef<Vec4T>(*gg,"<xmlattr>.quat","1 0 0 0");
    body._boxSize=PHYSICSMOTION::parsePtreeDef<Vec3T>(*gg,"<xmlattr>.size","0 0 0");
    printf("boxPos:%f,%f,%f, boxQuat:%f,%f,%f,%f, boxSize:%f,%f,%f\n",
           body._boxPos[0], body._boxPos[1],body._boxPos[2],
           body._boxQuat.w(), body._boxQuat.x(), body._boxQuat.y(), body._boxQuat.z(),
           body._boxSize[0], body._boxSize[1], body._boxSize[2]);
  } else {
    body._type = ShapeType::Unknown;
    body._isValid=false;
  }

  //read joints, temporarily ignore all angular constraints.
  if(g->FirstChildElement("joint")) {
    if(body._isValid && ShapeType::Capsule == body._type) {
      body._j._isValid=true;
      // parent:
      Body<T>& p = bodies[body._parent];
      Vec3T pC1 = Vec3T(p._ft[0], p._ft[1], p._ft[2]);
      Vec3T pC2 = Vec3T(p._ft[3], p._ft[4], p._ft[5]);
      Vec3T pX = (pC1+pC2)/2;
      QuatT pQ = QuatT::FromTwoVectors(Vec3T::UnitX(),(pC2-pC1).normalized());
      body._j._cB=body._parent;
      body._j._cBPos=pQ.inverse().toRotationMatrix()*(body._x-pX);
      // son:
      body._j._cA=_currId;
      Vec3T sC1 = Vec3T(body._ft[0], body._ft[1], body._ft[2]);
      Vec3T sC2 = Vec3T(body._ft[3], body._ft[4], body._ft[5]);
      Vec3T sX = (sC1+sC2)/2;
      QuatT sQ = QuatT::FromTwoVectors(Vec3T::UnitX(),(sC2-sC1).normalized());
      body._j._cAPos=-sQ.inverse().toRotationMatrix()*(sX);
    } else if (body._isValid && ShapeType::Box == body._type) {
      body._j._isValid=true;
      // parent is capsule
      Body<T>& p = bodies[body._parent];
      Vec3T pC1 = Vec3T(p._ft[0], p._ft[1], p._ft[2]);
      Vec3T pC2 = Vec3T(p._ft[3], p._ft[4], p._ft[5]);
      Vec3T pX = (pC1+pC2)/2;
      QuatT pQ = QuatT::FromTwoVectors(Vec3T::UnitX(),(pC2-pC1).normalized());
      body._j._cB=body._parent;
      body._j._cBPos=pQ.inverse().toRotationMatrix()*(body._x-pX);
      // son:
      body._j._cA=_currId;
      body._j._cAPos=-body._boxQuat.inverse().toRotationMatrix()*(body._boxPos);
    } else {
      body._j._isValid=false;
    }
  } else {
    body._j._isValid=false;
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
      c1 = body._x + body._q.toRotationMatrix() * c1;
      c2 = body._x + body._q.toRotationMatrix() * c2;
      if(body._parent>=0) {
        const auto& p = bodies[body._parent];
        c1 = p._x + p._q.toRotationMatrix()*c1;
        c2 = p._x + p._q.toRotationMatrix()*c2;
        body._x = p._x + p._q.toRotationMatrix()*body._x;
        body._q = (p._q * body._q).normalized();
      }
      body._c._radius=body._radius;
      body._c._len=(c2-c1).norm();
      body._c._x=(c1+c2)/2;
      body._c._q=QuatT::FromTwoVectors(Vec3T::UnitX(),(c2-c1).normalized());
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
      body._c._len = body._boxSize[0]*2;
      body._c._width = body._boxSize[1]*2;
      body._c._height = body._boxSize[2]*2;
    }
  }
}


int main(int argc,char** argv) {
  typedef LSCALAR T;
  DECL_MAT_VEC_MAP_TYPES_T
  std::vector<Body<T>> bodies;
  readMJCF(bodies, "/data/GPU-PBD/SKParser/SK_Mannequin_PhysicsAsset_ABFB4_MJCF.xml");
  std::cout<<"==========================original info==========================" << std::endl;
  for(int i=0; i<bodies.size(); i++) {
    Body<T>& b = bodies[i];
    std::string indent(b._depth * 2, ' ');
    std::cout << indent  << b._name <<": x:" << b._x[0] << " " << b._x[1] << " " << b._x[2]
              << ", q:" << b._q.w() << " " << b._q.x() << " " << b._q.y() << " " << b._q.z()
              << ", radius:" << b._radius
              << ", fromto:" << b._ft[0] << " " << b._ft[1] << " "<< b._ft[2] << " "<< b._ft[3] << " "<< b._ft[4] << " "<< b._ft[5]
              << ", parent Name: " << bodies[b._parent>-1?b._parent:0]._name
              << ", joint._isValid: " << b._j._isValid
              << ", joint._cAName: " << bodies[b._j._cA>-1?b._j._cA:0]._name
              << ", joint._cBName: " << bodies[b._j._cB>-1?b._j._cB:0]._name
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
    c._force = Vec3T(0, -9.8f*c._mass,0);
    c._isDynamic = true;
    ps.push_back(c);
  }

  // boundary
  Shape<T> b_1;
  b_1._type = ShapeType::Box;
  b_1._len=10;
  b_1._width=1;
  b_1._height=10;
  b_1._x = Vec3T(0,-4,0);
  b_1._q = QuatT(1,0,0,0);
  b_1.initInertiaTensor();
  b_1._isDynamic = false;
  ps.push_back(b_1);

  std::shared_ptr<Geometry<T>> geometry(new Geometry<T>);
  geometry->resize(ps.size());
  geometry->setShape(ps);
  XPBD<T> xpbd(geometry, 1.0f/60);
  // addJoint
  std::cout<<"==========================joint info==========================" << std::endl;
  for(auto& b : bodies) {
    auto& j = b._j;
    if(j._isValid) {
      xpbd.addJoint(j._cA,j._cB,j._cAPos,j._cBPos);
      std::cout<<"Parent: " << bodies[j._cB]._name << ", Self: " << bodies[j._cA]._name
               << ", Parent Pos: " << j._cBPos[0] << ", "  << j._cBPos[2] << ", " << j._cBPos[2]
               << ", Self Pos: " << j._cAPos[0] << ", "  << j._cAPos[2] << ", " << j._cAPos[2]
               << std::endl;
    }
  }
  DRAWER::Drawer drawer(argc,argv);
  drawer.addPlugin(std::shared_ptr<DRAWER::Plugin>(new DRAWER::CameraExportPlugin(GLFW_KEY_2,GLFW_KEY_3,"camera.dat")));
  drawer.addPlugin(std::shared_ptr<DRAWER::Plugin>(new DRAWER::CaptureGIFPlugin(GLFW_KEY_1,"record.gif",drawer.FPS())));
  auto shapeGeometry=visualizeOrUpdateGeometry(*geometry);
  auto shapeCollision=visualizeOrUpdateCollision(*geometry,xpbd.getDetector());
  drawer.addShape(shapeGeometry);
  drawer.addShape(shapeCollision);
  // drawer.addCamera3D(90,Eigen::Matrix<GLfloat,3,1>(0,1,0),Eigen::Matrix<GLfloat,3,1>(0,0,5),Eigen::Matrix<GLfloat,3,1>(0,0,-1));
  // drawer.getCamera3D()->setManipulator(std::shared_ptr<DRAWER::CameraManipulator>(new DRAWER::FirstPersonCameraManipulator(drawer.getCamera3D())));
  // drawer.addPlugin(std::shared_ptr<DRAWER::Plugin>(new DRAWER::ImGuiPlugin([&]() {
  //   drawer.getCamera3D()->getManipulator()->imGuiCallback();
  // })));
  bool sim=false;
  drawer.setFrameFunc([&](std::shared_ptr<DRAWER::SceneNode>& root) {
    if(sim) {
      xpbd.step();
      visualizeOrUpdateGeometry(*geometry,shapeGeometry);
      visualizeOrUpdateCollision(*geometry,xpbd.getDetector(),shapeCollision);
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
