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
struct Body {
  DECL_MAT_VEC_MAP_TYPES_T
  int _parent;
  int _depth;
  T _radius;
  Vec3T _x;
  QuatT _q;
  Vec6T _ft;
  std::string _name;
  Capsule<T> _c;
};
template <typename T>
void readBodies(std::vector<Body<T>>& bodies, int parentId, const tinyxml2::XMLElement* g) {
  DECL_MAT_VEC_MAP_TYPES_T
  Body<T> body;
  body._parent=parentId;
  body._name=g->Attribute("name");
  //compute depth
  body._depth=parentId==-1?0:bodies[parentId]._depth+1;
  //read trans
  {
    body._x=PHYSICSMOTION::parsePtreeDef<Vec3T>(*g,"<xmlattr>.pos","0 0 0");
    Vec4T tmpQ=PHYSICSMOTION::parsePtreeDef<Vec4T>(*g,"<xmlattr>.quat","1 0 0 0");
    body._q=QuatT(tmpQ[0],tmpQ[1],tmpQ[2],tmpQ[3]);
  }
  //read geometry
  if(g->FirstChildElement("geom")->FindAttribute("type") == NULL) {
    //capsule
    const tinyxml2::XMLElement* gg=g->FirstChildElement("geom");
    body._ft=PHYSICSMOTION::parsePtreeDef<Vec6T>(*gg,"<xmlattr>.fromto","0 0 0 0 0 0");
    body._radius=PHYSICSMOTION::get<T>(*gg,"<xmlattr>.size");
    int _parentId=(int)bodies.size();
    bodies.push_back(body);
    for(const tinyxml2::XMLElement* gc=g->FirstChildElement(); gc; gc=gc->NextSiblingElement())
      if(std::string(gc->Name()) == "body")
        readBodies(bodies,_parentId,gc);
  }
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
void updateCapsule(std::vector<Body<T>>& bodies) {
  DECL_MAT_VEC_MAP_TYPES_T
  for(auto& body : bodies) {
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
              << std::endl;
  }

  updateCapsule(bodies);
  std::cout<<"==========================updated info==========================" << std::endl;
  for(int i=0; i<bodies.size(); i++) {
    Body<T>& b = bodies[i];
    std::string indent(b._depth * 2, ' ');
    auto tmpx = b._c._q.toRotationMatrix() * Vec3T(b._c._len, 0, 0);
    std::cout << indent  << b._name <<": x:" << b._x[0] << " " << b._x[1] << " " << b._x[2]
              << ", q:" << b._q.w() << " " << b._q.x() << " " << b._q.y() << " " << b._q.z()
              << ", capsule radius: " << b._c._radius << ", len:" << b._c._len
              <<", x:" << b._c._x[0] << " " << b._c._x[1] << " " << b._c._x[2]
              << ", q:" << b._c._q.w() << " " << b._c._q.x() << " " << b._c._q.y() << " " << b._c._q.z()
              << ", parent Name: " << bodies[b._parent>-1?b._parent:0]._name
              << std::endl;
  }

  std::vector<Capsule<T>> ps;
  for(auto& b : bodies) {
    Capsule<T> c = b._c;
    ps.push_back(c);
  }

  std::shared_ptr<Geometry<T>> geometry(new Geometry<T>);
  geometry->resize(ps.size());
  geometry->setCapsule(ps);
  XPBD<T> xpbd(geometry, 1.0f/60);
  // TODO addJoint
  DRAWER::Drawer drawer(argc,argv);
  drawer.addPlugin(std::shared_ptr<DRAWER::Plugin>(new DRAWER::CameraExportPlugin(GLFW_KEY_2,GLFW_KEY_3,"camera.dat")));
  drawer.addPlugin(std::shared_ptr<DRAWER::Plugin>(new DRAWER::CaptureGIFPlugin(GLFW_KEY_1,"record.gif",drawer.FPS())));
  auto shapeGeometry=visualizeOrUpdateGeometry(*geometry);
  auto shapeCollision=visualizeOrUpdateCollision(*geometry,xpbd.getDetector());
  drawer.addShape(shapeGeometry);
  drawer.addShape(shapeCollision);
  // drawer.addCamera3D(90,Eigen::Matrix<GLfloat,3,1>(0,1,0),Eigen::Matrix<GLfloat,3,1>(0,0,2),Eigen::Matrix<GLfloat,3,1>(0,0,-1));
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
