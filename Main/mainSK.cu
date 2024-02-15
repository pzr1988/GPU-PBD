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
void readJoints(std::vector<Joint<T>>& joints, int parentId, const tinyxml2::XMLElement* g) {
  DECL_MAT_VEC_MAP_TYPES_T
  Joint<T> joint;
  joint._parent=parentId;
  joint._name=g->Attribute("name");
  //compute depth
  joint._depth=parentId==-1?0:joints[parentId]._depth+1;
  //read trans
  {
    joint._x=PHYSICSMOTION::parsePtreeDef<Vec3T>(*g,"<xmlattr>.pos","0 0 0");
    Vec4T tmpQ=PHYSICSMOTION::parsePtreeDef<Vec4T>(*g,"<xmlattr>.quat","1 0 0 0");
    joint._q=QuatT(tmpQ[0],tmpQ[1],tmpQ[2],tmpQ[3]);
  }
  //read geometry
  if(g->FirstChildElement("geom")->FindAttribute("type") == NULL) {
    //capsule
    const tinyxml2::XMLElement* gg=g->FirstChildElement("geom");
    joint._ft=PHYSICSMOTION::parsePtreeDef<Vec6T>(*gg,"<xmlattr>.fromto","0 0 0 0 0 0");
    joint._radius=PHYSICSMOTION::get<T>(*gg,"<xmlattr>.size");
    int _parentId=(int)joints.size();
    joints.push_back(joint);
    for(const tinyxml2::XMLElement* gc=g->FirstChildElement(); gc; gc=gc->NextSiblingElement())
      if(std::string(gc->Name()) == "body")
        readJoints(joints,_parentId,gc);
  }
}

template <typename T>
void readMJCF(std::vector<Joint<T>>& joints, const std::string& file) {
  tinyxml2::XMLDocument pt;
  pt.LoadFile(file.c_str());
  tinyxml2::XMLElement* link=pt.RootElement();
  for(const tinyxml2::XMLElement* g=link->FirstChildElement(); g; g=g->NextSiblingElement())
    if(std::string(g->Name()) == "worldbody")
      readJoints(joints,-1,g->FirstChildElement("body"));
}

template <typename T>
void updateCapsule(std::vector<Joint<T>>& joints) {
  DECL_MAT_VEC_MAP_TYPES_T
  for(auto& joint : joints) {
    Vec3T c1 = Vec3T(joint._ft[0], joint._ft[1], joint._ft[2]);
    Vec3T c2 = Vec3T(joint._ft[3], joint._ft[4], joint._ft[5]);
    c1 = joint._x + joint._q.toRotationMatrix() * c1;
    c2 = joint._x + joint._q.toRotationMatrix() * c2;
    if(joint._parent>=0) {
      const auto& p = joints[joint._parent];
      c1 = p._x + p._q.toRotationMatrix()*c1;
      c2 = p._x + p._q.toRotationMatrix()*c2;
      joint._x = p._x + p._q.toRotationMatrix()*joint._x;
      joint._q = (p._q * joint._q).normalized();
    }
    joint._c._radius=joint._radius;
    joint._c._len=(c2-c1).norm();
    joint._c._x=(c1+c2)/2;
    joint._c._q=QuatT::FromTwoVectors(Vec3T::UnitX(),(c2-c1).normalized());
  }
}


int main(int argc,char** argv) {
  typedef LSCALAR T;
  DECL_MAT_VEC_MAP_TYPES_T
  std::vector<Joint<T>> joints;
  readMJCF(joints, "/data/GPU-PBD/SKParser/SK_Mannequin_PhysicsAsset_ABFB4_MJCF.xml");
  std::cout<<"==========================original info==========================" << std::endl;
  for(int i=0; i<joints.size(); i++) {
    Joint<T>& j = joints[i];
    std::string indent(j._depth * 2, ' ');
    std::cout << indent  << j._name <<": x:" << j._x[0] << " " << j._x[1] << " " << j._x[2]
              << ", q:" << j._q.w() << " " << j._q.x() << " " << j._q.y() << " " << j._q.z()
              << ", radius:" << j._radius
              << ", fromto:" << j._ft[0] << " " << j._ft[1] << " "<< j._ft[2] << " "<< j._ft[3] << " "<< j._ft[4] << " "<< j._ft[5]
              << ", parent Name: " << joints[j._parent>-1?j._parent:0]._name
              << std::endl;
  }

  updateCapsule(joints);
  std::cout<<"==========================updated info==========================" << std::endl;
  for(int i=0; i<joints.size(); i++) {
    Joint<T>& j = joints[i];
    std::string indent(j._depth * 2, ' ');
    auto tmpx = j._c._q.toRotationMatrix() * Vec3T(j._c._len, 0, 0);
    std::cout << indent  << j._name <<": x:" << j._x[0] << " " << j._x[1] << " " << j._x[2]
              << ", q:" << j._q.w() << " " << j._q.x() << " " << j._q.y() << " " << j._q.z()
              << ", capsule radius: " << j._c._radius << ", len:" << j._c._len
              <<", x:" << j._c._x[0] << " " << j._c._x[1] << " " << j._c._x[2]
              << ", q:" << j._c._q.w() << " " << j._c._q.x() << " " << j._c._q.y() << " " << j._c._q.z()
              << ", parent Name: " << joints[j._parent>-1?j._parent:0]._name
              << std::endl;
  }

  std::vector<Capsule<T>> ps;
  for(auto& j : joints) {
    Capsule<T> c = j._c;
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
