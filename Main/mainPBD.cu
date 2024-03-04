#include <PBD/Geometry.h>
#include <PBD/Collision.h>
#include <PBD/XPBD.h>
#include <PBD/Visualizer.h>
#include <TinyVisualizer/FirstPersonCameraManipulator.h>
#include <TinyVisualizer/CameraExportPlugin.h>
#include <TinyVisualizer/CaptureGIFPlugin.h>
#include <TinyVisualizer/ImGuiPlugin.h>
#include <TinyVisualizer/Camera3D.h>
#include <random>

using namespace GPUPBD;

int main(int argc,char** argv) {
  typedef LSCALAR T;
  DECL_MAT_VEC_MAP_TYPES_T
  std::vector<Shape<T>> ps;

  std::mt19937 mt(123456789);
  std::uniform_real_distribution<T> uni(0, 1);

  T len = uni(mt);
  T radius = uni(mt)/3.f;
  for(int i=0; i<5; i++) {
    Shape<T> p;
    p._type = ShapeType::Capsule;
    p._len=len;
    p._radius=radius;
    QuatT q(uni(mt),0,0,uni(mt));
    q.normalize();
    p._q = q;
    Vec3T trans;
    trans.setRandom();
    p._x = trans*3;
    p._x.z() = 0;
    p.initInertiaTensor();
    p._force = Vec3T(0, -9.8f*p._mass,0);
    p._isDynamic = true;
    ps.push_back(p);
  }

  len = uni(mt);
  for(int i=0; i<5; i++) {
    Shape<T> p;
    p._type = ShapeType::Box;
    p._len=len;
    p._width=len;
    p._height=len;
    // p._radius=radius;
    QuatT q(uni(mt),0,0,uni(mt));
    q.normalize();
    p._q = q;
    Vec3T trans;
    trans.setRandom();
    p._x = trans*3;
    p._x.z() = 0;
    p.initInertiaTensor();
    p._force = Vec3T(0, -9.8f*p._mass,0);
    p._isDynamic = true;
    ps.push_back(p);
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

  DRAWER::Drawer drawer(argc,argv);
  drawer.addPlugin(std::shared_ptr<DRAWER::Plugin>(new DRAWER::CameraExportPlugin(GLFW_KEY_2,GLFW_KEY_3,"camera.dat")));
  drawer.addPlugin(std::shared_ptr<DRAWER::Plugin>(new DRAWER::CaptureGIFPlugin(GLFW_KEY_1,"record.gif",drawer.FPS())));
  auto shapeGeometry=visualizeOrUpdateGeometry(*geometry);
  auto shapeCollision=visualizeOrUpdateCollision(*geometry,xpbd.getDetector());
  drawer.addShape(shapeGeometry);
  drawer.addShape(shapeCollision);
  drawer.addCamera3D(90,Eigen::Matrix<GLfloat,3,1>(0,1,0),Eigen::Matrix<GLfloat,3,1>(0,0,5),Eigen::Matrix<GLfloat,3,1>(0,0,-1));
  drawer.getCamera3D()->setManipulator(std::shared_ptr<DRAWER::CameraManipulator>(new DRAWER::FirstPersonCameraManipulator(drawer.getCamera3D())));
  drawer.addPlugin(std::shared_ptr<DRAWER::Plugin>(new DRAWER::ImGuiPlugin([&]() {
    drawer.getCamera3D()->getManipulator()->imGuiCallback();
  })));
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
