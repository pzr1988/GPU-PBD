#include <PBD/Utils.h>
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
  constexpr std::size_t N=10;
  std::vector<Capsule<T>> ps(N);

  std::mt19937 mt(123456789);
  std::uniform_real_distribution<T> uni(0.0, 1.0);

  for(auto& p:ps) {
    p._len=uni(mt);
    p._radius=uni(mt)/3.;
    p._mass = 3.14*p._radius*p._radius*p._len+3.14*4.0/3.0*p._radius*p._radius*p._radius;
    p._force = Eigen::Matrix<T,3,1>(0, -9.8*p._mass,0);
    Eigen::Quaternion<T> q(uni(mt),uni(mt),uni(mt),uni(mt));
    q.normalize();
    Eigen::AngleAxis<T> euler(q);
    euler.axis().x()=0;
    euler.axis().y()=0;
    euler.axis().normalize();
    q = Eigen::Quaternion<T>(euler);
    q.normalize();
    p._q = q;
    Eigen::Matrix<T,3,1> trans;
    trans.setRandom();
    p._x = trans;
    p._x.z() = 0;
    p.initInertiaTensor();
    p._isDynamic = true;
  }

  // boundary
  Capsule<T> b_1;
  b_1._len = 10.;
  b_1._radius = 1.0;
  b_1._mass = 1.0;
  b_1._x = Eigen::Matrix<T,3,1>(0,-4,0);
  b_1._q = Eigen::Quaternion<T>(1,0,0,0);
  b_1.initInertiaTensor();
  b_1._isDynamic = false;
  ps.push_back(b_1);

  Capsule<T> b_2;
  b_2._len = 8.;
  b_2._radius = 1.0;
  b_2._mass = 1.0;
  b_2._x = Eigen::Matrix<T,3,1>(5,0,0);
  b_2._q = Eigen::Quaternion<T>(0.7071,0,0,0.7071);
  b_2.initInertiaTensor();
  b_2._isDynamic = false;
  ps.push_back(b_2);

  Capsule<T> b_3;
  b_3._len = 8;
  b_3._radius = 1.0;
  b_3._mass = 1.0;
  b_3._x = Eigen::Matrix<T,3,1>(-5,0,0);
  b_3._q = Eigen::Quaternion<T>(0.7071,0,0,0.7071);
  b_3.initInertiaTensor();
  b_3._isDynamic = false;
  ps.push_back(b_3);

  std::shared_ptr<GPUPBD::Geometry<T>> geometry(new GPUPBD::Geometry<T>);
  geometry->reserve(ps.size());
  geometry->resize(ps.size());
  geometry->setCapsule(ps);

  // GPUPBD::CollisionDetector<T> detector(geometry);
  // detector.detectCollisions();
  DRAWER::Drawer drawer(argc,argv);
  drawer.addPlugin(std::shared_ptr<DRAWER::Plugin>(new DRAWER::CameraExportPlugin(GLFW_KEY_2,GLFW_KEY_3,"camera.dat")));
  drawer.addPlugin(std::shared_ptr<DRAWER::Plugin>(new DRAWER::CaptureGIFPlugin(GLFW_KEY_1,"record.gif",drawer.FPS())));
  GPUPBD::XPBD<T> xpbd(geometry, 1.0/drawer.FPS());
  auto shapeGeometry=visualizeOrUpdateGeometry(*geometry);
  auto shapeCollision=visualizeOrUpdateCollision(*geometry,xpbd.getDetector());
  drawer.addShape(shapeGeometry);
  drawer.addShape(shapeCollision);
  drawer.addCamera3D(90,Eigen::Matrix<GLfloat,3,1>(0,1,0),
    Eigen::Matrix<GLfloat,3,1>(0,0,5),
    Eigen::Matrix<GLfloat,3,1>(0,0,-1));
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
