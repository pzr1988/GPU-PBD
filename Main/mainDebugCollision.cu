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
    Eigen::Quaternion<T> q(uni(mt),uni(mt),uni(mt),uni(mt));
    q.normalize();
    p._q = q;
    Eigen::Matrix<T,3,1> trans;
    trans.setRandom();
    p._x = trans;
    p.initInertiaTensor();
    p._isDynamic = true;
  }

  std::shared_ptr<GPUPBD::Geometry<T>> geometry(new GPUPBD::Geometry<T>);
  geometry->reserve(ps.size());
  geometry->resize(ps.size());
  geometry->setCapsule(ps);

  // GPUPBD::CollisionDetector<T> detector(geometry);
  // detector.detectCollisions();
  GPUPBD::XPBD<T> xpbd(geometry, 0.1);
  xpbd.step();

  DRAWER::Drawer drawer(argc,argv);
  drawer.addPlugin(std::shared_ptr<DRAWER::Plugin>(new DRAWER::CameraExportPlugin(GLFW_KEY_2,GLFW_KEY_3,"camera.dat")));
  drawer.addPlugin(std::shared_ptr<DRAWER::Plugin>(new DRAWER::CaptureGIFPlugin(GLFW_KEY_1,"record.gif",drawer.FPS())));
  auto shapeGeometry=visualizeOrUpdateGeometry(*geometry);
  auto shapeCollision=visualizeOrUpdateCollision(*geometry,xpbd.getDetector());
  drawer.addShape(shapeGeometry);
  drawer.addShape(shapeCollision);
  drawer.addCamera3D(90,Eigen::Matrix<GLfloat,3,1>(0,1,0));
  drawer.getCamera3D()->setManipulator(std::shared_ptr<DRAWER::CameraManipulator>(new DRAWER::FirstPersonCameraManipulator(drawer.getCamera3D())));
  drawer.addPlugin(std::shared_ptr<DRAWER::Plugin>(new DRAWER::ImGuiPlugin([&]() {
    drawer.getCamera3D()->getManipulator()->imGuiCallback();
  })));
  drawer.mainLoop();
  return 0;
}
