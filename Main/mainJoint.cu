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
  constexpr std::size_t N=5;
  std::vector<Capsule<T>> ps;

  std::mt19937 mt(123456789);
  std::uniform_real_distribution<T> uni(0, 1);

  T len = uni(mt);
  T radius = uni(mt)/5.f;
  // Joint, which includes 5 capsules.
  QuatT fixedRotation(sqrt(24)/5,0,0,1/5.0);
  fixedRotation.normalize();
  QuatT initQ(1,0,0,0);
  auto prevX = Vec3T(-2,1,0);
  for(int i=0; i<N; i++) {
    Capsule<T> c;
    c._len = len;
    c._radius = radius;
    QuatT q = fixedRotation * initQ;
    q.normalize();
    initQ = q;
    c._q = q;
    Vec3T trans;
    trans = prevX + q.toRotationMatrix()*Vec3T(len/2,0,0);
    prevX = prevX + q.toRotationMatrix()*Vec3T(len,0,0);
    c._x = trans;
    c._x.z() = 0;
    c.initInertiaTensor();
    c._force = Vec3T(0, -9.8f*c._mass,0);
    c._isDynamic = true;
    ps.push_back(c);
  }

  prevX = Vec3T(2, 1, 0);
  for(int i=0; i<N; i++) {
    Capsule<T> c;
    c._len = len;
    c._radius = radius;
    QuatT q(uni(mt),0,0,uni(mt));
    q.normalize();
    c._q = q;
    Vec3T trans;
    trans = prevX + q.toRotationMatrix()*Vec3T(len/2,0,0);
    prevX = prevX + q.toRotationMatrix()*Vec3T(len,0,0);
    c._x = trans;
    c._x.z() = 0;
    c.initInertiaTensor();
    c._force = Vec3T(0, -9.8f*c._mass,0);
    c._isDynamic = true;
    ps.push_back(c);
  }

  // boundary
  Capsule<T> b_1;
  b_1._len = 20;
  b_1._radius = 1;
  b_1._mass = 1;
  b_1._x = Vec3T(0,-4,0);
  b_1._q = QuatT(1,0,0,0);
  b_1.initInertiaTensor();
  b_1._isDynamic = false;
  ps.push_back(b_1);

  Capsule<T> b_2;
  b_2._len = 18.;
  b_2._radius = 1;
  b_2._mass = 1;
  b_2._x = Vec3T(5,0,0);
  b_2._q = QuatT(0.7071f,0,0,0.7071f);
  b_2.initInertiaTensor();
  b_2._isDynamic = false;
  ps.push_back(b_2);

  Capsule<T> b_3;
  b_3._len = 18;
  b_3._radius = 1;
  b_3._mass = 1;
  b_3._x = Vec3T(-5,0,0);
  b_3._q = QuatT(0.7071f,0,0,0.7071f);
  b_3.initInertiaTensor();
  b_3._isDynamic = false;
  ps.push_back(b_3);

  std::shared_ptr<Geometry<T>> geometry(new Geometry<T>);
  geometry->resize(ps.size());
  geometry->setCapsule(ps);
  XPBD<T> xpbd(geometry, 1.0f/60);
  xpbd.addJoint(0,1,ps[0].maxCorner(),ps[1].minCorner());
  xpbd.addJoint(1,2,ps[1].maxCorner(),ps[2].minCorner());
  xpbd.addJoint(2,3,ps[2].maxCorner(),ps[3].minCorner());
  xpbd.addJoint(3,4,ps[3].maxCorner(),ps[4].minCorner());
  // Because the fixedRotation is from A to B, we need to inverse it.
  xpbd.addJointAngular(0, 1, fixedRotation.conjugate());
  xpbd.addJointAngular(1, 2, fixedRotation.conjugate());
  xpbd.addJointAngular(2, 3, fixedRotation.conjugate());
  xpbd.addJointAngular(3, 4, fixedRotation.conjugate());


  xpbd.addJoint(5,6,ps[5].maxCorner(),ps[6].minCorner());
  xpbd.addJoint(6,7,ps[6].maxCorner(),ps[7].minCorner());
  xpbd.addJoint(7,8,ps[7].maxCorner(),ps[8].minCorner());
  xpbd.addJoint(8,9,ps[8].maxCorner(),ps[9].minCorner());

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
