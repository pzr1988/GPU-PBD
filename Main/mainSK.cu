#include <PBD/Geometry.h>
#include <PBD/Collision.h>
#include <PBD/XPBD.h>
#include <PBD/Visualizer.h>
#include <TinyVisualizer/FirstPersonCameraManipulator.h>
#include <TinyVisualizer/CameraExportPlugin.h>
#include <TinyVisualizer/CaptureGIFPlugin.h>
#include <TinyVisualizer/ImGuiPlugin.h>
#include <TinyVisualizer/Camera3D.h>
#include <SKParser/MJCFParser.h>
#include <SKParser/AnimationParser.h>
#include <fstream>

using namespace GPUPBD;

int main(int argc,char** argv) {
  typedef LSCALAR T;
  DECL_MAT_VEC_MAP_TYPES_T

  //MJCF Info
  auto mjcfParser=PHYSICSMOTION::MJCFParser<T>("/data/GPU-PBD/SKParser/MarathonCharacter_PhysicsAsset2.xml");
  //Animation Data
  auto animationData=PHYSICSMOTION::AnimationData<T>("/data/GPU-PBD/SKParser/animation.data", "/data/GPU-PBD/SKParser/root_translation.data");
  mjcfParser.modifyInitPosByAnimation(animationData);
  //Shapes list
  std::vector<Shape<T>> ps;
  mjcfParser.getShape(ps);

  // Add floor to shapes list
  Shape<T> b_1;
  b_1._type = ShapeType::Box;
  b_1._len=30;
  b_1._width=30;
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
  std::vector<PHYSICSMOTION::PositionConstraint<T>> pc;
  mjcfParser.getPositionConstraint(pc);
  for(auto& j: pc) {
    xpbd.addJoint(j._cA,j._cB,j._cAPos,j._cBPos);
  }

  // addAngularConstraint
  std::vector<PHYSICSMOTION::AngularConstraint<T>> ac;
  mjcfParser.getAngularConstraint(ac);
  for(auto& j: ac) {
    xpbd.addJointAngular(j._cA,j._cB,j._psQ, 0.0001f, j._sQ, j._pQ);
  }
  //AddAnimationData
  xpbd.addAnimation(animationData._frameNum, animationData._animation.begin(), animationData._animation.end(), animationData._rootQ.begin(), animationData._rootQ.end(), animationData._rootX.begin(), animationData._rootX.end());

  // addGoupInfo, shapes in the same group will not collide.
  // TODO remove hardcode
  std::vector<std::pair<int, int>> groupLinks= {{0,7},{7,8},{8,9},{9,10},{10,11},{9,14},{14,15},{9,16},{16,17},{1,4}};
  for(const auto& g : groupLinks) {
    xpbd.addGroupLink(g.first, g.second);
  }

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
