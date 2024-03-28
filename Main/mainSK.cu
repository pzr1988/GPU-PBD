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
  auto animationData=PHYSICSMOTION::AnimationData<T>("/data/GPU-PBD/SKParser/animation.data", "/data/GPU-PBD/SKParser/root_translation.data", "/data/GPU-PBD/SKParser/parent_indices.data");
  QuatT rootLocalQ;
  Vec3T rootLocalX;
  mjcfParser.getRootPos(rootLocalQ, rootLocalX);
  animationData.moveToMJCFRootPos(rootLocalQ, rootLocalX);
  mjcfParser.modifyInitGestureByAnimation(animationData);
  //Shapes list
  std::vector<Shape<T>> ps;
  mjcfParser.getShape(ps);

  //add floor
  Shape<T> floor;
  floor._type=ShapeType::Box;
  floor._len=30;
  floor._width=5;
  floor._height=1;
  floor._x=Vec3T(0,0,-0.5);
  floor._q=QuatT(1,0,0,0);
  floor.initInertiaTensor();
  floor._isDynamic=false;
  ps.push_back(floor);

  std::shared_ptr<Geometry<T>> geometry(new Geometry<T>);
  geometry->resize(ps.size());
  geometry->setShape(ps);
  XPBD<T> xpbd(geometry, 1.0f/60);

  //addPositionConstraint
  std::vector<PHYSICSMOTION::PositionConstraint<T>> pc;
  mjcfParser.getPositionConstraint(pc);
  for(auto& j:pc)
    xpbd.addJoint(j._cA,j._cB,j._cAPos,j._cBPos);

  //addAngularConstraint
  std::vector<PHYSICSMOTION::AngularConstraint<T>> ac;
  mjcfParser.getAngularConstraint(ac);
  for(auto& j:ac)
    xpbd.addJointAngular(j._cA,j._cB,j._psQ,1e-4f,j._sQ,j._pQ);
  //AddAnimationData
  xpbd.addAnimation(animationData._frameNum, animationData._animation.begin(), animationData._animation.end(), animationData._rootQ.begin(), animationData._rootQ.end(), animationData._rootX.begin(), animationData._rootX.end());

  //addGoupInfo, shapes in the same group will not collide
  std::vector<std::pair<int,int>> groupLinks= {{0,7},{7,8},{8,9},{9,10},{10,11},{9,14},{14,15},{9,16},{16,17},{1,4}};
  for(const auto& g:groupLinks)
    xpbd.addGroupLink(g.first,g.second);

  DRAWER::Drawer drawer(argc,argv);
  drawer.addPlugin(std::shared_ptr<DRAWER::Plugin>(new DRAWER::CameraExportPlugin(GLFW_KEY_2,GLFW_KEY_3,"camera.dat")));
  drawer.addPlugin(std::shared_ptr<DRAWER::Plugin>(new DRAWER::CaptureGIFPlugin(GLFW_KEY_1,"record.gif",drawer.FPS())));
  auto shapeGeometry=visualizeOrUpdateGeometry(*geometry);
  auto shapeCollision=visualizeOrUpdateCollision(*geometry,xpbd.getDetector(),xpbd.getJointPositions());
  drawer.addShape(shapeGeometry);
  drawer.addShape(shapeCollision);
  drawer.addCamera3D(90,Eigen::Matrix<GLfloat,3,1>(0,0,1),Eigen::Matrix<GLfloat,3,1>(0,2,1),Eigen::Matrix<GLfloat,3,1>(0,-1,0));
  drawer.getCamera3D()->setManipulator(std::shared_ptr<DRAWER::CameraManipulator>(new DRAWER::FirstPersonCameraManipulator(drawer.getCamera3D())));
  drawer.addPlugin(std::shared_ptr<DRAWER::Plugin>(new DRAWER::ImGuiPlugin([&]() {
    drawer.getCamera3D()->getManipulator()->imGuiCallback();
  })));
#define USE_LIGHT
#ifdef USE_LIGHT
  drawer.addLightSystem(2048,20);
  drawer.getLight()->lightSz(10);
  drawer.getLight()->addLight(Eigen::Matrix<GLfloat,3,1>(0,0,3),
                              Eigen::Matrix<GLfloat,3,1>(1,1,1),
                              Eigen::Matrix<GLfloat,3,1>(1,1,1),
                              Eigen::Matrix<GLfloat,3,1>(0,0,0));
#endif
  bool sim=false;
  drawer.setFrameFunc([&](std::shared_ptr<DRAWER::SceneNode>& root) {
    if(sim) {
      xpbd.step();
      visualizeOrUpdateGeometry(*geometry,shapeGeometry);
      visualizeOrUpdateCollision(*geometry,xpbd.getDetector(),xpbd.getJointPositions(),shapeCollision);
    }
  });
  drawer.setKeyFunc([&](GLFWwindow* wnd,int key,int scan,int action,int mods,bool captured) {
    if(captured)
      return;
    else if(key==GLFW_KEY_R && action==GLFW_PRESS)
      sim=!sim;
  });
  drawer.mainLoop();

  return 0;
}
