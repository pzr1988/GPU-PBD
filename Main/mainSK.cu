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
#include <fstream>

using namespace GPUPBD;

int main(int argc,char** argv) {
  typedef LSCALAR T;
  DECL_MAT_VEC_MAP_TYPES_T

  std::vector<Shape<T>> ps; //shapes list
  auto mjcfParser=PHYSICSMOTION::MJCFParser<T>("/data/GPU-PBD/SKParser/MarathonCharacter_PhysicsAsset2.xml");
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

  // addGoupInfo, shapes in the same group will not collide.
  std::vector<std::pair<int, int>> groupLinks= {{0,7},{7,8},{8,9},{9,10},{10,11},{9,14},{14,15},{9,16},{16,17},{1,4}};
  for(const auto& g : groupLinks) {
    xpbd.addGroupLink(g.first, g.second);
  }

  std::cout<<"==========================animation info==========================" << std::endl;
  std::ifstream inFile("/data/GPU-PBD/SKParser/animation.data");
  if (!inFile) {
    std::cerr << "Unable to open file" << std::endl;
    return 1;
  }
  std::vector<std::vector<double>> rawData;
  std::vector<QuatT> animation;
  std::vector<QuatT> rootQ;
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
  int frameNum = rawData.size() / numJoints;
  animation.resize(rawData.size());
  rootQ.resize(frameNum);
  // TODO fix hardcode of root local rotation.
  QuatT rootLocalQ = QuatT::FromTwoVectors(Vec3T::UnitX(),Vec3T(0.0237-0.0335, -0.0861-0.0849, -0.0278-(-0.0278)));
  for(int i=0; i<rawData.size(); i++) {
    const auto& v = rawData.at(i);
    QuatT q(v[0],v[1],v[2],v[3]);
    animation[i]=q;
    if(i%numJoints==0) {
      rootQ[i/numJoints]=q*rootLocalQ;
    }
  }
  std::ifstream inFileRoot("/data/GPU-PBD/SKParser/root_translation.data");
  if (!inFileRoot) {
    std::cerr << "Unable to open file" << std::endl;
    return 1;
  }
  std::vector<Vec3T> rootX;
  while (!inFileRoot.eof()) {
    std::vector<double> row;
    for (int i = 0; i < 3; ++i) {
      if (inFileRoot >> value) {
        row.push_back(value);
      }
    }
    if (!row.empty()) {
      // TODO fix hardcode of root local translation.
      rootX.push_back(Vec3T(row[0],row[1],row[2])-Vec3T(3.1954-0.0208, -0.0175, 0.9931-1.0399));
    }
  }
  inFileRoot.close();

  xpbd.addAnimation(frameNum, animation.begin(), animation.end(), rootQ.begin(), rootQ.end(), rootX.begin(), rootX.end());

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
