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
#include <queue>
#include <cmath>
#include <SKParser/MujocoXmlParser.h>

using namespace GPUPBD;

int main(int argc,char** argv) {
  typedef LSCALAR T;
  DECL_MAT_VEC_MAP_TYPES_T
  MujocoXmlParser parser;
  BodyParam* pRootBody = parser.Load("/data/GPU-PBD/SKParser/SK_Mannequin_PhysicsAsset_ABFB4_MJCF.xml");

  std::vector<Capsule<T>> ps;
  std::queue<BodyParam*> q;
  q.push(pRootBody);
  while(0 != q.size()) {
    BodyParam* pNow = q.front();
    q.pop();
    if(pNow->pChild) {
        q.push(pNow->pChild);
    }
    if(pNow->pSibling) {
        q.push(pNow->pSibling);
    }
    if(pNow->pGeoms){
        if(0 != strcmp(pNow->pGeoms->strType, "capsule")){
            continue;
        }
        Capsule<T> b;
        b._x = Vec3T(pNow->posGlobal[0], pNow->posGlobal[1], pNow->posGlobal[2]);
        b._len = pNow->pGeoms->fHalfHeight*2;
        b._radius = pNow->pGeoms->fRadius;
        T sinValue = sin(pNow->pGeoms->qOrient[0]*3.1415/360.0);
        T cosValue = cos(pNow->pGeoms->qOrient[0]*3.1415/360.0);
        b._q = QuatT(cosValue,sinValue*pNow->pGeoms->qOrient[1],sinValue*pNow->pGeoms->qOrient[2], sinValue*pNow->pGeoms->qOrient[3]);
        printf("Capusle pos:%f,%f,%f, radius:%f, len:%f, q:%f,%f,%f,%f \n",
            b._x[0], b._x[1],b._x[2],
            b._radius, b._len,
            b._q.w(), b._q.x(),b._q.y(),b._q.z());
        ps.push_back(b);
    }
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
//   ps.push_back(b_1);

  Capsule<T> b_2;
  b_2._len = 18.;
  b_2._radius = 1;
  b_2._mass = 1;
  b_2._x = Vec3T(5,0,0);
  b_2._q = QuatT(0.7071f,0,0,0.7071f);
  b_2.initInertiaTensor();
  b_2._isDynamic = false;
//   ps.push_back(b_2);

  Capsule<T> b_3;
  b_3._len = 18;
  b_3._radius = 1;
  b_3._mass = 1;
  b_3._x = Vec3T(-5,0,0);
  b_3._q = QuatT(0.7071f,0,0,0.7071f);
  b_3.initInertiaTensor();
  b_3._isDynamic = false;
//   ps.push_back(b_3);

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
//   drawer.addCamera3D(90,Eigen::Matrix<GLfloat,3,1>(0,1,0),Eigen::Matrix<GLfloat,3,1>(0,0,2),Eigen::Matrix<GLfloat,3,1>(0,0,-1));
//   drawer.getCamera3D()->setManipulator(std::shared_ptr<DRAWER::CameraManipulator>(new DRAWER::FirstPersonCameraManipulator(drawer.getCamera3D())));
//   drawer.addPlugin(std::shared_ptr<DRAWER::Plugin>(new DRAWER::ImGuiPlugin([&]() {
//     drawer.getCamera3D()->getManipulator()->imGuiCallback();
//   })));
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

  delete pRootBody;
  return 0;
}
