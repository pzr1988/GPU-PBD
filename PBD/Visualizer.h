#ifndef VISUALIZER_H
#define VISUALIZER_H
#include "Geometry.h"
#include "Collision.h"
#include <TinyVisualizer/MakeMesh.h>
#include <TinyVisualizer/MeshShape.h>
#include <TinyVisualizer/Bullet3DShape.h>

namespace GPUPBD {
template <typename T>
std::shared_ptr<DRAWER::CompositeShape> visualizeOrUpdateGeometry(const Geometry<T>& g,std::shared_ptr<DRAWER::CompositeShape> s=NULL,int RES=8,bool fill=false) {
  //geometry
  std::vector<Capsule<T>> cpuGeometry(g.getCapsules().size());
  thrust::copy(g.getCapsules().begin(),g.getCapsules().end(),cpuGeometry.begin());
  //update
  std::shared_ptr<DRAWER::CompositeShape> ret=s;
  if(!ret) {
    ret.reset(new DRAWER::CompositeShape());
    for(int i=0; i<(int)cpuGeometry.size(); i++) {
      std::shared_ptr<DRAWER::Bullet3DShape> c(new DRAWER::Bullet3DShape);
      c->addShape(DRAWER::makeSphericalBox(RES,fill,cpuGeometry[i]._radius,Eigen::Matrix<GLfloat,3,1>(cpuGeometry[i]._len/2.0,0,0)));
      c->setColorDiffuse(GL_LINES,.7,.7,.7);
      c->setLineWidth(5);
      ret->addShape(c);
    }
  }
  //visualize
  for(int i=0; i<(int)cpuGeometry.size(); i++) {
    auto c=std::dynamic_pointer_cast<DRAWER::Bullet3DShape>(ret->getChild(i));
    c->setLocalRotate(cpuGeometry[i]._q.toRotationMatrix());
    c->setLocalTranslate(cpuGeometry[i]._x);
  }
  return ret;
}
template <typename T>
std::shared_ptr<DRAWER::MeshShape> visualizeOrUpdateCollision(const Geometry<T>& g,const CollisionDetector<T>& cd,int width=5) {
  //geometry
  std::vector<Capsule<T>> cpuGeometry(g.getCapsules().size());
  thrust::copy(g.getCapsules().begin(),g.getCapsules().end(),cpuGeometry.begin());
  //collision
  std::vector<Collision<T>> cpuCollision(cd.getCollisions().size());
  thrust::copy(cd.getCollisions().begin(),cd.getCollisions().end(),cpuCollision.begin());
  //visualize
  std::shared_ptr<DRAWER::MeshShape> ret(new DRAWER::MeshShape);
  for(int i=0; i<(int)cpuCollision.size(); i++) {
    ret->addVertex(cpuCollision[i]._gobalPointA);
    ret->addVertex(cpuCollision[i]._gobalPointB);
    ret->addIndexSingle(i*2+0);
    ret->addIndexSingle(i*2+1);
  }
  ret->setMode(GL_LINES);
  ret->setLineWidth(width);
  ret->setColorDiffuse(GL_LINES,1,0,0);
  return ret;
}
}
#endif
