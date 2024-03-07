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
  std::vector<Shape<T>> cpuGeometry(g.size());
  thrust::copy(g.begin(),g.end(),cpuGeometry.begin());
  //update
  std::shared_ptr<DRAWER::CompositeShape> ret=s;
  if(!ret) {
    ret.reset(new DRAWER::CompositeShape());
    for(int i=0; i<(int)cpuGeometry.size(); i++) {
      std::shared_ptr<DRAWER::Bullet3DShape> c(new DRAWER::Bullet3DShape);
      if(cpuGeometry[i].isCapsule()) {
        c->addShape(DRAWER::makeSphericalBox(RES,fill,cpuGeometry[i]._radius,Eigen::Matrix<GLfloat,3,1>((GLfloat)cpuGeometry[i]._len/2.0f,0,0)));
      } else if(cpuGeometry[i].isBox()) {
        c->addShape(DRAWER::makeBox(RES,fill,Eigen::Matrix<GLfloat,3,1>((GLfloat)cpuGeometry[i]._len/2.0f,(GLfloat)cpuGeometry[i]._width/2.0f,(GLfloat)cpuGeometry[i]._height/2.0f)));
      } else if(cpuGeometry[i].isSphere()) {
        c->addShape(DRAWER::makeSphere(RES,fill,(GLfloat)cpuGeometry[i]._radius));
      }
      c->setColorDiffuse(GL_LINES,.7f,.7f,.7f);
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
std::shared_ptr<DRAWER::MeshShape> visualizeOrUpdateCollision(const Geometry<T>& g,const CollisionDetector<T>& cd,std::shared_ptr<DRAWER::MeshShape> s=NULL,int width=5) {
  //geometry
  std::vector<Shape<T>> cpuGeometry(g.size());
  thrust::copy(g.begin(),g.end(),cpuGeometry.begin());
  //collision
  std::vector<Constraint<T>> cpuCollision(cd.size());
  thrust::copy(cd.begin(),cd.end(),cpuCollision.begin());
  //visualize
  std::shared_ptr<DRAWER::MeshShape> ret=s;
  if(!ret)
    ret.reset(new DRAWER::MeshShape);
  else ret->clear();
  for(int i=0; i<(int)cpuCollision.size(); i++) {
    const auto& xA=cpuGeometry[cpuCollision[i]._shapeIdA]._x;
    const auto& pA=cpuGeometry[cpuCollision[i]._shapeIdA]._q;
    const auto& xB=cpuGeometry[cpuCollision[i]._shapeIdB]._x;
    const auto& pB=cpuGeometry[cpuCollision[i]._shapeIdB]._q;
    ret->addVertex(pA.toRotationMatrix()*cpuCollision[i]._localPointA+xA);
    ret->addVertex(pB.toRotationMatrix()*cpuCollision[i]._localPointB+xB);
    ret->addIndexSingle(i*2+0);
    ret->addIndexSingle(i*2+1);
  }
  ret->setMode(GL_LINES);
  ret->setLineWidth((GLfloat)width);
  ret->setColorDiffuse(GL_LINES,1,0,0);
  return ret;
}
template <typename T>
std::shared_ptr<DRAWER::MeshShape> visualizeOrUpdateCollision(const Geometry<T>& g,const CollisionDetector<T>& cd,const thrust::device_vector<Constraint<T>>& joint,std::shared_ptr<DRAWER::MeshShape> s=NULL,int width=5) {
  //geometry
  std::vector<Shape<T>> cpuGeometry(g.size());
  thrust::copy(g.begin(),g.end(),cpuGeometry.begin());
  //collision
  std::vector<Constraint<T>> cpuCollision(cd.size()+joint.size());
  thrust::copy(cd.begin(),cd.end(),cpuCollision.begin());
  thrust::copy(joint.begin(),joint.end(),cpuCollision.begin()+cd.size());
  //visualize
  std::shared_ptr<DRAWER::MeshShape> ret=s;
  if(!ret)
    ret.reset(new DRAWER::MeshShape);
  else ret->clear();
  for(int i=0; i<(int)cpuCollision.size(); i++) {
    const auto& xA=cpuGeometry[cpuCollision[i]._shapeIdA]._x;
    const auto& pA=cpuGeometry[cpuCollision[i]._shapeIdA]._q;
    const auto& xB=cpuGeometry[cpuCollision[i]._shapeIdB]._x;
    const auto& pB=cpuGeometry[cpuCollision[i]._shapeIdB]._q;
    ret->addVertex(pA.toRotationMatrix()*cpuCollision[i]._localPointA+xA);
    ret->addVertex(pB.toRotationMatrix()*cpuCollision[i]._localPointB+xB);
    ret->addIndexSingle(i*2+0);
    ret->addIndexSingle(i*2+1);
  }
  ret->setMode(GL_LINES);
  ret->setLineWidth((GLfloat)width);
  ret->setColorDiffuse(GL_LINES,1,0,0);
  return ret;
}
}
#endif
