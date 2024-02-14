#ifndef UTILS_H
#define UTILS_H

#include <PBD/Pragma.h>
#include <filesystem>
#include "SKParser/tinyxml2.h"
#include <string>

namespace PHYSICSMOTION {
extern std::vector<std::string> split(const std::string& l,const std::string& sep=" ");
extern bool beginsWith(const std::string& l,const std::string& s);
extern bool endsWith(const std::string& l,const std::string& s);
extern std::string toUpper(const std::string& l);
inline void hash_combine(std::size_t& seed, size_t hash) {
  seed ^= hash + 0x9e3779b9 + (seed<<6) + (seed>>2);
}
template <typename T>
void sort2(T& a,T& b) {
  if(a>b)
    std::swap(a,b);
}
template <typename T>
void sort3(T& a,T& b,T& c) {
  if(a>b)
    std::swap(a,b);
  if(a>c)
    std::swap(a,c);
  if(b>c)
    std::swap(b,c);
}
template <typename T>
void sort4(T& a,T& b,T& c,T& d) {
  if(a>b)
    std::swap(a,b);
  if(a>c)
    std::swap(a,c);
  if(a>d)
    std::swap(a,d);

  if(b>c)
    std::swap(b,c);
  if(b>d)
    std::swap(b,d);

  if(c>d)
    std::swap(c,d);
};

//filesystem
bool exists(const std::filesystem::path& path);
void removeDir(const std::filesystem::path& path);
void create(const std::filesystem::path& path);
void recreate(const std::filesystem::path& path);
void setCurrentWorkingPath(const std::filesystem::path& path);
std::vector<std::filesystem::path> files(const std::filesystem::path& path);
std::vector<std::filesystem::path> directories(const std::filesystem::path& path);
void sortFilesByNumber(std::vector<std::filesystem::path>& files);
bool isDir(const std::filesystem::path& path);
size_t fileSize(const std::filesystem::path& path);

//property tree
template <typename SCALAR>
struct PtreeSafeType {
  typedef int SAFE_SCALAR;
  static void put(tinyxml2::XMLElement& pt,const std::string& name,const SCALAR& val) {
    SAFE_SCALAR def=val;
    if(name.empty()) {
      std::string text=std::to_string(def);
      pt.SetText(text.c_str());
    } else pt.SetAttribute(name.c_str(),def);
  }
  static int get(const tinyxml2::XMLElement& pt,const std::string& name,const SCALAR& val) {
    SAFE_SCALAR def=val;
    if(name.empty())
      pt.QueryIntText(&def);
    else pt.QueryIntAttribute(name.c_str(),&def);
    return def;
  }
  static int get(const tinyxml2::XMLElement& pt,const std::string& name) {
    SAFE_SCALAR def=0;
    if(name.empty()) {
      assert(pt.QueryIntText(&def)==tinyxml2::XML_SUCCESS);
    } else {
      assert(pt.QueryIntAttribute(name.c_str(),&def)==tinyxml2::XML_SUCCESS);
    }
    return def;
  }
  static SAFE_SCALAR toSafe(SCALAR val) {
    return val;
  }
};
template <>
struct PtreeSafeType<float> {
  typedef float SAFE_SCALAR;
  //text
  static void put(tinyxml2::XMLElement& pt,const std::string& name,const float& val) {
    SAFE_SCALAR def=val;
    if(name.empty()) {
      std::string text=std::to_string(def);
      pt.SetText(text.c_str());
    } else pt.SetAttribute(name.c_str(),def);
  }
  static float get(const tinyxml2::XMLElement& pt,const std::string& name,const float& val) {
    SAFE_SCALAR def=val;
    if(name.empty())
      pt.QueryFloatText(&def);
    else pt.QueryFloatAttribute(name.c_str(),&def);
    return def;
  }
  static float get(const tinyxml2::XMLElement& pt,const std::string& name) {
    SAFE_SCALAR def=0;
    if(name.empty()) {
      assert(pt.QueryFloatText(&def)==tinyxml2::XML_SUCCESS);
    } else {
      assert(pt.QueryFloatAttribute(name.c_str(),&def)==tinyxml2::XML_SUCCESS);
    }
    return def;
  }
  static SAFE_SCALAR toSafe(float val) {
    return val;
  }
};
template <>
struct PtreeSafeType<double> {
  typedef float SAFE_SCALAR;
  static void put(tinyxml2::XMLElement& pt,const std::string& name,const double& val) {
    SAFE_SCALAR def=toSafe(val);
    if(name.empty()) {
      std::string text=std::to_string(def);
      pt.SetText(text.c_str());
    } else pt.SetAttribute(name.c_str(),def);
  }
  static double get(const tinyxml2::XMLElement& pt,const std::string& name,const double& val) {
    SAFE_SCALAR def=toSafe(val);
    if(name.empty())
      pt.QueryFloatText(&def);
    else pt.QueryFloatAttribute(name.c_str(),&def);
    return def;
  }
  static double get(const tinyxml2::XMLElement& pt,const std::string& name) {
    SAFE_SCALAR def=0;
    if(name.empty()) {
      assert(pt.QueryFloatText(&def)==tinyxml2::XML_SUCCESS);
    } else {
      assert(pt.QueryFloatAttribute(name.c_str(),&def)==tinyxml2::XML_SUCCESS);
    }
    return def;
  }
  static SAFE_SCALAR toSafe(double val) {
    return (SAFE_SCALAR)val;
  }
};
#ifdef QUADMATH_SUPPORT
template <>
struct PtreeSafeType<float128> {
  typedef float SAFE_SCALAR;
  static void put(tinyxml2::XMLElement& pt,const std::string& name,const float128& val) {
    SAFE_SCALAR def=val.convert_to<SAFE_SCALAR>();
    if(name.empty()) {
      std::string text=std::to_string(def);
      pt.SetText(text.c_str());
    } else pt.SetAttribute(name.c_str(),def);
  }
  static float128 get(const tinyxml2::XMLElement& pt,const std::string& name,const float128& val) {
    SAFE_SCALAR def=val.convert_to<SAFE_SCALAR>();
    if(name.empty())
      pt.QueryFloatText(&def);
    else pt.QueryFloatAttribute(name.c_str(),&def);
    return def;
  }
  static float128 get(const tinyxml2::XMLElement& pt,const std::string& name) {
    SAFE_SCALAR def=0;
    if(name.empty()) {
      assert(pt.QueryFloatText(&def)==tinyxml2::XML_SUCCESS);
    } else {
      assert(pt.QueryFloatAttribute(name.c_str(),&def)==tinyxml2::XML_SUCCESS);
    }
    return def;
  }
  static SAFE_SCALAR toSafe(float128 val) {
    return val.convert_to<SAFE_SCALAR>();
  }
};
#endif
template <>
struct PtreeSafeType<std::string> {
  typedef std::string SAFE_SCALAR;
  static void put(tinyxml2::XMLElement& pt,const std::string& name,const std::string& val) {
    if(name.empty()) {
      pt.SetText(val.c_str());
    } else pt.SetAttribute(name.c_str(),val.c_str());
  }
  static std::string get(const tinyxml2::XMLElement& pt,const std::string& name,const std::string& val) {
    if(name.empty())
      return pt.GetText();
    else {
      const char* ret=pt.Attribute(name.c_str());
      if(!ret)
        return val;
      else return ret;
    }
  }
  static std::string get(const tinyxml2::XMLElement& pt,const std::string& name) {
    const char* ret=NULL;
    if(name.empty())
      ret=pt.GetText();
    else ret=pt.Attribute(name.c_str());
    assert(ret!=0);
    return ret;
  }
  static SAFE_SCALAR toSafe(std::string val) {
    return val;
  }
};

//put/get
const tinyxml2::XMLElement* getChild(const tinyxml2::XMLElement& pt,const std::string& name);
tinyxml2::XMLElement* getChild(tinyxml2::XMLElement& pt,const std::string& name);
tinyxml2::XMLElement* addChild(tinyxml2::XMLElement& pt,const std::string& name);
const tinyxml2::XMLElement* getAttributeInfo(const tinyxml2::XMLElement& pt,std::string& name);
tinyxml2::XMLElement* getAttributeInfoPut(tinyxml2::XMLElement& pt,std::string& name);
bool hasAttribute(const tinyxml2::XMLElement& pt,const std::string& name);
template <typename T>
void put(tinyxml2::XMLElement& pt,const std::string& name,const T& val) {
  std::string nameProcessed=name;
  tinyxml2::XMLElement* e=getAttributeInfoPut(pt,nameProcessed);
  assert(e);
  PtreeSafeType<T>::put(*e,nameProcessed,val);
}
template <typename T>
void put(tinyxml2::XMLDocument& pt,const std::string& name,const T& val) {
  put<T>(*(pt.RootElement()),name,val);
}
template <typename T>
void putCond(tinyxml2::XMLElement& pt,const std::string& path,T val) {
  if(!hasAttribute(pt,path))
    put<T>(pt,path,val);
}
template <typename T>
void putCond(tinyxml2::XMLDocument& pt,const std::string& path,T val) {
  putCond<T>(*(pt.RootElement()),path,val);
}
template <typename T>
T get(const tinyxml2::XMLElement& pt,const std::string& name,const T& val) {
  std::string nameProcessed=name;
  const tinyxml2::XMLElement* e=getAttributeInfo(pt,nameProcessed);
  if(!e)
    return val;
  else return PtreeSafeType<T>::get(*e,nameProcessed,val);
}
template <typename T>
T get(const tinyxml2::XMLDocument& pt,const std::string& name,const T& val) {
  return get<T>(*(pt.RootElement()),name,val);
}
template <typename T>
T get(const tinyxml2::XMLElement& pt,const std::string& name) {
  std::string nameProcessed=name;
  const tinyxml2::XMLElement* e=getAttributeInfo(pt,nameProcessed);
  if(!e)
    throw std::runtime_error("Attribute "+name+" does not exist!");
  return PtreeSafeType<T>::get(*e,nameProcessed);
}
template <typename T>
T get(const tinyxml2::XMLDocument& pt,const std::string& name) {
  return get<T>(*(pt.RootElement()),name);
}
extern std::vector<std::string> getNames(std::string val,char sep=',');
extern std::vector<std::string> getNames(const tinyxml2::XMLDocument& pt,const std::string& name,char sep=',');

//parsePtree
std::vector<std::string> toParams(int argc,char** argv);
std::string parseProps(int argc,char** argv,tinyxml2::XMLElement& pt);
std::string parseProps(const std::vector<std::string>& params,tinyxml2::XMLElement& pt);
std::string parseProps(int argc,char** argv,tinyxml2::XMLDocument& pt);
std::string parseProps(const std::vector<std::string>& params,tinyxml2::XMLDocument& pt);
template <typename EIGEN_VEC>
EIGEN_VEC parsePtreeDef(const tinyxml2::XMLElement& pt,const std::string& path,const std::string& def,int r=-1,int c=-1) {
  EIGEN_VEC ret;
  if(r<0&&c<0)
    ret.setZero();
  else if(r>=0&&c<0)
    ret.setZero(r,ret.cols());
  else if(r<0&&c>=0)
    ret.setZero(ret.rows(),c);
  else ret.setZero(r,c);
  //read
  std::string str=get<std::string>(pt,path.c_str(),def.c_str());
  std::vector<std::string> strs=split(str," _,");
  int newSz=0;
  for(int j=0; j<(int)strs.size(); j++)
    if(!strs[j].empty())
      strs[newSz++]=strs[j];
  strs.resize(newSz);
  if(str.empty()) //fixing boost bug
    strs.clear();
  assert(ret.size()==0 || (int)strs.size() == ret.size());
  if(ret.size()==0)
    ret.resize((int)strs.size(),1);
  for(int i=0; i<ret.size(); i++) {
    typename PtreeSafeType<typename EIGEN_VEC::Scalar>::SAFE_SCALAR val;
    std::istringstream(strs[i]) >> val;
    ret.data()[i]=val;
  }
  return ret;
}
template <typename EIGEN_VEC>
EIGEN_VEC parsePtreeDef(const tinyxml2::XMLDocument& pt,const std::string& path,const std::string& def,int r=-1,int c=-1) {
  return parsePtreeDef<EIGEN_VEC>(*(pt.RootElement()),path,def,r,c);
}
template <typename EIGEN_VEC>
EIGEN_VEC parsePtreeDef(const tinyxml2::XMLElement& pt,const std::string& path,const EIGEN_VEC& ret,int r=-1,int c=-1) {
  std::string val;
  for(int i=0; i<ret.size(); i++)
    if(i == 0)
      val+=std::to_string(PtreeSafeType<typename EIGEN_VEC::Scalar>::toSafe(ret[i]));
    else val+=" "+std::to_string(PtreeSafeType<typename EIGEN_VEC::Scalar>::toSafe(ret[i]));
  return parsePtreeDef<EIGEN_VEC>(pt,path,val,r,c);
}
template <typename EIGEN_VEC>
EIGEN_VEC parsePtreeDef(const tinyxml2::XMLDocument& pt,const std::string& path,const EIGEN_VEC& ret,int r=-1,int c=-1) {
  return parsePtreeDef<EIGEN_VEC>(*(pt.RootElement()),path,ret,r,c);
}
template <typename EIGEN_VEC>
EIGEN_VEC parsePtreeDef(const tinyxml2::XMLElement& pt,const std::string& path,typename EIGEN_VEC::Scalar val,int r=-1,int c=-1) {
  EIGEN_VEC ret;
  if(r<0&&c<0)
    ret.setZero();
  else if(r>=0&&c<0)
    ret.setZero(r);
  else ret.setZero(r,c);
  ret.setConstant(val);
  return parsePtreeDef<EIGEN_VEC>(pt,path,ret,r,c);
}
template <typename EIGEN_VEC>
EIGEN_VEC parsePtreeDef(const tinyxml2::XMLDocument& pt,const std::string& path,typename EIGEN_VEC::Scalar val,int r=-1,int c=-1) {
  return parsePtreeDef<EIGEN_VEC>(*(pt.RootElement()),path,val,r,c);
}
template <typename EIGEN_VEC>
EIGEN_VEC parsePtree(const tinyxml2::XMLElement& pt,const std::string& path,int r=-1,int c=-1) {
  return parsePtreeDef<EIGEN_VEC>(pt,path,"",r,c);
}
template <typename EIGEN_VEC>
EIGEN_VEC parsePtree(const tinyxml2::XMLDocument& pt,const std::string& path,int r=-1,int c=-1) {
  return parsePtreeDef<EIGEN_VEC>(*(pt.RootElement()),path,"",r,c);
}
template <typename EIGEN_VEC>
void putPtree(tinyxml2::XMLElement& pt,const std::string& path,const EIGEN_VEC& ret) {
  std::string val;
  for(int i=0; i<ret.size(); i++) {
    if(i == 0)
      val+=std::to_string(PtreeSafeType<typename EIGEN_VEC::Scalar>::toSafe(ret.data()[i]));
    else val+=" "+std::to_string(PtreeSafeType<typename EIGEN_VEC::Scalar>::toSafe(ret.data()[i]));
  }
  put<std::string>(pt,path,val);
}
template <typename EIGEN_VEC>
void putPtree(tinyxml2::XMLDocument& pt,const std::string& path,const EIGEN_VEC& ret) {
  putPtree<EIGEN_VEC>(*(pt.RootElement()),path,ret);
}
template <typename EIGEN_VEC>
std::string toStringPtree(const EIGEN_VEC& v) {
  std::ostringstream oss;
  for(int i=0; i<v.size(); i++) {
    oss << v[i];
    if(i<v.size()-1)
      oss << ",";
  }
  return oss.str();
}
template <typename T>
std::string toStringPtree(const std::vector<T>& v) {
  std::ostringstream oss;
  for(int i=0; i<(int)v.size(); i++) {
    oss << v[i];
    if(i<(int)v.size()-1)
      oss << ",";
  }
  return oss.str();
}
}

#endif
