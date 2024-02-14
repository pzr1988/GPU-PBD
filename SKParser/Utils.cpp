#include "Utils.h"
#include <iostream>
#include <iomanip>

namespace PHYSICSMOTION {
std::vector<std::string> split(const std::string& l,const std::string& sep) {
  int i,last=-1;
  std::vector<std::string> ret;
  for(i=0; i<(int)l.size(); i++)
    if(std::find(sep.begin(),sep.end(),l[i])!=sep.end()) {
      if(last<i-1)
        ret.push_back(l.substr(last+1,i-last-1));
      last=i;
    }
  if(last<i-1)
    ret.push_back(l.substr(last+1,i-last-1));
  return ret;
}
bool beginsWith(const std::string& l,const std::string& s) {
  return l.length()>=s.length() && l.substr(0,s.length())==s;
}
bool endsWith(const std::string& l,const std::string& s) {
  return l.length()>=s.length() && l.substr(l.size()-s.length(),s.length())==s;
}
std::string toUpper(const std::string& l) {
  std::string ret=l;
  std::for_each(ret.begin(),ret.end(),[&](char c) {
    return std::toupper(c);
  });
  return ret;
}

//filesystem
bool notDigit(char c) {
  return !std::isdigit(c);
}
bool lessDirByNumber(std::filesystem::path A,std::filesystem::path B) {
  std::string strA=A.string();
  std::string::iterator itA=std::remove_if(strA.begin(),strA.end(),notDigit);
  strA.erase(itA,strA.end());

  std::string strB=B.string();
  std::string::iterator itB=std::remove_if(strB.begin(),strB.end(),notDigit);
  strB.erase(itB,strB.end());

  int idA,idB;
  std::istringstream(strA) >> idA;
  std::istringstream(strB) >> idB;
  return idA < idB;
}
bool exists(const std::filesystem::path& path) {
  return std::filesystem::exists(path);
}
void removeDir(const std::filesystem::path& path) {
  if(std::filesystem::exists(path))
    try {
      std::filesystem::remove_all(path);
    } catch(...) {}
}
void create(const std::filesystem::path& path) {
  if(!std::filesystem::exists(path))
    try {
      std::filesystem::create_directory(path);
    } catch(...) {}
}
void recreate(const std::filesystem::path& path) {
  removeDir(path);
  try {
    std::filesystem::create_directory(path);
  } catch(...) {}
}
void setCurrentWorkingPath(const std::filesystem::path& path) {
  if(isDir(path))
    std::filesystem::current_path(path);
  else std::filesystem::current_path(path.parent_path());
}
std::vector<std::filesystem::path> files(const std::filesystem::path& path) {
  std::vector<std::filesystem::path> ret;
  for(std::filesystem::directory_iterator beg(path),end; beg!=end; beg++)
    if(std::filesystem::is_regular_file(*beg))
      ret.push_back(*beg);
  return ret;
}
std::vector<std::filesystem::path> directories(const std::filesystem::path& path) {
  std::vector<std::filesystem::path> ret;
  for(std::filesystem::directory_iterator beg(path),end; beg!=end; beg++)
    if(std::filesystem::is_directory(*beg))
      ret.push_back(*beg);
  return ret;
}
void sortFilesByNumber(std::vector<std::filesystem::path>& files) {
  sort(files.begin(),files.end(),lessDirByNumber);
}
bool isDir(const std::filesystem::path& path) {
  return std::filesystem::is_directory(path);
}
size_t fileSize(const std::filesystem::path& path) {
  return (size_t)std::filesystem::file_size(path);
}

//put/get
const tinyxml2::XMLElement* getChild(const tinyxml2::XMLElement& pt,const std::string& name) {
  for(const tinyxml2::XMLElement* v=pt.FirstChildElement(); v; v=v->NextSiblingElement())
    if(v->Name()==name)
      return v;
  return NULL;
}
tinyxml2::XMLElement* getChild(tinyxml2::XMLElement& pt,const std::string& name) {
  for(tinyxml2::XMLElement* v=pt.FirstChildElement(); v; v=v->NextSiblingElement())
    if(v->Name()==name)
      return v;
  return NULL;
}
tinyxml2::XMLElement* addChild(tinyxml2::XMLElement& pt,const std::string& name) {
  //for(tinyxml2::XMLElement* v=pt.FirstChildElement(); v; v=v->NextSiblingElement())
  //  if(v->Name()==name)
  //    return v;
  tinyxml2::XMLElement* node=pt.GetDocument()->NewElement(name.c_str());
  pt.InsertEndChild(node);
  return node;
}
const tinyxml2::XMLElement* getAttributeInfo(const tinyxml2::XMLElement& pt,std::string& name) {
  const tinyxml2::XMLElement* curr=&pt;
  std::vector<std::string> paths=split(name,".");
  for(int i=0; i<(int)paths.size(); i++)
    if(paths[i]=="<xmlattr>") {
      assert(i==(int)paths.size()-2);
      name=paths[i+1];
      return curr;
    } else {
      curr=getChild(*curr,paths[i]);
      if(!curr)
        return curr;
    }
  name="";
  return curr;
}
tinyxml2::XMLElement* getAttributeInfoPut(tinyxml2::XMLElement& pt,std::string& name) {
  tinyxml2::XMLElement* curr=&pt;
  std::vector<std::string> paths=split(name,".");
  for(int i=0; i<(int)paths.size(); i++)
    if(paths[i]=="<xmlattr>") {
      assert(i==(int)paths.size()-2);
      name=paths[i+1];
      return curr;
    } else {
      tinyxml2::XMLElement* tmp=getChild(*curr,paths[i]);
      if(!tmp) {
        tinyxml2::XMLElement* node=pt.GetDocument()->NewElement(paths[i].c_str());
        curr->InsertEndChild(node);
        curr=node;
      } else {
        curr=tmp;
      }
    }
  name="";
  return curr;
}
bool hasAttribute(const tinyxml2::XMLElement& pt,const std::string& name) {
  std::string nameProcessed=name;
  const tinyxml2::XMLElement* e=getAttributeInfo(pt,nameProcessed);
  if(!e)
    return false;
  else if(nameProcessed.empty())
    return true;
  else return e->Attribute(nameProcessed.c_str())!=NULL;
}
std::vector<std::string> getNames(std::string val,char sep) {
  std::vector<std::string> names;
  while(true) {
    size_t pos=val.find_first_of(sep);
    if(pos==std::string::npos)
      break;
    names.push_back(val.substr(0,pos));
    if(pos==val.length()-1) {
      val="";
      break;
    } else val=val.substr(pos+1);
  }
  if(!val.empty())
    names.push_back(val);
  return names;
}
std::vector<std::string> getNames(const tinyxml2::XMLDocument& pt,const std::string& name,char sep) {
  std::string str=get<std::string>(pt,name,"");
  return getNames(str,sep);
}

//parsePtree
std::vector<std::string> toParams(int argc,char** argv) {
  std::vector<std::string> params;
  for(int i=0; i<argc; i++)
    params.push_back(argv[i]);
  return params;
}
std::string parseProps(int argc,char** argv,tinyxml2::XMLElement& pt) {
  return parseProps(toParams(argc,argv),pt);
}
std::string parseProps(const std::vector<std::string>& params,tinyxml2::XMLElement& pt) {
  std::string addParam;
  for(int i=0; i<(int)params.size(); i++) {
    const std::string& str=params[i];
    size_t pos=str.find("=");
    if(pos != std::string::npos) {
      std::string LHS=str.substr(0,pos);
      std::string RHS=str.substr(pos+1);
      put<std::string>(pt,LHS.c_str(),RHS.c_str());
      addParam+="_"+str;
    }
  }
  return addParam+"_";
}
std::string parseProps(int argc,char** argv,tinyxml2::XMLDocument& pt) {
  return parseProps(argc,argv,*(pt.RootElement()));
}
std::string parseProps(const std::vector<std::string>& params,tinyxml2::XMLDocument& pt) {
  return parseProps(params,*(pt.RootElement()));
}
}
