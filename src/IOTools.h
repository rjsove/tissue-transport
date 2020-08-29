#pragma once 

class std::string;

namespace ttgpu
{

enum streamType {PHYSICS,GEOMETRY,BOUNDARY_CONDITIONS,DISCRETIZATION,OUTPUT_SETTINGS};

struct inStream
{
  static std::string usr_directory; 
  virtual void readFile(){} 
};

struct physics: public inStream
{
  float D,k,Dpdms,kpdms,M0,Km,Pc,K;
  void readFile() override;
};

struct geometry: public inStream
{
  float L,H,W,l,h,xs,ys,th;
  void readFile() override;
};

struct boundaryConditions: public inStream 
{
  bool homogeneousIC,constantBC;
  std::string icFilename;
  float icValue,bcValue,bcSqWvHighValue,bcSqWvLowValue;
  void readFile() override;
};

struct discretization: public inStream
{
  int Nx,Ny,Nz;
  float dt;
  void readFile() override;
};

struct outputSettings: public inStream
{
  bool writeFullSolution;
  int numberOfSlices;
  int *sliceDim,*slice;
  std::string filename,*fileSuffix;
  void readFile() override;
};

inStream read(std::string filename,streamType type);
};