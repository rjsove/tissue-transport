#pragma once 

#include <string>

#define NUMBER_OF_PARAMETERS 35

namespace ttgpu
{  
struct istream
{
  istream();
  struct
  {
    float D,k,Dpdms,kpdms,M0,Km,Pc,K;
  } physics;
  struct
  {
    float L,H,W,l,h,xs,ys,th;
  } geometry;
  struct
  {
    bool homogeneousIC,constantBC;
    std::string icFilename;
    float icValue,bcValue,bcSqWvHighValue,bcSqWvLowValue;
  }ibvp;
  struct 
  {
    int Nx,Ny,Nz;
    float dt;
  } discretization;
  struct 
  {
    bool writeFullSolution;
    int numberOfSlices,*sliceDim,*slice;
    float sim_time,write_interval;
    std::string name,*fileSuffix;
  } output; 
private:
  void read(std::string filename);
  void setParameters(std::string line);
};
}