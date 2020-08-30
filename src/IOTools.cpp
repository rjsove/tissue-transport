#include "IOTools.h"
#include <iostream>
#include <fstream>
#include <sstream>

void ttgpu::istream::setParameters(std::string line)
{
  std::string name;
  size_t pos=0;
  while ((pos = line.find('=', pos)) != std::string::npos)
      line[pos] = ' ';
  std::istringstream ss(line);
  ss >> name;
  
  // Store Value
  // Physics
  if (name.compare("D")==0){
    ss >> this->physics.D; return;}
  else if (name.compare("k")==0){
    ss >> this->physics.k; return;}
  else if (name.compare("Dpdsm")==0){
    ss >> this->physics.Dpdms; return;}
  else if (name.compare("kpdms")==0){
    ss >> this->physics.kpdms; return;}
  else if (name.compare("M0")==0){
    ss >> this->physics.M0; return;}
  else if (name.compare("Km")==0){
    ss >> this->physics.Km; return;}
  else if (name.compare("Pc")==0){
    ss >> this->physics.Pc; return;}
  else if (name.compare("K")==0){
    ss >> this->physics.K; return;}
  // Geometry 
  else if (name.compare("L")==0){
    ss >> this->geometry.L; return;}
  else if (name.compare("H")==0){
    ss >> this->geometry.H; return;}
  else if (name.compare("W")==0){
    ss >> this->geometry.W; return;}
  else if (name.compare("l")==0){
    ss >> this->geometry.l; return;}
  else if (name.compare("h")==0){
    ss >> this->geometry.h; return;}
  else if (name.compare("xs")==0){
    ss >> this->geometry.xs; return;}
  else if (name.compare("ys")==0){
    ss >> this->geometry.ys; return;}
  else if (name.compare("th")==0){
    ss >> this->geometry.th; return;}
  // Boundary Conditions 
  else if (name.compare("homogeneousIC")==0){
    ss >> this->ibvp.homogeneousIC; return;}
  else if (name.compare("constantBC")==0){
    ss >> this->ibvp.constantBC; return;}
  else if (name.compare("icFilename")==0){
    ss >> this->ibvp.icFilename; return;}
  else if (name.compare("icValue")==0){
    ss >> this->ibvp.icValue; return;}
  else if (name.compare("bcValue")==0){
    ss >> this->ibvp.bcValue; return;}
  else if (name.compare("bcSqWvHighValue")==0){
    ss >> this->ibvp.bcSqWvHighValue; return;}
  else if (name.compare("bcSqWvLowValue")==0){
    ss >> this->ibvp.bcSqWvLowValue; return;}
  // Discretization  
  else if (name.compare("Nx")==0)
    ss >> this->discretization.Nx;
  else if (name.compare("Ny")==0)
    ss >> this->discretization.Ny;
  else if (name.compare("Nz")==0)
    ss >> this->discretization.Nz;
  else if (name.compare("dt")==0)
    ss >> this->discretization.dt;  
  // Output 
       
}

void ttgpu::istream::read(std::string filename)
{
  std::string line;
  std::ifstream file(filename);
  while (getline(file,line))
    setParameters(line);
  file.close();
}

ttgpu::istream::istream()
{
  // Read Data from Files 
  read("physics.ttgpu");
  read("geometry.ttgpu");
  read("boundary_conditions.ttgpu");
  read("discretization.ttgpu");
  read("output_settings.ttgpu");
}