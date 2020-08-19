#ifndef PDETOOLS_H
#define PDETOOLS_H
#include <string>
#include <fstream>

void stepIC(float* u,int Nx,int Ny,int Nz);
void constIC(float* u,float u0,int N);
void varIC(float* u,std::string filename,int N);
float squareWave(float t,float T,float ubsl,float uhigh,float ulow);
void print(float* u,int N,std::string filename);
void print(float* u,int Nx,int Ny,int Nz,int dim,int slc,std::string filename);

class print_scheduler
{
public:
  print_scheduler(float print_frequency);
  bool operator()(float t);
  void schedule(float start_time,float print_frequency);
private:
  float next_time;
  float* end_time;
  float* frequency;
  int count,index;
};

class time_writer
{
public:
  time_writer(std::string filename);
  ~time_writer();
  void operator()(double t)
  {
    out << t << std::endl;
  }
private:
  std::ofstream out;
};

class model 
{
public:
  model(float alpha,float beta,float gamma,float ub,float km,float lambda,float sigma)
    : alpha(alpha),beta(beta),gamma(gamma),ub(ub),km(km),lambda(lambda),sigma(sigma){}  
  model(float alpha,float beta,float gamma,float ub,float km)
    : model(alpha,beta,gamma,ub,km,1.0f,1.0f){}  
  model(float alpha,float beta,float ub,float km) 
    : model(alpha,beta,1.0f,ub,km,1.0f,1.0f){}
  
  float alpha,beta,gamma,ub,km,lambda,sigma;
};

class grid 
{
public:
  grid(int Nx,int Ny,int Nz,float dt,float ay,float az);
  int Nx,Ny,Nz,N;
  float dt,dx,dy,dz;
};

class geometry
{
public:
  geometry(float L,float H,float W,float w,float h);
  geometry(float L,float H,float W,float w,float h,float xs,float ys);
  float L,H,W,l,h,xs,ys;
};

#endif
