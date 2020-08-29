#pragma once 

#include <string>
#include <fstream>

void stepIC(float* u,int Nx,int Ny,int Nz);
void constIC(float* u,float u0,int size);
void varIC(float* u,std::string filename,int size);
float squareWave(float t,float endTime,float ubsl,float uhigh,float ulow);
void print(float* u,int size,std::string filename);
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

struct model 
{
  model(float alpha,float beta,float gamma,float ub,float km,float lambda,float sigma)
    : alpha(alpha),beta(beta),gamma(gamma),ub(ub),km(km),lambda(lambda),sigma(sigma){}  
  model(float alpha,float beta,float gamma,float ub,float km)
    : model(alpha,beta,gamma,ub,km,1.0f,1.0f){}  
  model(float alpha,float beta,float ub,float km) 
    : model(alpha,beta,1.0f,ub,km,1.0f,1.0f){}
  
  float alpha,beta,gamma,ub,km,lambda,sigma;
};

struct geometry 
{
  geometry(float L,float H,float W,float l,float h,float xs,float ys,float th)
    : L(L),H(H),W(W),l(l),h(h),xs(xs),ys(ys),th(th) {}
  geometry(float L,float H,float W,float l,float h,float xs,float ys)
    : geometry(L,H,W,l,h,xs,ys,-1.0) {}
  geometry(float L,float H,float W,float l,float h)
    : geometry(L,H,W,l,h,0.0,0.0,-1.0) {}
      
  float L,H,W,l,h,xs,ys,th;
};

struct grid 
{
  grid(int Nx,int Ny,int Nz,float dx,float dy,float dz,float dt);
  int Nx,Ny,Nz,size;
  float dt,dx,dy,dz;
};