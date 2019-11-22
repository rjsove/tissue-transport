#ifndef PDETOOLS_H
#define PDETOOLS_H
#include <string>
#include <fstream>

void stepIC(float* u,int Nx,int Ny,int Nz);
void constIC(float* u,float u0,int Nx,int Ny,int Nz);
void varIC(float* u,std::string filename,int N);
float squareWave(float t,float T,float ubsl,float uhigh,float ulow);
void print(float* u,int N,std::string filename,int nt);

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

#endif
