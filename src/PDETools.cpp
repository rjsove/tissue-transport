#include "PDETools.h"
#include <string>
#include <fstream>
#include <cmath>
#include <iostream>

using namespace std;

// Initial Conditions
// Cube centered in domain with side length of 1/4 set to 1
void stepIC(float* u,int Nx,int Ny,int Nz)
{
  for (int i = 3*Nx/8; i <= 5*Nx/8; i++)
  {
    for (int j = 3*Ny/8; j <= 5*Ny/8; j++)
    {
      for (int k = 3*Nz/8; k <= 5*Nz/8; k++)
      {
        u[i+Nx*(j+k*Ny)] = 1.0f;
      }
    }
  }
}

// Constant Value of u0 
void constIC(float* u,float u0,int Nx,int Ny,int Nz)
{
  for (int i = 0; i < Nx; i++)
  {
    for (int j = 0; j < Ny; j++)
    {
      for (int k = 0; k < Nz; k++)
      {
        u[i+Nx*(j+k*Ny)] = u0;
      }
    }
  }
}

// Variable I.C. set from file
void varIC(float* u,std::string filename,int N)
{
  ifstream file;
  file.open(filename);
  int i = 0;
  for (int i = 0; i < N; i++)
  {
    file >> u[i];
  }
  file.close();
}

// Square Wave Boundary Condition
float squareWave(float t,float T,float ubsl,float uhigh,float ulow)
{
  // Window PO2
  float u;
  
  // Baseline
  if (t < T/6 || t >= 5*T/6)
  {
    u = ubsl;
  }
  // High O2
  else if (t < T/2)
  {
    u = uhigh;
  }
  // Low O2
  else
  {
    u = ulow;
  }
  
  return u;
}

// I/O
void print(float* u,int N,string filename)
{
  ofstream file;
  file.open(filename);
  for (int i = 0; i < N; i++)
  {
    file << u[i] << endl;
  }
  file.close();
}

print_scheduler::print_scheduler(float print_frequency)
{
  // Set Schedule Count
  count = 1;
  
  // Intialize Schedules
  end_time = new float[100];
  end_time[0] = nanf("");
  frequency = new float[100];
  frequency[0] = print_frequency;
  
  // Set Next Time 
  next_time = print_frequency;
  index = 0;
}

bool print_scheduler::operator()(float t)
{
  bool output = false;
  if (t >= next_time)
  {
    output = true;
    if (next_time+frequency[index] > end_time[index])
    {
      next_time = end_time[index]+frequency[index+1];
      index++;
    }
    else
    {
      next_time += frequency[index];
    }
  }
  return output;
}

void print_scheduler::schedule(float start_time,float print_frequency)
{
  end_time[count-1] = start_time;
  end_time[count] = nanf("");
  frequency[count] = print_frequency;
  count++;
}

time_writer::time_writer(string filename)
{
  out.open(filename);
}

time_writer::~time_writer()
{
  out.close();
}

grid::grid(int Nx,int Ny,int Nz,float dt,float ay,float az)
{
  this->Nx = Nx;
  this->Ny = Ny;
  this->Nz = Nz;
  this->N = Nx*Ny*Nz;
  this->dt = dt;
  this->dx = 1.0f/(Nx-1);
  this->dy = ay/(Ny-1);
  this->dz = az/(Nz-1);
}

geometry::geometry(float L,float H,float W,float l,float h)
{
  this->L = L;
  this->H = H;
  this->W = W;
  this->l = l;
  this->h = h;
  this->xs = 0.0f;
  this->ys = 0.0f;
}

geometry::geometry(float L,float H,float W,float l,float h,float xs,float ys)
{
  this->L = L;
  this->H = H;
  this->W = W;
  this->l = l;
  this->h = h;
  this->xs = xs;
  this->ys = ys;
}