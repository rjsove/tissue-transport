// Solve Ach Transport on Tissue Volume Using GPGPU Acceleration
#include <iostream>
#include <string>
#include <cuda.h>
#include "PDETools.h"
#include "PDEsolve.h"
#include "configuration.h"
using namespace std;

int main(int argc,char** argv)
{
  // Boundary Condition (user input)
  float C0 = 1e-8; // <M>
  
  // Output Filename (user input)
  string dir = "out/Ach/";
  string filename = "test";
   
  // Initialize Physical Constants (user input)
  float D = 4.00e-6f; // diffusivity <cm^2/s>
  float q = 1.037e-1f; // blood vessel permeability <1/s>
  float Cb = 0.0f; // blood concentration <M>
  float Vmax = 4.83e-8f; // maximal consumption <M/s>
  float Km = 6.0e-3f; // Michaelis-Menton Constant <M>
  float L = 0.2f; // tissue length <cm>
  float H = 0.2f; // tissue height <cm>
  float W = 0.06f; // tissue depth <cm>
  float l = 0.06f; // window length <cm>
  float h = 0.03f; // window width <cm>
  
  // Simulation Time (user input)
  float sim_time = 181.0f; // simulation time <s>
  if (argc == 2)
    sim_time = atof(argv[1]);
  float print_frequency = 1.0f; // print frequency <s>
  
  // Write-Out Schedule
  // 0-10s: 1s, 10-30s: 5s, 30-180s: 30s
  print_scheduler print_time(print_frequency);
  print_time.schedule(10.0f,5.0f); // (start_time <s>,frequency <s>)
  print_time.schedule(30.0f,30.0f);
  
  // Initialize Computational Domain (user input)
  int Nx = 144;
  int Ny = 144;
  int Nz = 48;
  float dt = 1e-6;
  
  // Calculate Dimensionless Parameters
  float tau = L*L/D;
  float alpha = q*tau;
  float ub = Cb/C0;
  float beta = Vmax/C0*tau; 
  float km = Km/C0;
  float ay = H/L;
  float az = W/L;
  
  // Calculate Steady-State Solution
  float u0 = 0.0f;
  
  // Calculate Computational Parameters
  model mdl(alpha,beta,ub,km);
  grid grd(Nx,Ny,Nz,dt,ay,az);
  geometry geo(L,W,H,l,h);
  int N = Nx*Ny*Nz;
  size_t size = N*sizeof(float);
  float dx = 1.0f/(Nx-1.0f);
  float dy = ay/(Ny-1.0f);
  float dz = az/(Nz-1.0f);
  float T = sim_time/tau;
  
  // Print Parameters
  float ds = min(min(dx,dy),dz); 
  cout << "\n\n---Acetylcholine---\n\n";
  cout << "\n\nSimulation Parameters\n\n";
  cout << "Runtime: " << sim_time << " s \t\t[" << T/dt << " time steps]\n\n";
  cout << "dt/dx^2 = " << dt/ds/ds << endl << endl;
  cout << "tau = " << tau << " s\n";
  cout << "C0 = " << C0 << " M\n\n";
  cout << "u0 = " << u0 << endl;
  cout << "alpha = " << alpha << endl;
  cout << "ub = " << ub << endl;
  cout << "beta = " << beta << endl;
  cout << "km = " << km << endl << endl;
  
  // Allocate Memory on Host
  float* u_h = new float[N]();
  constIC(u_h,u0,Nx,Ny,Nz); 
  print(u_h,N,dir+filename+"0.csv");
  
  // Allocate Memory on Device 
  float *uold_d,*unew_d;
  cudaMalloc((void**) &uold_d,size);
  cudaMalloc((void**) &unew_d,size);
    
  // Copy Memory to Device
  cudaMemcpy(uold_d,u_h,size,cudaMemcpyHostToDevice);
  
  // Setup Configuration for GPU
  dim3 dimGrid(GRID_SIZE_X,GRID_SIZE_Y,GRID_SIZE_Z);
  dim3 dimBlock(BLOCK_SIZE_X,BLOCK_SIZE_Y,BLOCK_SIZE_Z);

  // Time Iteration
  float t = 0.0f; int np = 1; time_writer write_time(dir+"t.csv"); write_time(t*tau);
  float uwin = 1.0f; // Boundary Condition
  for (int nt = 1; t < T; nt++)
  { 
    // Call GPU Kernel
    step<<<dimGrid,dimBlock>>>(uold_d,unew_d,uwin,mdl,grd,geo);
    t += dt;
    
    // Print Solution
    if (print_time(t*tau))
    {
      cout << "Writing t = " << t*tau << " s...\n";
      cudaMemcpy(u_h,unew_d,size,cudaMemcpyDeviceToHost);
      print(u_h,N,dir+filename+to_string(np)+".csv");
      write_time(t*tau);
      np++;
    }
  }
  
  // Delete Pointers
  delete[] u_h;
  cudaFree(uold_d);
  cudaFree(unew_d);
  
  return 0;
}
