// Solve O2 Transport on Tissue Volume Using GPGPU Acceleration
#include <iostream>
#include <string>
#include <cuda.h>
#include "PDETools.h"
#include "PDEsolve.h"
#include "configuration.h"
using namespace std;

int main()
{
  // Output Filename (user input)
  string filename_out = "data/square_wave";
   
  // Initialize Physical Constants (user input)
  float D = 2.41e-5f; // diffusivity <cm^2/s>
  float k = 3.89e-5f; // solubility <mLO2/mL/mmHg>
  float VO2 = 1.50e-4f; // O2 consumption <mLO2/mL/s>
  float Pcrit = 0.5f; // critical PO2 <mmHg>
  float P0 = 48.0f; // capillary PO2 <mmHg>
  float K = 30.0f; // capillary source <mmHg/s>
  float L = 0.2f; // tissue length <cm>
  float H = 0.2f; // tissue height <cm>
  float W = 0.06f; // tissue depth <cm>
  
  // Square Wave Parameters
  float Pbsl = 38.0f; // Baseline PO2 <mmHg>
  float Phigh = 53.2f; // High PO2 <mmHg>
  float Plow = 15.2f; // Low PO2 <mmHg>
  
  // Simulation Time (user input)
  float sim_time = 360.0f; // simulation time <s>
  float print_frequency = 0.25f; // print frequency <s>
  
  // Initialize Computational Domain (user input)
  int Nx = 144;
  int Ny = 144;
  int Nz = 48;
  float dt = 1e-6;
  
  // Calculate Dimensionless Parameters
  float q = VO2*L*L/(P0*D*k);
  float ucrit = Pcrit/P0;
  float s = K*L*L/(P0*D);
  float tau = L*L/D;
  float ay = H/L;
  float az = W/L;
  
  // Calculate Computational Parameters
  int N = Nx*Ny*Nz;
  size_t size = N*sizeof(float);
  float dx = 1.0f/(Nx-1.0f);
  float dy = ay/(Ny-1.0f);
  float dz = az/(Nz-1.0f);
  float pfreq = print_frequency/tau; 
  float T = sim_time/tau;
  
  // Print Parameters
  float ds = min(min(dx,dy),dz); 
  cout << "\n\nSimulation Parameters\n\n";
  cout << "dt/dx^2 = " << dt/ds/ds << endl;
  cout << "Runtime: " << T << endl << endl;
  cout << "tau = " << tau << endl;
  cout << "s = " << s << endl;
  cout << "q = " << q << endl << endl;
  cout << "ucrit = " << ucrit << endl;
  cout << "uhigh = " << Phigh/P0 << endl;
  cout << "ulow = " << Plow/P0 << endl << endl;
  
  // Allocate Memory on Host
  float* u_h = new float[N]();
  //constIC(u_h,0.0f,Nx,Ny,Nz);
  varIC(u_h,"data/baseline_steady-state1.csv",Nx*Ny*Nz);
  print(u_h,N,filename_out,0);
  
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
  float t = 0.0f;
  int np = 1; 
  float uwin;
  for (int nt = 1; t < T; nt++)
  { 
    // Boundary Condition
    //uwin = squareWave(t,T,Pbsl/P0,Phigh/P0,Plow/P0); // square wave in time
    uwin = Plow/P0; // constant in time
   
    // Call GPU Kernel
    step<<<dimGrid,dimBlock>>>(uold_d,unew_d,q,ucrit,s,uwin,Nx,Ny,Nz,dx,dy,dz,dt);
    t += dt;
    
    // Print Solution
    if (t >= np*pfreq)
    {
      cout << "Writing t = " << t << "...\n";
      cudaMemcpy(u_h,unew_d,size,cudaMemcpyDeviceToHost);
      print(u_h,N,filename_out,np);
      np++;
    }
  }
  
  // Delete Pointers
  delete[] u_h;
  cudaFree(uold_d);
  cudaFree(unew_d);
  
  return 0;
}
