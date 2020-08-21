// Solve O2 Transport on Tissue Volume Using GPGPU Acceleration
#include <iostream>
#include <string>
#include <cuda.h>
#include "PDETools.h"
#include "PDEsolve.h"
#include "configuration.h"
#include "timer.h"
using namespace std;

int main(int argc,char** argv)
{
  // Start Main Timer 
  timer timer1("Total");
  
  // Square Wave Parameters
  float Pbsl = 38.0f; // Baseline PO2 <mmHg>
  float Phigh = 53.2f; // High PO2 <mmHg>
  float Plow = 15.2f; // Low PO2 <mmHg>
   
  // Initialize Physical Constants (user input)
  float D = 2.41e-5f; // diffusivity <cm^2/s>
  float k = 3.89e-5f; // solubility <mLO2/mL/mmHg>
  float Dpdms = 3.55e-5f; // diffusivity <cm^2/s>
  float kpdms = 1.32e-5f; // solubility <mLO2/mL/mmHg>
  float VO2 = 1.50e-4f; // O2 consumption <mLO2/mL/s>
  float Pcrit = 0.5f; // critical PO2 <mmHg>
  float P0 = 48.0f; // capillary PO2 <mmHg>
  float K = 30.0f; // capillary source <mmHg/s>
  float L = 0.2f; // tissue length <cm>
  float H = 0.2f; // tissue height <cm>
  float W = 0.06f; // tissue depth <cm>
  float l = 0.04f; // window length <cm>
  float h = 0.02f; // window height <cm>
  //float th = 0.004; // PDMS layer thickness <cm>
  
  // Simulation Time (user input)
  float sim_time = 750.0f; // simulation time <s> was 360.0f
  if (argc == 2)
    sim_time = atof(argv[1]);
  float print_frequency = 750.0f; // print frequency <s> was 0.25f
  
  // Write-Out Schedule
  // 0-10s: 1s, 10-30s: 5s, 30-180s: 30s
  print_scheduler print_time(print_frequency);
  //print_time.schedule(10.0f,5.0f); // (start_time <s>,frequency <s>) 
  //print_time.schedule(30.0f,30.0f);  
  
  // Initialize Computational Domain (user input)
  int Nx = 288; // 576, 288, 144
  int Ny = 288; // 576, 288, 144
  int Nz = 96; // 192, 96, 48 
  float dt = 1e-6; // was 1e-6
  
  // Output Filename (user input)
  string dir = "out/PO2/";
  string filename = "test";
  int dim0 = 1;
  int slc0 = Ny/2-1;
  int dim1 = 2;
  int slc1 = 17;
  int dim2 = 2;
  int slc2 = 22;
  
  // Calculate Dimensionless Parameters
  float tau = L*L/D;
  float alpha = K/P0*tau;
  float ub = 1;
  float beta = VO2/(P0*k)*tau; 
  float km = Pcrit/P0;
  float lambda = Dpdms/D;
  float sigma = lambda*kpdms/k;
  float ay = H/L;
  float az = W/L;
  
  // Calculate Computational Parameters
  model mdl(alpha,beta,1.0f,ub,km,lambda,sigma);
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
  cout << "\n\n---Oxygen---\n\n";
  cout << "\n\nSimulation Parameters\n\n";
  cout << "Runtime: " << sim_time << " s \t\t[" << T/dt << " time steps]\n\n";
  cout << "dt/dx^2 = " << dt/ds/ds << endl << endl;
  cout << "tau = " << tau << " s\n";
  cout << "P0 = " << P0 << " mmHg\n\n";
  cout << "alpha = " << alpha << endl;
  cout << "ub = " << ub << endl;
  cout << "beta = " << beta << endl;
  cout << "km = " << km << endl << endl;
  cout << "ubsl = " << Pbsl/P0 << endl;
  cout << "uhigh = " << Phigh/P0 << endl;
  cout << "ulow = " << Plow/P0 << endl << endl;
  
  // Allocate Memory on Host
  float* u_h = new float[N]();
  constIC(u_h,1.0f,N);
  //varIC(u_h,"data/baseline_steady-state1.csv",N);
  print(u_h,Nx,Ny,Nz,dim0,slc0,dir+filename+"0.csv");
  
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
  float uwin;
  for (int nt = 1; t < T; nt++)
  { 
    // Boundary Condition
    //uwin = squareWave(t,T,Pbsl/P0,Phigh/P0,Plow/P0); // square wave in time
    uwin = Plow/P0; // constant in time
   
    // Call GPU Kernel
    step<<<dimGrid,dimBlock>>>(uold_d,unew_d,uwin,mdl,grd,geo);
    t += dt;
    
    // Print Solution
    if (print_time(t*tau))
    {
      timer timer2("Write Out");
      cout << "Writing t = " << t << "...\n";
      cudaMemcpy(u_h,unew_d,size,cudaMemcpyDeviceToHost);
      print(u_h,Nx,Ny,Nz,dim0,slc0,dir+filename+"_1_"+to_string(np)+".csv");
      print(u_h,Nx,Ny,Nz,dim1,slc1,dir+filename+"_2_"+to_string(np)+".csv");
      print(u_h,Nx,Ny,Nz,dim2,slc2,dir+filename+"_3_"+to_string(np)+".csv");
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
