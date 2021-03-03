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
  float Pbsl = 38.0f; // Baseline PO2 <mmHg> 5%
  float Phigh = 91.2f; // High PO2 <mmHg> 12% 
  float Plow = 15.2f; // Low PO2 <mmHg> 2%
   
  // Initialize Physical Constants (user input)
  float D = 2.41e-5f; // diffusivity <cm^2/s>
  float k = 3.89e-5f; // solubility <mLO2/mL/mmHg>
  float Dpdms = 3.40e-5f; // diffusivity <cm^2/s> [Merkel 2000] 3.55e-5
  float kpdms = 3.68e-5f; // solubility <mLO2/mL/mmHg> [Shiku 2006]
  float VO2 = 1.50e-4f; // O2 consumption <mLO2/mL/s>
  float Pcrit = 0.5f; // critical PO2 <mmHg>
  float P0 = 48.0f; // capillary PO2 <mmHg>
  float K = 30.0f; // capillary source <mmHg/s>
  float L = 0.2f; // tissue length <cm> regular simulation: [0.2x0.2]
  float H = 0.2f; // tissue height <cm> five window simulation [0.52x0.44]
  float W = 0.06f; // tissue depth <cm>
  float l = 0.04f; // window length <cm>
  float h = 0.02f; // window height <cm>
  float xs = 0.1f; // horizontal window spacing <cm>
  float ys = 0.2f; // vertical window spacing <cm>
  //float th = 0.004; // PDMS layer thickness <cm>
  
  // Simulation Time (user input)
  float sim_time = 10.0f; // simulation time <s> was 360.0f
  if (argc == 2)
    sim_time = atof(argv[1]);
  float print_frequency = 0.25f; // print frequency <s> was 0.25f
  
  // Write-Out Schedule
  // 0-10s: 1s, 10-30s: 5s, 30-180s: 30s
  print_scheduler print_time(print_frequency);
  //print_time.schedule(10.0f,5.0f); // (start_time <s>,frequency <s>) 
  //print_time.schedule(30.0f,30.0f);  
  
  // Initialize Computational Domain (user input)
  int Nx = 288; // 576, 288, 144
  int Ny = 288; // 576, 288, 144
  int Nz = 96; // 192, 96, 48 
  float dt = 1e-8; // was 1e-6
  
  // Output Filename (user input)
  string dir = "out/PO2/";
  string filename = "step";
  // Plane Slices
  int dimY = 1;
  int dimZ = 2;
  int slc = 71;
  int slc0 = 5; // 25 um
  int slc1 = 9; // 50 um
  int slc2 = 13; // 75 um
  int slc3 = 17; // 100 um
  
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
  geometry geo(L,H,W,l,h,xs,ys);
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
  cout << "[dx,dy,dz,dt] = " << "[" << dx << "," << dy << "," << dz << "," << dt << "]\n\n";
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
  varIC(u_h,"out/PO2/steady-state.csv",N);
  //constIC(u_h,1.0f,N);/*
  print(u_h,Nx,Ny,Nz,dimY,slc,dir+filename+"_yPlane0.csv");
  print(u_h,Nx,Ny,Nz,dimZ,slc0,dir+filename+"_zPlane025um0.csv");
  print(u_h,Nx,Ny,Nz,dimZ,slc1,dir+filename+"_zPlane050um0.csv");
  print(u_h,Nx,Ny,Nz,dimZ,slc2,dir+filename+"_zPlane075um0.csv");
  print(u_h,Nx,Ny,Nz,dimZ,slc3,dir+filename+"_zPlane100um0.csv");//*/
  
  
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
  float t = 0.0f; int np = 1; time_writer write_time(dir+"t_"+filename+".csv"); write_time(t*tau);
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
      cout << u_h[0] << endl;
      //print(u_h,N,dir+filename+".csv");/*
      std::string id = (print_frequency==sim_time)? "":to_string(np);
      print(u_h,Nx,Ny,Nz,dimY,slc,dir+filename+"_yPlane"+id+".csv");
      print(u_h,Nx,Ny,Nz,dimZ,slc0,dir+filename+"_zPlane025um"+id+".csv");
      print(u_h,Nx,Ny,Nz,dimZ,slc1,dir+filename+"_zPlane050um"+id+".csv");
      print(u_h,Nx,Ny,Nz,dimZ,slc2,dir+filename+"_zPlane075um"+id+".csv");
      print(u_h,Nx,Ny,Nz,dimZ,slc3,dir+filename+"_zPlane100um"+id+".csv");//*/
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
