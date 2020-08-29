// Solve O2 Transport on Tissue Volume Using GPGPU Acceleration
#include <iostream>
#include <string>
#include <cuda.h>
#include "IOTools.h"
#include "PDETools.h"
#include "PDEsolve.h"
#include "configuration.h"
#include "timer.h"

// User Inputs
/*
Usage:
    $ runGPU <path to project directory>

Project Directory Contains Files 
    - physics.ttgpu
    - geometry.ttgpu
    - boundary_conditions.ttgpu
    - discretization.ttgpu
    - output_settings.ttgpu
*/

/*
// Physical Constants (user input)
float D = 2.41e-5f; // diffusivity <cm^2/s>
float k = 3.89e-5f; // solubility <mLO2/mL/mmHg>
float Dpdms = 3.40e-5f; // diffusivity <cm^2/s> [Merkel 2000]
float kpdms = 2.37e-4f; // solubility <mLO2/mL/mmHg> [Merkel 2000]
float M0 = 1.50e-4f; // O2 consumption <mLO2/mL/s>
float Pcrit = 0.5f; // critical PO2 <mmHg>
float Pc = 48.0f; // capillary PO2 <mmHg>
float K = 30.0f; // capillary source <mmHg/s>

// Geometry (user input)
float L = 0.52f; // tissue length <cm> regular simulation: [0.2x0.2]
float H = 0.44f; // tissue height <cm> five window simulation [0.52x0.44]
float W = 0.06f; // tissue depth <cm>
float l = 0.04f; // window length <cm>
float h = 0.02f; // window height <cm>
float xs = 0.1f; // horizontal window spacing <cm>
float ys = 0.2f; // vertical window spacing <cm>
float th = 0.002; // PDMS layer thickness <cm>

// Boundary Condition Parameters (user input)
float Pbsl = 53.2.0f; // Baseline PO2 <mmHg> 7%
float Phigh = 91.2f; // High PO2 <mmHg> 12% 
float Plow = 15.2f; // Low PO2 <mmHg> 2%

// Computational Domain (user input)
int Nx = 288; // 576, 288, 144
int Ny = 288; // 576, 288, 144
int Nz = 96; // 192, 96, 48 
float dt = 1e-8; // was 1e-6

// Simulation Time (user input)
float sim_time = 750.0f; // simulation time <s> was 360.0f
if (argc == 2)
  sim_time = atof(argv[1]);
float print_frequency = sim_time; // print frequency <s> was 0.25f

// Output Filename (user input)
string dir = "out/PO2/";
string filename = "five-windows";
// z-Plane Slice
//int dimY = 1;
int dimZ = 2;
//int slc = round(Ny-1)/2
int slc0 = 5; // 25 um
int slc1 = 9; // 50 um
int slc2 = 13; // 75 um
int slc3 = 17; // 100 um
*/

int main(int argc,char** argv)
{
  // Start Main Timer 
  timer timer1("Total");
  
  // Read User Inputs
  ttgpu::inStream consts = ttgpu::read(ttgpu::PHYSICS)
  ttgpu::inStream geom   = ttgpu::read(ttgpu::GEOMETRY)
  ttgpu::inStream ivbp   = ttgpu::read(ttgpu::BOUNDARY_CONDITIONS)
  ttgpu::inStream discr  = ttgpu::read(ttgpu::DISCRETIZATION)
  ttgpu::inStream outset = ttgpu::read(ttgpu::OUTPUT_SETTINGS)
   
  // Physical Constants 
  float D     = consts.D;     // diffusivity <cm^2/s>
  float k     = consts.k;     // solubility <mLO2/mL/mmHg>
  float Dpdms = consts.Dpdms; // diffusivity <cm^2/s> [Merkel 2000]
  float kpdms = consts.kpdms; // solubility <mLO2/mL/mmHg> [Merkel 2000]
  float M0    = consts.M0;    // O2 consumption <mLO2/mL/s>
  float Pcrit = consts.Pcrit; // critical PO2 <mmHg>
  float Pc    = consts.Pc;    // capillary PO2 <mmHg>
  float K     = consts.K;     // capillary source <mmHg/s>
  
  // Geometry 
  // regular simulation:     [LxW]=[0.20x0.20]
  // five window simulation: [LxW]=[0.52x0.44]
  float L  = geom.L;  // tissue length <cm> 
  float H  = geom.H;  // tissue height <cm> 
  float W  = geom.W;  // tissue depth <cm>
  float l  = geom.l;  // window length <cm>
  float h  = geom.h;  // window height <cm>
  float xs = geom.xs; // horizontal window spacing <cm>
  float ys = geom.ys; // vertical window spacing <cm>
  float th = geom.th; // PDMS layer thickness <cm>
  
  // Boundary Condition Parameters 
  float Pbsl = ibvp.bcValue; // Baseline Window PO2 <mmHg>
  if (!ibvp.constantBC)
  {
    float Plow  = ibvp.bcSqWvLowValue;  // Low PO2 <mmHg> 
    float Phigh = ibvp.bcSqWvHighValue; // High PO2 <mmHg> 
  }
  
  // Computational Domain
  int   Nx = discr.Nx;
  int   Ny = discr.Ny;
  int   Nz = discr.Nz; 
  float dt = discr.dt;
  
  // Simulation Time 
  float sim_time        = ouset.sim_time;        // simulation time <s> 
  float write_interval  = outset.write_interval; // print interval <s> 
  
  // Output Filename 
  string outputFilename = outset.filename;
  
  // Write-Out Schedule
  print_scheduler print_time(write_interval);
  
  // Calculate Dimensionless Parameters
  float tau = L*L/D;
  float alpha = K/Pc*tau;
  float ub = 1;
  float beta = M0/(Pc*k)*tau; 
  float km = Pcrit/Pc;
  float lambda = Dpdms/D;
  float sigma = lambda*kpdms/k;
  float ay = H/L;
  float az = W/L;
  
  // Calculate Computational Parameters
  model mdl(alpha,beta,1.0f,ub,km,lambda,sigma);
  grid grd(Nx,Ny,Nz,dx,dy,dz,dt);
  geometry geo(L,H,W,l,h,xs,ys);
  int numNodes = Nx*Ny*Nz;
  size_t size = numNodes*sizeof(float);
  float dx = 1.0f/(Nx-1.0f);
  float dy = ay/(Ny-1.0f);
  float dz = az/(Nz-1.0f);
  float endTime = sim_time/tau;
  
  // Print Simulation Info
  float ds = min(min(dx,dy),dz); 
  std::cout << "\n\n---Oxygen---\n\n";
  std::cout << "\n\nSimulation Parameters\n\n";
  std::cout << "Runtime: " << sim_time << " s \t\t[" << endTime/dt << " time steps]\n\n";
  std::cout << "dt/dx^2 = " << dt/ds/ds << "\n\n";
  std::cout << "[dx,dy,dz,dt] = " << "[" << dx << "," << dy << "," << dz << "," << dt << "]\n\n";
  std::cout << "Length Scale: " << L << " cm\n";
  std::cout << "Time Scale: " << tau << " s\n";
  std::cout << "O2 Scale: " << Pc << " mmHg\n\n";
  if (ibvp.constantBC)
    std::cout << "Window O2: " << Pbsl/760 << "%\n\n";
  else
  {
    std::cout << "Window O2: Square-Wave\n";
    std::cout << "Low O2: "   << Plow/760  << "%\n";
    std::cout << "Baseline: " << Pbsl/760  << "%\n";
    std::cout << "High O2: "  << Phigh/760 << "%\n\n";
  }
  std::cout << "PDMS Layer Thickness: " << round(geom.th/dz)*dz*1e4 << " um\n\n";
  
  // Allocate Memory on Host
  float* u_h = new float[numNodes];
  if (ibvp.homogeneousIC)
    constIC(u_h,ibvp.icValue,numNodes);
  else
    varIC(u_h,ibvp.icFilename,numNodes);
  
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
  float t = 0.0f; int np = 1; 
  time_writer write_time(dir+"t.csv"); write_time(t*tau);
  float uwin;
  for (int nt = 1; t < endTime; nt++)
  { 
    // Boundary Condition
    if (ibvp.constantBC)
      uwin = Pbsl/Pc;
    else
      uwin = squareWave(t,endTime,Pbsl/Pc,Phigh/Pc,Plow/Pc); 
   
    // Call GPU Kernel
    step<<<dimGrid,dimBlock>>>(uold_d,unew_d,uwin,mdl,grd,geo);
    t += dt;
    
    // Print Solution
    if (print_time(t*tau))
    {
      // Time Write Out
      timer timer2("Done.");
      std::cout << "Writing t = " << t << "...\n";
      
      // Copy Memory from GPU
      cudaMemcpy(u_h,unew_d,size,cudaMemcpyDeviceToHost);
      
      // Call Prints
      if (outset.writeFullSolution)
        print(u_h,size,filname+"csv");
      for (int slice = 0; slice < outset.numberOfSlices; slice++)
      {
        print(u_h,Nx,Ny,Nz,ouset.sliceDim[slice],outset.slice[slice],filename+outset.fileSuffix[slice]);
      }
      write_time(t*tau);
      
      // Update Print Counter
      np++;
    }
  }
  
  // Delete Pointers
  delete[] u_h;
  cudaFree(uold_d);
  cudaFree(unew_d);
  
  return 0;
}