#include "PDEsolve.h"

#define WINDOWBC windowBC(i,j,k)
//#define WINDOWBC fiveWindowsBC(i,j,k)

// Comment out if not using PDMS layer
#define PDMS
#define INTERFACE 6 // make sure this corresponds to number of node in PDMS layer

// Solves PDE du/ds = lapl(u) + alpha*(ub-u) - beta*u/(km+u)

// GPU Index Evaluation
__device__ int at(int i,int j,int k)
{
  if (i == -1) i+=2;
  else if (i == Nx) i-=2;
  if (j == -1) j+=2;
  else if (j == Ny) j-=2;
  if (k == -1) k+=2;
  else if (k == Nz) k-=2; 
  return i+Nx*(j+Ny*k);
}

// GPU Central Difference
__device__ float CDM(float* u,int i,int j,int k)
{ 
  return gam*(u[at(i+1,j,k)]-2*u[at(i,j,k)]+u[at(i-1,j,k)])/dx/dx 
       + (u[at(i,j+1,k)]-2*u[at(i,j,k)]+u[at(i,j-1,k)])/dy/dy
       + (u[at(i,j,k+1)]-2*u[at(i,j,k)]+u[at(i,j,k-1)])/dz/dz;
}

// Central Difference at Interface (isotropic only: gamma==1.0)
__device__ float CDMi(float* u,int i,int j,int k)
{
  return (u[at(i+1,j,k)]-2*u[at(i,j,k)]+u[at(i-1,j,k)])/dx/dx/2   
       + (u[at(i+1,j,k)]-2*u[at(i,j,k)]+u[at(i-1,j,k)])/dx/dx*sigma/2 
       + (u[at(i,j+1,k)]-2*u[at(i,j,k)]+u[at(i,j-1,k)])/dy/dy/2
       + (u[at(i,j+1,k)]-2*u[at(i,j,k)]+u[at(i,j-1,k)])/dy/dy*sigma/2
       + (u[at(i,j,k+1)]-u[at(i,j,k)])/dz/dz*lambda;
       + (u[at(i,j,k-1)]-u[at(i,j,k)])/dz/dz; 
}

// GPU Implicit-Explicit
__device__ void imex(float* u_old,float* u_new,float BC,int i,int j,int k)
{
  // Apply Boundary Conditions
  int n = at(i,j,k);
  // Fixed Value BC at Window 
  if (k==0)
  {
    u_new[n] = 1.0f;
  }
  //if (WINDOWBC)
  //{
  //  u_new[n] = BC;
  //}
  //#ifdef PDMS 
  // Interface B.C.
  else if (k==INTERFACE)
  {
    u_new[n] = 2.0f;//(u_old[n] + dt*(2*(CDMi(u_old,i,j,k)) + alpha*ub - beta*u_old[at(i,j,k)]/(km+u_old[at(i,j,k)])))/(2+alpha*dt);    
  } 
  // PDMS and Zero Flux Boundaries
  else if (k<INTERFACE&&k>0)
  {
    u_new[n] = 3.0f;//(u_old[n] + dt*(lambda*CDM(u_old,i,j,k)))/(1+alpha*dt);
  } 
  //#endif
  else // Tissue and Zero Flux Boundaries
  {
    u_new[n] = (u_old[n] + dt*(CDM(u_old,i,j,k) + alpha*ub - beta*u_old[at(i,j,k)]/(km+u_old[at(i,j,k)])))/(1+alpha*dt);
  }
}

// GPU Time Step Kernel
__global__ void step(float* u_old,float* u_new,float BC,model mdl,grid grd,geometry geo)
{ 
  // Set GPU Variables
  alpha = mdl.alpha; ub = mdl.ub;
  beta = mdl.gamma; km = mdl.km;
  gam = mdl.gamma; lambda = mdl.lambda;
  sigma = mdl.sigma;
  Nx = grd.Nx; dx = grd.dx;
  Ny = grd.Ny; dy = grd.dy;
  Nz = grd.Nz; dz = grd.dz;
  L = geo.L; l = geo.l;
  H = geo.H; h = geo.h;
  
  // Determine Position within array
  int ATOM_SIZE_X = Nx/(blockDim.x*gridDim.x);
  int ATOM_SIZE_Y = Ny/(blockDim.y*gridDim.y);
  int ATOM_SIZE_Z = Nz/(blockDim.z*gridDim.z);
  int i0 = (threadIdx.x+blockIdx.x*blockDim.x)*ATOM_SIZE_X;
  int j0 = (threadIdx.y+blockIdx.y*blockDim.y)*ATOM_SIZE_Y;
  int k0 = (threadIdx.z+blockIdx.z*blockDim.z)*ATOM_SIZE_Z;

  int i,j,k,n; 
   
  for (i = i0; i < i0+ATOM_SIZE_X; i++) 
  {
    for (j = j0; j < j0+ATOM_SIZE_Y; j++)
    {
      for (k = k0; k < k0+ATOM_SIZE_Z; k++)
      {
        // Call Time Scheme
        imex(u_old,u_new,BC,i,j,k);
      }
    }
  }
  
  __syncthreads(); 
  
  for (i = i0; i < i0+ATOM_SIZE_X; i++) 
  {
    for (j = j0; j < j0+ATOM_SIZE_Y; j++)
    {
      for (k = k0; k < k0+ATOM_SIZE_Z; k++)
      {
        n = i+Nx*(j+Ny*k);
        u_old[n] = u_new[n];
      }
    }
  }
  
  __syncthreads(); 
}

// Boundary Conditions
// Window Boundary Condition
__device__ bool windowBC(int i,int j,int k)
{   
  // Relative Window Dimensions 
  float l_ = l/L,h_ = h/H;
   
  return abs(2*i*dx-1)<=l_ & abs(2*j*dy-1)<= h_ & k==0;
}
// Five Windows Boundary Condition
__device__ bool fiveWindowsBC(int i,int j,int k)
{
  // Domain Dimensions (cm)
  // *must be consistant with simulation
  // float W = 0.52f,H = 0.44f; 

  // Window Dimensions (cm)
  // float w = 0.04f,h = 0.02f; 
  
  // Window Spacing (cm)
  // float xs = 0.1f,ys = 0.2f;
  
  // Relative Window Dimensions 
  float l_ = l/L,h_ = h/H;
  float xs_ = xs/L,ys_ = ys/H;
  
  // Windows
  float x0,y0;
  x0 = (1+l_+xs_)/2; y0 = (1+h_+ys_)/2;
  bool window1 = abs(2*(i*dx-x0))<=l_ & abs(2*(j*dy-y0))<=h_ & k==0; 
  x0 = (1-l_-xs_)/2; y0 = (1+h_+ys_)/2;
  bool window2 = abs(2*(i*dx-x0))<=l_ & abs(2*(j*dy-y0))<=h_ & k==0; 
  x0 = 0.5-l_-xs_; y0 = (1-h_-ys_)/2;
  bool window3 = abs(2*(i*dx-x0))<=l_ & abs(2*(j*dy-y0))<=h_ & k==0; 
  x0 = 0.5; y0 = (1-h_-ys_)/2;
  bool window4 = abs(2*(i*dx-x0))<=l_ & abs(2*(j*dy-y0))<=h_ & k==0; 
  x0 = 0.5+l_+xs_; y0 = (1-h_-ys_)/2;
  bool window5 = abs(2*(i*dx-x0))<=l_ & abs(2*(j*dy-y0))<=h_ & k==0; 
  
  return window1 | window2 | window3 | window4 | window5;
}