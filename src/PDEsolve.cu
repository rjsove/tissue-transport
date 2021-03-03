#include "PDEsolve.h"

#define WINDOWBC windowBC(i,j,k)
//#define WINDOWBC fiveWindowsBC(i,j,k) // need to change domain dimensions and geometry for window spacing

// Comment out if not using PDMS layer
#define PDMS
#define INTERFACE 5 // make sure this corresponds to number of nodes in PDMS layer

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
  return (1+sigma)*((u[at(i+1,j,k)]-2*u[at(i,j,k)]+u[at(i-1,j,k)])/dx/dx/2    
       + (u[at(i,j+1,k)]-2*u[at(i,j,k)]+u[at(i,j-1,k)])/dy/dy/2)
       + (u[at(i,j,k+1)]-u[at(i,j,k)])/dz/dz
       + sigma*(u[at(i,j,k-1)]-u[at(i,j,k)])/dz/dz; 
}

// GPU Implicit-Explicit
__device__ void imex(float* u_old,float* u_new,float BC,int i,int j,int k)
{
  // Apply Boundary Conditions
  int n = at(i,j,k);
  // Fixed Value BC at Window 
  if (WINDOWBC)
  {
    u_new[n] = BC;
  }
  #ifdef PDMS 
  // Interface B.C.
  else if (k==INTERFACE)
  {
    u_new[n] = (u_old[n] + dt*lambda/(lambda+sigma)*(2*CDMi(u_old,i,j,k) 
             + alpha*ub - beta*u_old[at(i,j,k)]/(km+u_old[at(i,j,k)])))/(1+alpha*dt*lambda/(lambda+sigma));    
  } 
  // PDMS and Zero Flux Boundaries
  else if (k<INTERFACE&&k>0)
  {
    u_new[n] = (u_old[n] + dt*(lambda*CDM(u_old,i,j,k)))/(1+alpha*dt);
  } 
  #endif
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
  dt  = grd.dt;
  Nx = grd.Nx; dx = grd.dx;
  Ny = grd.Ny; dy = grd.dy;
  Nz = grd.Nz; dz = grd.dz;
  L = geo.L; l = geo.l;
  H = geo.H; h = geo.h;
  xs = geo.xs; ys = geo.ys;
  
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
    printf("(%d,%d,%d,%d)\n",blockIdx.x,blockIdx.y,threadIdx.x,threadIdx.y);
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
  float l_ = l/L,h_ = h/L;
  float a = H/L;
   
  return abs(2*i*dx-1)<=l_ & abs(2*j*dy-a)<= h_ & k==0;
}
// Five Windows Boundary Condition
__device__ bool fiveWindowsBC(int i,int j,int k)
{
  // Domain Dimensions (cm)
  // float W = 0.52f,H = 0.44f; 
  // Window Dimensions (cm)
  // float w = 0.04f,h = 0.02f; 
  // Window Spacing (cm)
  // float xs = 0.1f,ys = 0.2f;
  
  // Relative Window Dimensions 
  float l_ = l/L,h_ = h/L; 
  float xs_ = xs/L,ys_ = ys/L; 
  float a = H/L; 
  
  // Windows
  float x0[] = {(1+l_+xs_)/2,(1-l_-xs_)/2,0.5f-l_-xs_,0.5f,0.5f+l_+xs_}; 
  float y0[] = {(a+h_+ys_)/2,(a+h_+ys_)/2,(a-h_-ys_)/2,(a-h_-ys_)/2,(a-h_-ys_)/2}; 
  
  bool window = false;
  for (int win = 0; win < 5; win++)
  {
    window = window | abs(2*(i*dx-x0[win]))<=l_ & abs(2*(j*dy-y0[win]))<=h_ & k==0; 
  }
  
  return window;
}
