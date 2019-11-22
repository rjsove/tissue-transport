#include "PDEsolve.h"

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
  return (u[at(i+1,j,k)]-2*u[at(i,j,k)]+u[at(i-1,j,k)])/dx/dx 
       + (u[at(i,j+1,k)]-2*u[at(i,j,k)]+u[at(i,j-1,k)])/dy/dy
       + (u[at(i,j,k+1)]-2*u[at(i,j,k)]+u[at(i,j,k-1)])/dz/dz;
}

// GPU Forward Difference
__device__ void forwardEuler(float* u_old,float* u_new,float alpha,float ub,float beta,float km,float BC,float h,int i,int j,int k)
{
  // Apply Boundary Conditions
  int n = at(i,j,k);
  // Fixed Value BC at Window 
  if (windowBC(i,j,k))
  {
    u_new[n] = BC;
  }
  // Zero Flux  
  else
  {
    u_new[n] = u_old[n] + h*(CDM(u_old,i,j,k)+alpha*(ub-u_old[at(i,j,k)])-beta*u_old[at(i,j,k)]/(km+u_old[at(i,j,k)]));
  }
}

// GPU Implicit-Explicit
__device__ void imex(float* u_old,float* u_new,float alpha,float ub,float beta,float km,float BC,float h,int i,int j,int k)
{
  // Apply Boundary Conditions
  int n = at(i,j,k);
  // Fixed Value BC at Window 
  if (windowBC(i,j,k))
  {
    u_new[n] = BC;
  }
  // Zero Flux  
  else
  {
    u_new[n] = (u_old[n] + h*(CDM(u_old,i,j,k) + alpha*ub - beta*u_old[at(i,j,k)]/(km+u_old[at(i,j,k)])))/(1+alpha*h);
  }
}

// GPU Time Step Kernel
__global__ void step(float* u_old,float* u_new,float alpha,float ub,float beta,float km,float BC,
                     int nx,int ny,int nz,float hx,float hy,float hz,float h)
{ 
  // Set GPU Variables
  Nx = nx; dx = hx;
  Ny = ny; dy = hy;
  Nz = nz; dz = hz;
  
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
        imex(u_old,u_new,alpha,ub,beta,km,BC,h,i,j,k);
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
  // Domain Dimensions (cm)
  // *must be consistant with simulation
  float W = 0.2f,H = 0.2f; 

  // Window Dimensions (cm)
  // float w = 0.04f,h = 0.02f; // my window
  float w = 0.06f,h = 0.03f; // Graham's window
  
  // Relative Window Dimensions 
  float w_ = w/W,h_ = h/H;
   
  return abs(2*i*dx-1)<=w_ & abs(2*j*dy-1)<= h_ & k==0;
}
// Five Windows Boundary Condition
__device__ bool fiveWindowsBC(int i,int j,int k)
{
  // Domain Dimensions (cm)
  // *must be consistant with simulation
  float W = 0.52f,H = 0.44f; 

  // Window Dimensions (cm)
  float w = 0.04f,h = 0.02f; 
  
  // Window Spacing (cm)
  float xs = 0.1f,ys = 0.2f;
  
  // Relative Window Dimensions 
  float w_ = w/W,h_ = h/H;
  float xs_ = xs/W,ys_ = ys/H;
  
  // Windows
  float x0,y0;
  x0 = (1+w_+xs_)/2; y0 = (1+h_+ys_)/2;
  bool window1 = abs(2*(i*dx-x0))<=w_ & abs(2*(j*dy-y0))<=h_ & k==0; 
  x0 = (1-w_-xs_)/2; y0 = (1+h_+ys_)/2;
  bool window2 = abs(2*(i*dx-x0))<=w_ & abs(2*(j*dy-y0))<=h_ & k==0; 
  x0 = 0.5-w_-xs_; y0 = (1-h_-ys_)/2;
  bool window3 = abs(2*(i*dx-x0))<=w_ & abs(2*(j*dy-y0))<=h_ & k==0; 
  x0 = 0.5; y0 = (1-h_-ys_)/2;
  bool window4 = abs(2*(i*dx-x0))<=w_ & abs(2*(j*dy-y0))<=h_ & k==0; 
  x0 = 0.5+w_+xs_; y0 = (1-h_-ys_)/2;
  bool window5 = abs(2*(i*dx-x0))<=w_ & abs(2*(j*dy-y0))<=h_ & k==0; 
  
  return window1 | window2 | window3 | window4 | window5;
}
