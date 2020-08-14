// PDE Solver for the PDE:
//
//         du/dt = lapl(u) alpha*(ub-u) - beta*u/(km+u) 
//         u(x,0) = u0(x)
//         u(x_d,t) = ub(t)
//         du(x_n,t)/dn = 0
//

#ifndef PDESOLVE_H
#define PDESOLVE_H

#include "PDETools.h"

__device__ int Nx,Ny,Nz;
__device__ float dt,dx,dy,dz;
__device__ float L,H,W,l,h,xs,ys;
__device__ float alpha,beta,gam,ub,km,lambda,sigma; 
__device__ float* swap;
__device__ int at(int i,int j,int k);
__device__ float CDM(float* u,int i,int j,int k);
__device__ float CDMi(float* u,int i,int j,int k);
__device__ void imex(float* u_old,float* u_new,float BC,int i,int j,int k);
__global__ void step(float* uold_d,float* unew_d,float BC,model,grid,geometry);
__device__ bool windowBC(int i,int j,int k);
__device__ bool fiveWindowsBC(int i,int j,int k);

#endif