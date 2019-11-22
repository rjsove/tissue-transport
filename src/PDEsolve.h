// PDE Solver for the PDE:
//
//         du/dt = lapl(u) alpha*(ub-u) - beta*u/(km+u) 
//         u(x,0) = u0(x)
//         u(x_d,t) = ub(t)
//         du(x_n,t)/dn = 0
//

#ifndef PDESOLVE_H
#define PDESOLVE_H

__device__ int Nx,Ny,Nz;
__device__ float dx,dy,dz; 
__device__ float* swap;
__device__ int at(int i,int j,int k);
__device__ float CDM(float* u,int i,int j,int k);
__device__ void forwardEuler(float* u_old,float* u_new,float alpha,float ub,float beta,float km,float BC,float h,int i,int j,int k);
__device__ void imex(float* u_old,float* u_new,float alpha,float ub,float beta,float km,float BC,float h,int i,int j,int k);
__global__ void step(float* uold_d,float* unew_d,float alpha,float ub,float beta,float km,float BC,
                    int nx,int ny,int nz,float hx,float hy,float hz,float dt);
__device__ bool windowBC(int i,int j,int k);
__device__ bool fiveWindowsBC(int i,int j,int k);

#endif