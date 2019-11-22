#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
using namespace std;

int main()
{
  // Probe GPU
  cout << endl << endl;
  cout << "GPU Specifications:\n";
  cudaDeviceProp dev;
  cudaGetDeviceProperties(&dev,0);
  cout << "Total Device Memory: "
    << dev.totalGlobalMem/1073741824 << " GB\n"; 
  cout << "Shared Memory Per Block: "
    << dev.sharedMemPerBlock/1024 << " kB\n";
  cout << "Maximum Threads Per Block: " 
    << dev.maxThreadsPerBlock << "\n";
  cout << "Maximum Block Dimensions: ("
    << dev.maxThreadsDim[0] << ","
    << dev.maxThreadsDim[1] << ","
    << dev.maxThreadsDim[2] << ")\n";
  cout << "Maximum Grid Size: (" 
    << dev.maxGridSize[0] << ","
    << dev.maxGridSize[1] << "," 
    << dev.maxGridSize[2] << ")\n";
  cout << endl << endl;

  return 0;
}
