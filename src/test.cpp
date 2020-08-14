#include <iostream>
#include "PDETools.h"

int main()
{
  // Create Test Data on the Heap
  size_t size_x = 5;
  size_t size_y = 4;
  size_t size_z = 3;
  size_t size = size_x*size_y*size_z;  
  float* u = new float[size];
  
  // Assign Values
  int k = 0; float val = 1.0;
  for (int j = 0; j < size_y; j++)
  {
    for (int i = 0; i < size_x; i++)
    {
      u[i+size_x*(j+size_y*k)] = val;
      val++;
    }
  }
  k = 1; val = 20;
  for (int j = 0; j < size_y; j++)
  {
    for (int i = 0; i < size_x; i++)
    {
      u[i+size_x*(j+size_y*k)] = val;
      val+=20;
    }
  }
  k = 2; val = 11;
  for (int j = 0; j < size_y; j++)
  {
    for (int i = 0; i < size_x; i++)
    {
      u[i+size_x*(j+size_y*k)] = val;
      val+=11;
    }
  }
  
  
  for (int k = 0; k < size_z; k++)
  {
    for (int j = 0; j < size_y; j++)
    {
      for (int i = 0; i < size_x; i++)
      {
        std::cout << u[i+size_x*(j+size_y*k)] << "  ";
  
      }
      std::cout << "\n";
    }
    std::cout << "\n\n";
  }
  
  return 0;
}