#include "IOTools.h"
#include <iostream>
#include <string>
#include <fstream>

void ttgpu::physics::readFile()
{
  // Initialize File Stream
  ifstream file("physics.ttpgu");
  
  // Read in Parameters
  temp        << file;
  this->D     << file;
  temp        << file;
  this->k     << file;
  temp        << file;
  this->Dpdms << file;
  temp        << file;
  this->kpdms << file;
  temp        << file;
  this->M0    << file;
  temp        << file;
  this->Km    << file;
  temp        << file;
  this->Pc    << file;
  temp        << file;
  this->K     << file;
  file.close();
}

void ttgpu::geometry::readFile()
{
  std::String temp;
  
  // Initialize File Stream
  ifstream file("geometry.ttpgu");
  
  // Read in Parameters
  temp        << file;
  this->L     << file;
  temp        << file;
  this->H     << file;
  temp        << file;
  this->W     << file;
  temp        << file;
  this->l     << file;
  temp        << file;
  this->h     << file;
  temp        << file;
  this->xs    << file;
  temp        << file;
  this->ys    << file;
  temp        << file;
  this->th    << file;
  file.close();
}

ttgpu::boundaryConditions::readFile()
{
  std::string temp,icTemp;
  
  // Initialize File Stream
  ifstream file("physics.ttpgu");
  temp          << file;
  icTemp        << file;
  temp          << file;
  this->bcValue << file;
  
  // Check if bcValue is at eof 
  if ()
  {
    this->bcSqWvLowValue  << file;
    this->bcSqWvHighValue << file;  
  }
  
  // Close File
  file.close();
  
  // Check if icTemp is Filename or Value 
  if ()
  {
    this->icValue = icTemp;
  }
  else
  {
    this->icFilename = icTemp;
  }
}

ttgpu::discretization::readFile()
{
  ifstream file("physics.ttpgu");
  temp        << file;
  this->D     << file;
  temp        << file;
  this->k     << file;
  temp        << file;
  this->Dpdms << file;
  temp        << file;
  this->kpdms << file;
  temp        << file;
  this->M0    << file;
  temp        << file;
  this->Km    << file;
  temp        << file;
  this->Pc    << file;
  temp        << file;
  this->K     << file;
  file.close();
}

ttgpu::outputSettings::readFile()
{
  // Initialize File Stream
  ifstream file("output_settings.ttpgu");
  
  // Read in Parameters
  std::string temp;
  temp                    << file;
  this->filename          << file;
  temp                    << file;
  this->writeFullSolution << file;
  temp                    << file;
  this->numberOfSlices    << file;
  
  // Initialize Arrays
  this->sliceDim = new int[numberOfSlices];
  this->slice = new int[numberOfSlices];
  this->fileSuffix = new std::string[numberOfSlices];
  
  // Loop Through the Rest of the File
  for (int i = 0; i < 3; i++)
  {
    temp << file;
    for (int j = 0; j < this->numberOfSlices; j++)
    {
      switch (i)
      {
        case 0:
          this->slideDim[j] << file;
          break;
        case 1:
          this->slice[j] << file;
          break;
        case 2:
          this->fileSuffix[j] << file;
          break;
      }
    }
    file.close();
  }
  
  this->Dpdms << file;
  temp        << file;
  this->kpdms << file;
  temp        << file;
  this->M0    << file;
  temp        << file;
  this->Km    << file;
  temp        << file;
  this->Pc    << file;
  temp        << file;
  this->K     << file;
  file.close();
}

ttgpu::inStream ttgpu::read(streamType type)
{
  std::string filename;
  ifstream file;
  
  // Determine Which File to Read
  switch (type)
  {
    case ttgpu::PHYSICS:
      filename = "physics.ttgpu";
      break;
    case ttgpu::GEOMETRY:
      filename = "geometry.ttgpu";
      break;
    case ttgpu::BOUNDARY_CONDITIONS:
      filename = "boundary_conditions.ttgpu";
      break;
    case ttgpu::DISCRETIZATION:
      filename = "discretization.ttgpu";
      break;
    case ttgpu::OUTPUT_SETTINGS:
      filename = "output_settings.ttgpu";
      break;
    default:
      std::cerr << "ttgpu: incorrect usage of read\n";
      break; 
  }
}