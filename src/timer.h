#pragma once 

#include <chrono>
#include <string>

class timer 
{
private:
  // Private Member
  std::string name;
  std::chrono::time_point<std::chrono::high_resolution_clock> start;
  
public:
  // Constructor
  timer();
  timer(const std::string& timer_name);
  
  // Destructor
  ~timer();
}; 

void display_runtime(long runtime);