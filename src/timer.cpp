#include "timer.h"
#include <iostream>

timer::timer() : name("\b")
{
  start = std::chrono::high_resolution_clock::now();
}

timer::timer(const std::string& name) : name(name)
{
  start = std::chrono::high_resolution_clock::now();
}

timer::~timer()
{
  // Calculate runtime
  auto stop = std::chrono::high_resolution_clock::now();
  auto runtime = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
  
  // Display Runtime
  std::cout << name << " Runtime: ";
  display_runtime(runtime.count());
}

void display_runtime(long runtime)
{
  long remainder = runtime;
  
  // Calculate Milliseconds
  int ms = remainder%1000;
  remainder = remainder/1000;
  
  // Calculate Seconds
  int s = remainder%60;
  remainder = remainder/60;
  
  // Calculate Minutes
  int mins = remainder%60;
  remainder = remainder/60;
  
  // Calculate Hours and Days
  int hrs = remainder%24;
  int days = remainder/24;
  
  // Output Runtime 
  if (days != 0)
  {
    if (days == 1)
      std::cout << "1 day ";
    else
      std::cout << days << " days "; 
  }
  if (hrs != 0)
  {
    if (hrs == 1)
      std::cout << "1 hour ";
    else
      std::cout << hrs << " hours ";
  }
  if (mins != 0)
  {
    if (mins == 1)
      std::cout << "1 minute ";
    else
      std::cout << mins << " minutes ";
  }
  if (s != 0)
  {
    if (s == 1)
      std::cout << "1 second ";
    else
      std::cout << s << " seconds ";
  }
  if (days == 0 && hrs == 0 && mins == 0 && s == 0)
  {
    if (ms == 1)
      std::cout << "1 millisecond";
    else
      std::cout << ms << " milliseconds";
  }
  std::cout <<  "\n";
}