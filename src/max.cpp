#include <iostream>
#include <string>
#include <fstream>
using namespace std;

int main(int argc,char** argv)
{
  string filename = argv[1];
  double max = 0.0, temp;
  ifstream file;
  file.open(filename);
  while (!file.eof())
  {
    file >> temp;
    if (temp > max)
      max = temp;
  }
  
  cout << "\nmax = " << max << endl << endl; 
  
  return 0;
}