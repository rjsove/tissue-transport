#include <iostream>
#include <string>
#include <fstream>
using namespace std;

int main(int argc,char** argv)
{
  string filename = argv[1];
  double min = 1e25, temp;
  ifstream file;
  file.open(filename);
  while (!file.eof())
  {
    file >> temp;
    if (temp < min)
      min = temp;
  }
  
  cout << "\nmin = " << min << endl << endl; 
  
  return 0;
}