#include <iostream>
#include <string>
#include <fstream>
using namespace std;

int main(int argc,char** argv)
{
  string filename = argv[1];
  double sum = 0.0,temp;
  int N = 0;
  ifstream file;
  file.open(filename);
  while (!file.eof())
  {
    N++;
    file >> temp;
    sum+=temp;
  }
  
  cout << "\nmean = " << sum/N << endl << endl; 
  
  return 0;
}