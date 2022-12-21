// This is ROOT macro.
// root -l -q test.cpp+O
#include <stdlib.h>
#include "TSorter.cpp"
using namespace std;

int main(int argc, char *argv[])
{
  TSorter sorter;
  sorter.LoadFiles(atoi(argv[1]), 0, 250);
  // sorter.LoadFiles(2267, 0, 47);
  //  sorter.LoadFiles(2223, 0, 49);  
  //  sorter.LoadFiles(2139, 0);
  sorter.SortData();
  sorter.WriteData();  
  return 0;
}
