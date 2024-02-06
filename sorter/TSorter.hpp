#ifndef TSorter_hpp
#define TSorter_hpp 1

#include <vector>

#include <TString.h>

#include "TreeData.h"


class TSorter
{
public:
  TSorter();
  ~TSorter();

  void SetServerName(TString server) { fServerName = server; };
  void SetInputDirectory(TString inputDir) { fInputDir = inputDir; };
  
  // When the stop < 0, read all files more than start
  void LoadFiles(Int_t runNumber, Int_t start = 0, Int_t stop = -1);
  void SortData();
  // void WriteMultipleFile();
  // void WriteOneFile(TString fileName = "");
  void WriteData(TString fileName = "");

  void TestOffset();// same as WriteData.  But subtract time offset from Si
  
private:
  std::vector<TreeData *> fDataVec;
  TString fServerName = "ssgant1";
  TString fInputDir = "/home/sakra/data/Bucharest2022/sorter";
  //  TString fServerName = "";
  //  TString fInputDir = "/home/sakra/data/Bucharest2022/data";
  Int_t fRunNumber;
  
};


#endif
