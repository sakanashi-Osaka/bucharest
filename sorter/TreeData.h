#ifndef TreeData_hpp
#define TreeData_hpp 1

#include <vector>

//#include <CAENDigitizerType.h>

class TreeData
{ // no getter setter.  using public member variables.
public:
  TreeData(){};

  TreeData(uint32_t nSamples)
  {};

  ~TreeData(){};

  int domain;
  ULong64_t TS; // for readinng tree
  double FineTS;
  //  double TDC;
  int ADC;
  float Amax;
  float Energy;
  uint32_t Extras;

  static const uint16_t OneHitSize = sizeof(domain) + sizeof(FineTS) + sizeof(ADC) + sizeof(Energy) + sizeof(Amax);
};
typedef TreeData TreeData_t;

#endif
