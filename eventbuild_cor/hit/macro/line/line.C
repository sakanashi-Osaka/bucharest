#include <fstream>
#include <iostream>
#include <string>

int line(string name){
  double ene[101]={}; 
  double theta[101]={}; 
  
  std::ifstream ifs(name);
  for(int i=0; i<101; i++){
    ifs >> theta[i] >> ene[i];
    cout << theta[i] << " " << ene[i] << endl;
  }

  TGraph *g = new TGraph(101,theta,ene);
  g->Draw();
  
  return 0;
}
