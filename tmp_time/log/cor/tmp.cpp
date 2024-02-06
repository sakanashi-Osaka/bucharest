#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

int main(int argc, char* argv[]){
  std::ifstream ifs(argv[1]);

  std::vector< std::vector<double> > array(101);
  std::vector< std::vector<double> > time(101);
  std::string str = "";
  int i = 0;

  while(getline(ifs, str)){
    std::string tmp = "";
    std::istringstream stream(str);
    while(getline(stream, tmp, ' ')){
      array.at(i).push_back(std::stod(tmp));
    }
    i++;
  }

  for(int j=1; j<array.at(0).size(); j++){
    std::cout << array.at(0).at(j) << " ";
  }
  std::cout << std::endl;

  for(int i=0; i<array.size(); i++){
    std::cout << 64+i-2 << " " << (i+14)%16 << " ";
    for(int j=1; j<array.at(i).size(); j++){
      if(array.at(i).size()!=array.at(0).size()) continue;
      time.at(i).push_back(array.at(i).at(j)-array.at(0).at(j));
      //      std::cout << time.at(i).at(j-1) << " ";
      
      bool flag = false;
      for(int k=-100; k<100; k++){
        if(std::abs((double)time.at(i).at(j-1)-2.147483*(double)k)<1.01){
	  std::cout << k << " ";
          flag = true;
        }
      }
      if(flag==false) std::cout << "error" << " ";

    }
    std::cout << std::endl;
  }

  return 0;
}
