#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <TCanvas.h>
#include <TEventList.h>
#include <cstdint>
#include <vector>
#include <signal.h>
#include <map>
#include <algorithm>
#include <regex>
#include <utility>
#include <TH1F.h>
#include <TH2.h>
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>

#define PI 3.14159265
#define N_BOARD 10
#define N_CH 16
#define N_HIT_MAX 10
#define DtoR 0.0174533

const float s1_r1 = 48/2.0;  // inner radius      
const float s1_r2 = 96/2.0;  // outer radius
const float s1_dist = 40;

float get_front_theta(int front_ch);
std::pair<float,float> get_angle(float x, float y, int ch_f, int ch_r);

//analysis
double enea[3]={4.2,3,3.4};
int chfa[3]={5,6,12};
int chra[3]={9,13,13};

double gomi_culc(){

  double m4He=931.494*4+2.4249;
  double m12C=931.494*12;
  double tha[3]={};
  double pha[3]={};
  double moma[3][3]={};
  double mom12C[3]={0,0,0};
  double total12C=-1;
  double ex12C=-1;
  
  for(int n=0; n<3; n++){
    tha[n] = get_angle(0,2,chfa[n],chra[n]).first;
    pha[n] = get_angle(0,2,chfa[n],chra[n]).second;
    //    cout << tha[n] << " " << pha[n] << endl;
    //    cout << tha[n]*180/3.14 << " " << pha[n]*180/3.14 << endl;
    
    moma[n][0] = pow(2*m4He*enea[n]+enea[n]*enea[n],0.5) * sin(tha[n]) * cos(pha[n]);
    moma[n][1] = pow(2*m4He*enea[n]+enea[n]*enea[n],0.5) * sin(tha[n]) * sin(pha[n]);
    moma[n][2] = pow(2*m4He*enea[n]+enea[n]*enea[n],0.5) * cos(tha[n]);
    //    cout << enea[n] << " " << moma[n][0]  << " " << moma[n][1] << " "  << moma[n][2] << endl;
    //    cout << m4He+enea[n] << " " << moma[n][0]  << " " << moma[n][1] << " "  << moma[n][2] << endl;

    mom12C[0] += moma[n][0];
    mom12C[1] += moma[n][1];
    mom12C[2] += moma[n][2];
  }
  //  cout << mom12C[0]  << " " << mom12C[1] << " "  << mom12C[2] << endl;

  total12C = pow(mom12C[0]*mom12C[0]+mom12C[1]*mom12C[1]+mom12C[2]*mom12C[2]+m12C*m12C,0.5);
  ex12C = enea[0]+enea[1]+enea[2]+m4He*3 - total12C;
  
  return ex12C;
}



float get_front_theta(int front_ch){
  
  float strp_wid = (s1_r2 - s1_r1)/16.0;
  float r=0;  
  r = s1_r1 + strp_wid*(0.5) + strp_wid*(15-front_ch);
  
  return atan(r/s1_dist);
}

std::pair<float,float> get_angle(float x, float y, int ch_f, int ch_r){

  float tmp_x, tmp_y, tmp_z, tmp_theta, tmp_phi;
  float front_theta = get_front_theta(ch_f);
  float rear_phi = (90 + 22.5/2.0 + 22.5*ch_r)*PI/180;

  tmp_x = s1_dist * tan(front_theta)*cos(rear_phi) - x;
  tmp_y = s1_dist * tan(front_theta)*sin(rear_phi) - y;
  tmp_z = s1_dist;
  tmp_theta = acos(tmp_z/pow(pow(tmp_x,2)+pow(tmp_y,2)+pow(tmp_z,2),0.5));
  if(tmp_y>=0) tmp_phi = acos(tmp_x/pow(pow(tmp_x,2)+pow(tmp_y,2),0.5));
  if(tmp_y<0) tmp_phi = -1 * acos(tmp_x/pow(pow(tmp_x,2)+pow(tmp_y,2),0.5));
  
  return std::make_pair(tmp_theta,tmp_phi);
}




