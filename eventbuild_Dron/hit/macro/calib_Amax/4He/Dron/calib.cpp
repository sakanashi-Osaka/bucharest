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
#include <TH1F.h>
#include "TF1.h"
#include <TH2.h>
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>

#define PI 3.14159265
#define N_BOARD 10
#define N_CH 16
#define N_HIT_MAX 10

using namespace std;

const float s1_r1 = 48/2.0;  // inner radius      
const float s1_r2 = 96/2.0;  // outer radius
const float s1_dist = 40;

float get_front_theta(int front_ch);
std::pair<float,float> get_angle(float x, float y, int ch_f, int ch_r);
float get_front_ex(float theta, float ene);
double detect3a(double *enea, int *chfa, int *chra);


int main(int argc, char *argv[]){
  //int calib(int run, int fch, int rch){

  int run = atoi(argv[1]);
  int fch = atoi(argv[2]);
  int rch = atoi(argv[3]);
  
  ofstream ofs(Form("test%d.txt",run),std::ios::app);
  
  TFile *fin =new TFile(Form("../../../rootfile/Dron%d.root",run));
  TTree *tree = (TTree*)fin->Get("tree");

  TH1D *hAmax = new TH1D("hAmax","hAmax",100,-0.05,0.10);

  double chf4He;
  double chr4He;
  double Ex4He;
  double Amax4He;

  tree->SetBranchAddress("chf4He",&chf4He);
  tree->SetBranchAddress("chr4He",&chr4He);
  tree->SetBranchAddress("Ex4He",&Ex4He);
  tree->SetBranchAddress("Amax4He",&Amax4He);

  ULong64_t N=tree->GetEntries();
  cout << "Total entry: " << N << endl;
  for(ULong64_t evtn=0; evtn<N; evtn++){
    if(evtn%100000==0) cout << "\rAnalyzed entry:" << evtn; std::cout << flush;

    tree->GetEntry(evtn);
    if(chf4He==fch && chr4He==rch && Amax4He>-0.05){
      if(abs(Ex4He-9.64)<0.25) hAmax->Fill(Amax4He);
    }
  }

  int nbin = hAmax->GetMaximumBin();
  double m = -0.05+0.15/100*(double)nbin;
  int n = hAmax->GetBinContent(nbin);

  TF1 *f = new TF1("f","[0]*exp(-pow(x-[1],2)/pow([2],2)/2)",-0.05,0.1);
  f->SetParameter(0,n);
  f->SetParameter(1,m);
  f->SetParameter(2,0.005);

  
  hAmax->Fit(f,"","",m-0.007,m+0.005);
 
  double p0 = f->GetParameter(0);
  double p1 = f->GetParameter(1);
  double p2 = f->GetParameter(2);

  ofs << fch << " " << rch << " " << p0 << " " << p1 << " " << p2 << endl;
  
  return 0;
}  




float get_front_theta(int front_ch){
  
  float strp_wid = (s1_r2 - s1_r1)/16.0;
  float r=0;  
  r = s1_r1 + strp_wid*(0.5) + strp_wid*(15-front_ch);
  
  return atan(r/s1_dist);
}

float get_front_ex(float theta, float ene){
  const float beam_ene = 25.0;
  const float AMU = 931.4943;
  const float mass_ex_12c = 0;
  const float mass_ex_4he = 2.425;
  const float mass_ex_13c = 3.125;    

  float mass_12c = AMU*12 + mass_ex_12c;
  float mass_4he = AMU*4  + mass_ex_4he;
  float mass_13c = AMU*13 + mass_ex_13c;  
  
  float m1, m2, m3, m4;
  float E1, E3;
  float p1, p3;

  float s,t,u;
  float total_m4;
  float ex4;
  
  m1 = mass_4he;
  m2 = mass_12c;
  m3 = mass_4he;
  m4 = mass_12c;

  E1 = m1 + beam_ene;
  E3 = m3 + ene;  

  p1 = sqrt(E1*E1-m1*m1);
  p3 = sqrt(E3*E3-m3*m3);  
  
  s = m1*m1 + m2*m2 + 2*m2*E1;
  t = m1*m1 + m3*m3 + 2*(p1*p3*cos(theta) - E1*E3);
  u = m2*m2 + m3*m3 - 2*m2*E3;

  total_m4 = sqrt(s+t+u - m1*m1 -m2*m2 - m3*m3);
  ex4 = total_m4 - m4;

  return ex4;
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




