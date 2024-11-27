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



int calib(int run, int fch, int rch){
  
  ofstream ofs("test.txt",std::ios::app);
  
  TFile *fin =new TFile(Form("../rootfile/test%d.root",run));
  TTree *tree = (TTree*)fin->Get("tree");

  TH1D *htsr = new TH1D("htsr","htsr",100,0,100);
  TH1D *htsf = new TH1D("htsf","htsf",100,0,100);

  double chf12C;
  double chr12C;
  double tsf12C;
  double tsr12C;

  tree->SetBranchAddress("chf12C",&chf12C);
  tree->SetBranchAddress("chr12C",&chr12C);
  tree->SetBranchAddress("tsf12C",&tsf12C);
  tree->SetBranchAddress("tsr12C",&tsr12C);

  ULong64_t N=tree->GetEntries();
  cout << "Total entry: " << N << endl;
  for(ULong64_t evtn=0; evtn<N; evtn++){
    if(evtn%100000==0) cout << "\rAnalyzed entry:" << evtn; std::cout << flush;

    tree->GetEntry(evtn);
    if(chf12C==fch && chr12C==rch){
      htsr->Fill(tsr12C);
      htsf->Fill(tsf12C);
    }
  }

  int mr = htsr->GetMaximumBin();
  int mf = htsf->GetMaximumBin();
  int nr = htsr->GetBinContent(mr);
  int nf = htsf->GetBinContent(mf);

  TF1 *fr = new TF1("fr","[0]*exp(-pow(x-[1],2)/(2*pow([2],2)))",mr-10,mr+10);
  fr->SetParameter(0,nr);
  fr->SetParameter(1,mr);
  fr->SetParameter(2,5);
  htsr->Fit(fr,"","",mr-10,mr+10);

  TF1 *ff = new TF1("ff","[0]*exp(-pow(x-[1],2)/(2*pow([2],2)))",mr-10,mr+10);
  ff->SetParameter(0,nf);
  cout << endl;
  cout << nf << endl;
  ff->SetParameter(1,mf);
  ff->SetParameter(2,5);
  htsf->Fit(ff,"","",mr-10,mr+10);
  
  double pr2 = fr->GetParameter(1);
  double pf2 = ff->GetParameter(1);
  
  ofs << fch << " " << rch << " " << pr2 << " " << pf2 << endl;
 
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




