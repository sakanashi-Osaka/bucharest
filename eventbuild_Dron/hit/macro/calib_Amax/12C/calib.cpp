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
#include <TGraph.h>
#include <TF1.h>
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
  
  ofstream ofs(Form("gomi%d.txt",atoi(argv[1])),std::ios::app);
  
  TFile *fin;
  TTree *tree;
  TH1D *hAmax2[16][16];
  TH1D *hAmax3[16][16];
  TH1D *hAmax4[16][16];
  TH1D *hAmax5[16][16];
  
  fin = new TFile(Form("../../rootfile/gomi%d.root",atoi(argv[1])));
  tree = (TTree*)fin->Get("tree");

  for(int i=0; i<16; i++){
    for(int j=0; j<16; j++){
      hAmax2[i][j] = new TH1D(Form("hAmax2_%d_%d",i,j),Form("hAmax2_%d_%d",i,j),100,-0.2,0.5);
      hAmax3[i][j] = new TH1D(Form("hAmax3_%d_%d",i,j),Form("hAmax3_%d_%d",i,j),100,-0.2,0.5);
      hAmax4[i][j] = new TH1D(Form("hAmax4_%d_%d",i,j),Form("hAmax4_%d_%d",i,j),100,-0.2,0.5);
      hAmax5[i][j] = new TH1D(Form("hAmax5_%d_%d",i,j),Form("hAmax5_%d_%d",i,j),100,-0.2,0.5);
    }
  }

  double chf12C;
  double chr12C;
  double K12C;
  double Amax12C;

  tree->SetBranchAddress("chf12C",&chf12C);
  tree->SetBranchAddress("chr12C",&chr12C);
  tree->SetBranchAddress("K12C",&K12C);
  tree->SetBranchAddress("Amax12C",&Amax12C);

  ULong64_t N=tree->GetEntries();
  cout << "Total entry: " << N << endl;
  for(ULong64_t evtn=0; evtn<N; evtn++){
  //  for(ULong64_t evtn=0; evtn<1e7; evtn++){
    if(evtn%100000==0) cout << "\rAnalyzed entry:" << evtn; std::cout << flush;

    tree->GetEntry(evtn);
    if(chf12C>-1 && chr12C>-1 && Amax12C>-0.2){
      if(abs(K12C-2)<0.2) hAmax2[(int)chf12C][(int)chr12C]->Fill(Amax12C);
      if(abs(K12C-3)<0.2) hAmax3[(int)chf12C][(int)chr12C]->Fill(Amax12C);
      if(abs(K12C-4)<0.2) hAmax4[(int)chf12C][(int)chr12C]->Fill(Amax12C);
      if(abs(K12C-5)<0.2) hAmax5[(int)chf12C][(int)chr12C]->Fill(Amax12C);
    }
  }

  double m2[16][16];
  double m3[16][16];
  double m4[16][16];
  double m5[16][16];
  for(int i=0; i<16; i++){
    for(int j=0; j<16; j++){
      m2[i][j] = -0.2+0.7/100*hAmax2[i][j]->GetMaximumBin();
      m3[i][j] = -0.2+0.7/100*hAmax3[i][j]->GetMaximumBin();
      m4[i][j] = -0.2+0.7/100*hAmax4[i][j]->GetMaximumBin();
      m5[i][j] = -0.2+0.7/100*hAmax5[i][j]->GetMaximumBin();
    }
  }

  double x[16][16][4];
  double y[16][16][4];
  for(int i=0; i<16; i++){
    for(int j=0; j<16; j++){
      x[i][j][0]=2;
      x[i][j][1]=3;
      x[i][j][2]=4;
      x[i][j][3]=5;
      y[i][j][0]=m2[i][j];
      y[i][j][1]=m3[i][j];
      y[i][j][2]=m4[i][j];
      y[i][j][3]=m5[i][j];
    }
  }

  TGraph *g[16][16];
  TF1 *f[16][16];
  for(int i=0; i<16; i++){
    for(int j=0; j<16; j++){
      g[i][j] = new TGraph(4,x[i][j],y[i][j]);
      g[i][j]->SetMarkerStyle(2);
      g[i][j]->Draw("AP");
      f[i][j] = new TF1(Form("f_%d_%d",i,j),"[0] + [1]*x + [2]/(x-[3])",2,5);
      g[i][j]->Fit(f[i][j]);
      double p0,p1,p2,p3,p4;
      p0 = f[i][j]->GetParameter(0);
      p1 = f[i][j]->GetParameter(1);
      p2 = f[i][j]->GetParameter(2);
      p3 = f[i][j]->GetParameter(3);
      ofs << p0 << " " << p1 << " " << p2 << " " << p3 << endl;
    }
  }

  std::cout << argv[1] << "finished!" << endl;
  fin->Close();

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




