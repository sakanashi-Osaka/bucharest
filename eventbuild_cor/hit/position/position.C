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
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>

#define PI 3.14159265
#define N_BOARD 10
#define N_CH 16
#define N_HIT_MAX 20

using namespace std;

const double s1_r1 = 48/2.0;  // inner radius      
const double s1_r2 = 96/2.0;  // outer radius
const double s1_dist = 40;

double get_front_theta(int front_ch);
std::pair<double,double> get_angle(double x, double y, int ch_f, int ch_r);
double get_front_ex(double theta, double ene);
double detect3a(double *enea, int *chfa, int *chra);


int hit_n[4]={0,0,0,0};
int ch_f[4][N_HIT_MAX]={};
int ch_r[4][N_HIT_MAX]={};
double Energy_f[4][N_HIT_MAX]={};
double Energy_r[4][N_HIT_MAX]={};

double cor_ex[4][40][40]={};
int ch_F[4]={};
int ch_R[4]={};


int position(int run){
  
  ofstream ofs(Form("log/log%d.txt",run),std::ios::app);
  
  TFile *fin =new TFile(Form("../rootfile/hit%d.root",run));
  TTree *hit = (TTree*)fin->Get("hit");
  
  TFile *fout =new TFile(Form("rootfile/position%d.root",run),"recreate");
  
  TH2D *matrix = new TH2D("matrix","matrix",40,0,40,40,0,40);
  TH1D *matrix_x = new TH1D("matrix_x","matrix_x",40,0,40);
  TH1D *matrix_y = new TH1D("matrix_y","matrix_y",40,0,40);
    
  TH1D *h[3][16][40][40]; //[3]...chf(1,8,14), [8]...ch_r, [40]...dx, [40]...dy
  for(Int_t i=0;i<3;i++){
    for(Int_t z=0;z<16;z++){
      for(Int_t j=0;j<40;j++){
	for(Int_t k=0;k<40;k++){   
	  h[i][z][j][k] = new TH1D(Form("h_%d_%d_%d_%d",i,z,j,k),Form("h_%d_%d_%d_%d",i,z,j,k),100,-2,2);
	}
      }
    }
  }
  
  hit->SetBranchAddress("hit_n",hit_n);
  hit->SetBranchAddress("ch_f",ch_f);
  hit->SetBranchAddress("ch_r",ch_r);
  hit->SetBranchAddress("Energy_f",Energy_f);
  hit->SetBranchAddress("Energy_r",Energy_r);
  
  
  ULong64_t N=hit->GetEntries();
  cout << "Total entry: " << N << endl;
  //  for(ULong64_t evtn=0; evtn<N; evtn++){
  for(ULong64_t evtn=0; evtn<100000; evtn++){
    if(evtn%1000==0) cout << "\rAnalyzed entry:" << evtn; std::cout << flush;
    
    for(int seg=0; seg<4; seg++){
      for(int l=0; l<40; l++){
	for(int m=0; m<40; m++){
	  cor_ex[seg][l][m] = -10;
	}
      }
    }
    
    hit->GetEntry(evtn); 

    for(int seg=0; seg<4; seg++){
      
      if(hit_n[seg]==11 && Energy_f[seg][0]>18){
	for(int l=0; l<40; l++){
	  for(int m=0; m<40; m++){
	    double dx = ((double)l-20)/5;
	    double dy = ((double)m-20)/5;
	    double tmp_theta1 = get_angle(dx, dy, ch_f[seg][0], ch_r[seg][0]).first;
	    cor_ex[seg][l][m] = get_front_ex(tmp_theta1, (double)Energy_f[seg][0]);

	    ch_F[seg]=ch_f[seg][0];
	    ch_R[seg]=ch_r[seg][0];

	    if(ch_F[seg]==1) h[0][ch_R[seg]][l][m]->Fill((double)cor_ex[seg][l][m]);
	    if(ch_F[seg]==8) h[1][ch_R[seg]][l][m]->Fill((double)cor_ex[seg][l][m]);
	    if(ch_F[seg]==14) h[2][ch_R[seg]][l][m]->Fill((double)cor_ex[seg][l][m]);
	  }
	}
      }
    } 
  }

  for(Int_t i=0;i<3;i++){
    for(Int_t z=0;z<16;z++){
      for(Int_t j=0;j<40;j++){
	for(Int_t k=0;k<40;k++){   
	  h[i][z][j][k]->Write();
	}
      }
    }
  }

  TSpectrum *s[3][16][40][40];
  double *xpeaks[3][16][40][40];
  double peak[3][16][40][40];
  double diff[3][4][40][40]={};
  double sum_diff[3][40][40]={};
  double sum_sum_diff[40][40]={};

  for(Int_t i=0;i<3;i++){
    for(Int_t z=0;z<16;z++){
      for(Int_t j=0;j<40;j++){
	for(Int_t k=0;k<40;k++){   
	  diff[i][z][j][k]=0;
	  sum_diff[i][j][k]=0;
	  sum_sum_diff[j][k]=0;
	}
      }
    }
  }
	  
  for(Int_t i=0;i<3;i++){
    for(Int_t z=0;z<16;z++){
      for(Int_t j=0;j<40;j++){
	for(Int_t k=0;k<40;k++){   
	  s[i][z][j][k] = new TSpectrum(1);
	  s[i][z][j][k]->Search(h[i][z][j][k],1,"n new",0.2);
	  xpeaks[i][z][j][k]=s[i][z][j][k]->GetPositionX();
	  peak[i][z][j][k]=xpeaks[i][z][j][k][0];
	}
      }
    }
  }
  
  for(Int_t i=0;i<3;i++){
    for(Int_t z=0;z<4;z++){
      for(Int_t j=0;j<40;j++){
	for(Int_t k=0;k<40;k++){
	  diff[i][z][j][k] = abs( peak[i][4*z+3][j][k] - peak[i][4*z+0][j][k] ); 
	}
      }
    }
  }
  
  for(Int_t i=0;i<3;i++){
    for(Int_t z=0;z<4;z++){
      for(Int_t j=0;j<40;j++){
	for(Int_t k=0;k<40;k++){
	  sum_diff[i][j][k] += diff[i][z][j][k]; 
	}
      }
    }
  }
  
  for(Int_t i=0;i<3;i++){
    for(Int_t j=0;j<40;j++){
      for(Int_t k=0;k<40;k++){
	sum_sum_diff[j][k] += sum_diff[i][j][k]; 
      }
    }
  }
  
  for(Int_t j=0;j<40;j++){
    for(Int_t k=0;k<40;k++){
      ofs << sum_sum_diff[j][k] << " ";
      for(Int_t tmp=0; tmp<(int)(sum_sum_diff[j][k]*100); tmp++){
	matrix->Fill(j,k);
	matrix_x->Fill(j);
	matrix_y->Fill(k);
      }
    }
    ofs << endl;
  }


  matrix->Write();
  matrix_x->Write();
  matrix_y->Write();
  fout->Close();

  return 0;
}




double get_front_theta(int front_ch){
  
  double strp_wid = (s1_r2 - s1_r1)/16.0;
  double r=0;  
  r = s1_r1 + strp_wid*(0.5) + strp_wid*(15-front_ch);
  
  return atan(r/s1_dist);
}

double get_front_ex(double theta, double ene){
  const double beam_ene = 25.0;
  const double AMU = 931.4943;
  const double mass_ex_12c = 0;
  const double mass_ex_4he = 2.425;
  const double mass_ex_13c = 3.125;    

  double mass_12c = AMU*12 + mass_ex_12c;
  double mass_4he = AMU*4  + mass_ex_4he;
  double mass_13c = AMU*13 + mass_ex_13c;  
  
  double m1, m2, m3, m4;
  double E1, E3;
  double p1, p3;

  double s,t,u;
  double total_m4;
  double ex4;
  
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

std::pair<double,double> get_angle(double x, double y, int ch_f, int ch_r){

  double tmp_x, tmp_y, tmp_z, tmp_theta, tmp_phi;
  double front_theta = get_front_theta(ch_f);
  double rear_phi = (90 + 22.5/2.0 + 22.5*ch_r)*PI/180;

  tmp_x = s1_dist * tan(front_theta)*cos(rear_phi) - x;
  tmp_y = s1_dist * tan(front_theta)*sin(rear_phi) - y;
  tmp_z = s1_dist;
  tmp_theta = acos(tmp_z/pow(pow(tmp_x,2)+pow(tmp_y,2)+pow(tmp_z,2),0.5));
  if(tmp_y>=0) tmp_phi = acos(tmp_x/pow(pow(tmp_x,2)+pow(tmp_y,2),0.5));
  if(tmp_y<0) tmp_phi = -1 * acos(tmp_x/pow(pow(tmp_x,2)+pow(tmp_y,2),0.5));
  
  return std::make_pair(tmp_theta,tmp_phi);
}




