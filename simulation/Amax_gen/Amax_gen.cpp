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

int Amax_gen(){
   
  ofstream ofs("prm/prm_Amax_2305.txt",std::ios::app);
  ofstream ofs1("prm/fit_Amax_mean2305.txt",std::ios::app);
  ofstream ofs2("prm/fit_Amax_sigma2305.txt",std::ios::app);
  
  TFile *fin;
  TTree *tree;
  TH1D *hAmax[4][16][16];
  TF1 *fAmax[4][16][16];
  
  fin = new TFile("rootfile/Amax_no_cor2305.root");
  tree = (TTree*)fin->Get("hit");
  
  for(int i=0; i<4; i++){
    for(int j=0; j<16; j++){
      for(int k=0; k<16; k++){
	hAmax[i][j][k] = new TH1D(Form("hAmax_%d_%d_%d",i,j,k),Form("hAmax_%d_%d_%d",i,j,k),100,0.2,0.7);
	fAmax[i][j][k] = new TF1(Form("fAmax_%d_%d_%d",i,j,k),"[0]*exp(-pow(x-[1],2)/(2*pow([2],2)))",0.2,0.7);
      }
    }
  }
  
  double Amax[4][10];
  double Energy_f[4][10];
  int ch_f[4][10];
  int ch_r[4][10];
  int hit_n[4];
	       
  tree->SetBranchAddress("Amax",&Amax);
  tree->SetBranchAddress("Energy_f",&Energy_f);
  tree->SetBranchAddress("ch_f",&ch_f);
  tree->SetBranchAddress("ch_r",&ch_r);
  tree->SetBranchAddress("hit_n",&hit_n);

  TFile *fout = new TFile("rootfile/gomi.root","recreate");
  
  ULong64_t N=tree->GetEntries();
  cout << "Total entry: " << N << endl;
  for(ULong64_t evtn=0; evtn<N; evtn++){
  //  for(ULong64_t evtn=0; evtn<5e6; evtn++){
    if(evtn%100000==0) cout << "\rAnalyzed entry:" << evtn; std::cout << flush;

    tree->GetEntry(evtn);
    for(int j=0; j<4; j++){
      if(ch_f[j][0]>-1 && ch_r[j][0]>-1 && Amax[j][0]>0.2 && hit_n[j]==11){
	if(abs(Energy_f[j][0]-1.5)<0.15) hAmax[0][ch_f[j][0]][ch_r[j][0]]->Fill(Amax[j][0]);
	if(abs(Energy_f[j][0]-2.5)<0.15) hAmax[1][ch_f[j][0]][ch_r[j][0]]->Fill(Amax[j][0]);
	if(abs(Energy_f[j][0]-3.5)<0.15) hAmax[2][ch_f[j][0]][ch_r[j][0]]->Fill(Amax[j][0]);
	if(abs(Energy_f[j][0]-4.5)<0.15) hAmax[3][ch_f[j][0]][ch_r[j][0]]->Fill(Amax[j][0]);
      }
    }
  }

  int bin[4][16][16];
  int n[4][16][16];
  double p0[4][16][16];
  double p1[4][16][16];
  double p2[4][16][16];
  for(int i=0; i<4; i++){
    for(int j=0; j<16; j++){
      for(int k=0; k<16; k++){
        bin[i][j][k] = hAmax[i][j][k]->GetMaximumBin();
	n[i][j][k] = hAmax[i][j][k]->GetBinContent(bin[i][j][k]);

	fAmax[i][j][k]->SetParameter(0,n[i][j][k]);
	fAmax[i][j][k]->SetParameter(1,0.2+(double)bin[i][j][k]*0.005);
	fAmax[i][j][k]->SetParameter(2,0.01);
	fAmax[i][j][k]->SetParLimits(0,(double)n[i][j][k]*0.5,n[i][j][k]*2);
	fAmax[i][j][k]->SetParLimits(1,0.2+(double)bin[i][j][k]*0.005-0.05,0.2+(double)bin[i][j][k]*0.005+0.05);
	fAmax[i][j][k]->SetParLimits(2,0.005,0.05);

	hAmax[i][j][k]->Fit(fAmax[i][j][k],"","",0.2,0.7);
	p0[i][j][k] = fAmax[i][j][k]->GetParameter(0);
	p1[i][j][k] = fAmax[i][j][k]->GetParameter(1);
	p2[i][j][k] = fAmax[i][j][k]->GetParameter(2);
	ofs << i << " " << j << " " << k << " " << p0[i][j][k] << " " << p1[i][j][k] << " " << p2[i][j][k] <<  " " << endl;
      }
    }
  }
  
  
  double x[16][16][4];
  double y[16][16][4];
  for(int i=0; i<16; i++){
    for(int j=0; j<16; j++){
      x[i][j][0]=1.5;
      x[i][j][1]=2.5;
      x[i][j][2]=3.5;
      x[i][j][3]=4.5;
      y[i][j][0]=p1[0][i][j];
      y[i][j][1]=p1[1][i][j];
      y[i][j][2]=p1[2][i][j];
      y[i][j][3]=p1[3][i][j];
    }
  }
  TGraph *g[16][16];
  TF1 *f[16][16];
  for(int i=0; i<16; i++){
    for(int j=0; j<16; j++){
      g[i][j] = new TGraph(4,x[i][j],y[i][j]);
      g[i][j]->SetMarkerStyle(2);
      g[i][j]->Draw("AP");
      f[i][j] = new TF1(Form("f_%d_%d",i,j),"[0] + [1]*x + [2]/(x-[3])",1.5,4.5);
      g[i][j]->Fit(f[i][j]);
      double p0,p1,p2,p3,p4;
      p0 = f[i][j]->GetParameter(0);
      p1 = f[i][j]->GetParameter(1);
      p2 = f[i][j]->GetParameter(2);
      p3 = f[i][j]->GetParameter(3);
      ofs1 << i << " " << j << " " << p0 << " " << p1 << " " << p2 << " " << p3 << " " << endl;
    }
  }

  double x1[16][16][4];
  double y1[16][16][4];
  for(int i=0; i<16; i++){
    for(int j=0; j<16; j++){
      x1[i][j][0]=1.5;
      x1[i][j][1]=2.5;
      x1[i][j][2]=3.5;
      x1[i][j][3]=4.5;
      y1[i][j][0]=p2[0][i][j];
      y1[i][j][1]=p2[1][i][j];
      y1[i][j][2]=p2[2][i][j];
      y1[i][j][3]=p2[3][i][j];
    }
  }
  TGraph *g1[16][16];
  TF1 *f1[16][16];
  for(int i=0; i<16; i++){
    for(int j=0; j<16; j++){
      g1[i][j] = new TGraph(4,x1[i][j],y1[i][j]);
      g1[i][j]->SetMarkerStyle(2);
      g1[i][j]->Draw("AP");
      f1[i][j] = new TF1(Form("f_%d_%d",i,j),"[0] + [1]*x + [2]/(x-[3])",1.5,4.5);
      g1[i][j]->Fit(f1[i][j]);
      double p0,p1,p2,p3,p4;
      p0 = f1[i][j]->GetParameter(0);
      p1 = f1[i][j]->GetParameter(1);
      p2 = f1[i][j]->GetParameter(2);
      p3 = f1[i][j]->GetParameter(3);
      ofs2 << i << " " << j << " " << p0 << " " << p1 << " " << p2 << " " << p3 << endl;
    }
  }


  //  fout->Close();

  return 0;
}  


