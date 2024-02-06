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

#define N_BOARD 10
#define N_CH 16

int entry_type(int board, int ch);



int tmp_time(int run){
  
  TFile *fin =new TFile(Form("../sorter_cor/rootfile/run%d_-1.root",run));
  TTree *tree = (TTree*)fin->Get("tree");

  std::ofstream ofs(Form("log/log%d.txt",run),std::ios::app);
  
  int tmp_domain;
  double tmp_ts;  
  tree->SetBranchAddress("domain", &tmp_domain);
  tree->SetBranchAddress("FineTS", &tmp_ts);  

  
  double last;
  double ave[100]={};
  int ave_n[100]={};
  double ave51=0;
  int ave_n51=0;
  int count[100][8000]={}; //ch63-162, 0s-4000s(1/2s dev.)
  int count51[8000]={}; //ch51, 0s-4000s(1/2s dev.)
  for(int i=0; i<100; i++){
    ave[i]=0;
    ave_n[i]=0;
    for(int j=0; j<8000; j++){
      count[i][j]=0;
    }
  }
  ULong64_t N=tree->GetEntries();
  cout << "Total entry: " << N << endl;
  for(ULong64_t evtn=0; evtn<N; evtn++){
  //  for(ULong64_t evtn=0; evtn<2e8; evtn++){
    if(evtn%100000==0) cout << "\rAnalyzed entry:" << evtn; std::cout << flush;

    tree->GetEntry(evtn);
    last = tmp_ts/5e11;
    int bin = (int)(tmp_ts/5e11);
    if(tmp_domain<51 || tmp_domain>=163 || bin<0 || bin>=8000) continue;
    if(tmp_domain>62) count[tmp_domain-63][bin]+=1;
    if(tmp_domain==51) count51[bin]+=1;
  }

  for(int i=0; i<100; i++){
    for(int j=0; j<8000; j++){
      if(count[i][j]>0){
	ave[i]+=count[i][j];
	ave51+=count51[j];
	ave_n[i]++;
	ave_n51++;
      }
    }
    ave[i]=ave[i]/ave_n[i];
    cout << i << " " << ave[i] << endl;
  }
  ave51=ave51/ave_n51;
  cout << "51 " << ave51 << endl;


  ofs << 51 << " ";
  for(int j=2; j<8000-10; j++){
    if(count51[j-2]<ave51/3 && count51[j-1]<ave51/3 && count51[j]<ave51/3 && count51[j+1]>ave51/3 && count51[j+2]>ave51/3 && count51[j+3]>ave51/3 && count51[j+4]>ave51/3 && count51[j+5]>ave51/3){
      ofs << (double)j/2 << " ";
    }
  }
  ofs << endl;
	  
  for(int i=0; i<100; i++){
    ofs << i+63 << " ";
    for(int j=1; j<8000-10; j++){
      if(count[i][j]<ave[i]/5 && count[i][j+1]>ave[i]/5 && count[i][j+2]>ave[i]/5 && count[i][j+3]>ave[i]/5 && count[i][j+4]>ave[i]/5 && count[i][j+5]>ave[i]/5 && count[i][j+6]>ave[i]/5 && count[i][j+7]>ave[i]/5 && count[i][j+8]>ave[i]/5){
	  ofs << (double)j/2 << " ";
      }
    }
    ofs << endl;
  }


  
  return 0;
}



int entry_type(int board, int ch){
  int type_n=-1; // 0:trigger, 1:Gamma_signal, 2:Si_f_signal, 3:Si_r_signal, 4:scaler

  if(board==0 && ch==0) type_n = 0;
  if(board==0 && ch>1) type_n = 1;
  if(board==1 || board==2 || board==3){
    if(ch<12) type_n = 1;
    if(ch>11) type_n = 4;
  }
  if(board==4 || board==5 || board==6 || board==7){
    //    if(ch<1) type_n = 0;
    if(ch>1) type_n = 2;
  }
  //  if(board==8 && ch<1) type_n = 0;
  if(board==8 && ch>1 && ch<10) type_n = 2;
  if(board==8 && ch>9) type_n = 3;
  //  if(board==9 && ch<1) type_n = 0;
  if(board==9 && ch>1 && ch<10) type_n = 3;  
  
  return type_n;
}; 

