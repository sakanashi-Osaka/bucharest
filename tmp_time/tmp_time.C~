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
  
  TFile *fin =new TFile(Form("../sorter_cor/rootfile/marged%d_-1.root",run));
  TTree *tree = (TTree*)fin->Get("tree");

  std::ofstream ofs(Form("log/log%d.txt",run),std::ios::app);
  
  int tmp_domain;
  double tmp_ts;  
  tree->SetBranchAddress("domain", &tmp_domain);
  tree->SetBranchAddress("FineTS", &tmp_ts);  

  double last;
  int count[160][8000]={};
  double ave[160]={};
  for(int i=0; i<160; i++){
    ave[i]=0;
    for(int j=0; j<8000; j++){
      count[i][j]=0;
    }
  }
  
  ULong64_t N=tree->GetEntries();
  cout << "Total entry: " << N << endl;
  for(ULong64_t evtn=0; evtn<N; evtn++){
  //    for(ULong64_t evtn=0; evtn<0.5e8; evtn++){
    if(evtn%10000==0) cout << "\rAnalyzed entry:" << evtn; std::cout << flush;

    tree->GetEntry(evtn);
    last = tmp_ts/5e11;
    int bin = (int)(tmp_ts/5e11);
    if(tmp_domain<0 || tmp_domain>=160 || bin<0 || bin>=8000) continue;
    count[tmp_domain][bin]+=1;
  }

  for(int i=0; i<160; i++){
    for(int j=0; j<8000; j++){
      ave[i]+=(double)count[i][j]/last;
    }
    ofs << i << " ";
    for(int j=0; j<8000-8; j++){
      if(i<64 && entry_type((i-i%16)/16,i%16)!=1){
	if(count[i][j]<ave[i]/100 && count[i][j+1]>ave[i]/100 && count[i][j+2]>ave[i]/100 && count[i][j+3]>ave[i]/100 && count[i][j+4]>ave[i]/100 && count[i][j+5]>ave[i]/100){
	  ofs << (double)j/2 << " ";
	}
      }
      if(i<64 && entry_type((i-i%16)/16,i%16)==1){
	if(count[i][j]<ave[i]/1.5 && count[i][j+1]>ave[i]/1.5 && count[i][j+2]>ave[i]/1.5 && count[i][j+3]>ave[i]/1.5 && count[i][j+4]>ave[i]/1.5 && count[i][j+5]>ave[i]/1.5){
	  ofs << (double)j/2 << " ";
	}
      }
      if(i>63){
	if(count[i][j]<10 && count[i][j+1]<5 && count[i][j+2]<5 && count[i][j+3]>5 && count[i][j+4]>5 && count[i][j+5]>5 && count[i][j+6]>5 && count[i][j+7]>5){
	  ofs << (double)j/2 << " ";
	}
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

