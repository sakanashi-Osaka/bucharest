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


int scaler(int run){
  
  ofstream ofs(Form("log/log%d.txt",run),std::ios::out);
 
  TFile *fin =new TFile(Form("../../sorter_cor/rootfile/run%d_-1.root",run));
  TTree *tree = (TTree*)fin->Get("tree");
  
  int tmp_count[160];
  double ref_ts;
  for(int i=0; i<160; i++) tmp_count[i]=0;

  int domain;
  double FineTS;
  tree->SetBranchAddress("domain",&domain);
  tree->SetBranchAddress("FineTS",&FineTS);

  ULong64_t N=tree->GetEntries();
  cout << "Total entry: " << N << endl;
  for(ULong64_t evtn=0; evtn<N; evtn++){
    
    
    if(evtn%100000==0) cout << "\rAnalyzed entry:" << evtn; std::cout << flush;
    
    tree->GetEntry(evtn);

    tmp_count[domain]++;
    ref_ts = FineTS;
  }
  for(int i=0; i<160; i++) ofs << tmp_count[i] << std::endl;
  ofs << ref_ts << std::endl;
  
  return 0;
}


