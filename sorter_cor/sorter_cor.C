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


int sorter_cor(int run){

  ifstream ifs(Form("../time_cor/log/log%d.txt",run));  
  int num[200]={};
  double p0[200]={};
  double p1[200]={};
  int tmp[200]={};
  double time_diff=pow(2,31)/1e9;

  for(int i=0; i<200; i++){
    ifs >> num[i] >> p0[i] >> p1[i];
    for(int j=0; j<50; j++){
      if(abs(p1[i] - time_diff * j)<0.2) tmp[i]=j;
    }
    //    cout << num[i] << " " << p1[i] << " "<< tmp[i] << endl;  
  }

  int array[200]={};
  for(int i=0; i<200; i++){
    if(i>0) array[i]=array[i-1];
    for(int j=0; j<200; j++){
      if(num[j]==i) array[i]=tmp[j];
    }
    cout << array[i] << endl;
  }


  
  TFile *fin =new TFile(Form("../sorter/rootfile/run%d_-1_ssgant1.root",run));
  TTree *tree = (TTree*)fin->Get("tree");

  int tmp_domain;
  int tmp_adc;
  double tmp_ts;
  float tmp_energy;
  float tmp_amax;
  
  tree->SetBranchAddress("domain", &tmp_domain);
  tree->SetBranchAddress("ADC", &tmp_adc);
  tree->SetBranchAddress("FineTS", &tmp_ts);
  tree->SetBranchAddress("Energy", &tmp_energy);
  tree->SetBranchAddress("Amax", &tmp_amax);
  
  
  TFile *fout =new TFile(Form("rootfile/run%d_0_ssgant1.root",run),"recreate");
  TTree *tout = new TTree("tree","tree");  
  
  int domain;
  int ADC;
  double FineTS;
  float Energy;
  float Amax;
  
  tout->Branch("domain",&domain,"domain/I");
  tout->Branch("FineTS",&FineTS,"FineTS/D");
  tout->Branch("ADC",&ADC,"ADC/I");
  tout->Branch("Energy",&Energy,"Energy/F");
  tout->Branch("Amax",&Amax,"Amax/F");


  ULong64_t N=tree->GetEntries();
  cout << "Total entry: " << N << endl;
  for(ULong64_t evtn=0; evtn<N; evtn++){

    int loop = (int)((evtn-evtn%(N/200))/(N/200));
      
    domain=-1;
    ADC=-1;
    FineTS=-1;
    Energy=-1;
    Amax=-1;

    if(evtn%10000==0) cout << "\rAnalyzed entry:" << evtn << " loop:" << loop; std::cout << flush;

    tree->GetEntry(evtn);

    domain = tmp_domain;
    ADC = tmp_adc;
    FineTS = tmp_ts;
    if(domain>=64) FineTS = tmp_ts - array[loop]*pow(2,31)*1e3;
    Energy = tmp_energy;
    Amax = tmp_amax;

    tout->Fill();
  }
  
  tout->AutoSave();
  fout->Close();
  
  return 0;
}

