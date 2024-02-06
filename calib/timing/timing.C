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

#define ch 160
#define peak_n 1

int timing(int run){

  TSpectrum *s[peak_n];
  TH1F *h[ch];
  TF1 *f[ch][peak_n];
  double *xpeaks[ch];
  double peak[peak_n];
  double mean[ch][peak_n];
  double sigma[ch][peak_n];
  int n;

  ofstream ofs(Form("time%d.txt",run),std::ios::app);

  TFile *fin =new TFile(Form("../../eventbuild_cor/rootfile/run%d.root",run));
  TTree *event = (TTree*)fin->Get("event");

  for(Int_t i=0;i<ch;i++){
    int Board = (int)((i-i%16)/16);
    int Ch = i%16;

    cout << "board:" << Board << " ch:" << Ch << "   ";

    s[i] = new TSpectrum(peak_n);
    if(i<64) h[i] = new TH1F(Form("h_%d",i),Form("h_%d",i),500,-600e3,-100e3);
    if(i>63) h[i] = new TH1F(Form("h_%d",i),Form("h_%d",i),500,400e3,900e3);
    event->Draw(Form("ts_diff[%d][%d]>>h_%d",Board,Ch,i),Form("Energy[%d][%d]>1",Board,Ch),"",5e6);
    s[i]->Search(h[i],1,"new",0.1);
    xpeaks[i]=s[i]->GetPositionX();

    for(Int_t j=0;j<peak_n;j++){
      peak[j]=xpeaks[i][j];
    }
    n = sizeof(peak)/sizeof(peak[0]);
    std::sort(peak,peak+n);

    cout << i << " ";
    for(Int_t j=0;j<peak_n;j++){
      cout << peak[j] << " ";
    }
    cout << endl;

    for(Int_t j=0;j<peak_n;j++){
      f[i][j]=new TF1(Form("f%d_%d",i,j),"gaus",peak[j]-5e3,peak[j]+5e3);
      f[i][j]->SetParameters(100,peak[j],5);
    }

    h[i]->Draw();
    for(Int_t j=0;j<peak_n;j++){
      h[i]->Fit(Form("f%d_%d",i,j),"+","",peak[j]-5e3,peak[j]+5e3);
      mean[i][j]=f[i][j]->GetParameter(1);
      sigma[i][j]=f[i][j]->GetParameter(2);
      cout << mean[i][j] << " " << sigma[i][j] << endl;
      ofs << i << " " << j << " " << mean[i][j] << " " << sigma[i][j] << endl;
    }

  }

  ofs.close();

  return 0;
}
