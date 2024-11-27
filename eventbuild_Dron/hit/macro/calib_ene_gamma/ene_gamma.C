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

#define ch 50
#define peak_n 1

int ene_gamma(int run){

  TSpectrum *s[peak_n];
  TH1F *h[ch];
  TF1 *f[ch][peak_n];
  double *xpeaks[ch];
  double peak[peak_n];
  double mean[ch][peak_n];
  double sigma[ch][peak_n];
  int n;

  ofstream ofs(Form("log/ene%d.txt",run),std::ios::app);

  TFile *fin =new TFile(Form("../rootfile/single%d.root",run));
  TTree *tree = (TTree*)fin->Get("tree");

  for(Int_t i=0;i<ch;i++){
  //  for(Int_t i=0;i<2;i++){

    s[i] = new TSpectrum(peak_n);
    h[i] = new TH1F(Form("h_%d",i),Form("h_%d",i),100,4.1,5.1);
    tree->Draw(Form("Kg_cor>>h_%d",i),Form("chg[0]==%d",i),"");
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
      f[i][j]=new TF1(Form("f%d_%d",i,j),"gaus",peak[j]-0.08,peak[j]+0.08);
      f[i][j]->SetParameters(150,peak[j],0.05);
    }

    h[i]->Draw();
    for(Int_t j=0;j<peak_n;j++){
      h[i]->Fit(Form("f%d_%d",i,j),"+","",peak[j]-0.08,peak[j]+0.08);
      mean[i][j]=f[i][j]->GetParameter(1);
      sigma[i][j]=f[i][j]->GetParameter(2);
      cout << mean[i][j] << " " << sigma[i][j] << endl;
      ofs << i << " " << j << " " << mean[i][j] << " " << sigma[i][j] << endl;
    }

  }

  ofs.close();

  return 0;
}
