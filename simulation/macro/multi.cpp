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

using namespace std; 

#define AMU 931.494
#define PI 3.14159265
#define N_BOARD 10
#define N_CH 16
#define N_HIT_MAX 10

double culc_Amax(double ene, int fch, int rch);

const double s1_r1 = 48/2.0;  // inner radius      
const double s1_r2 = 96/2.0;  // outer radius
double s1_dist = 40;

int fch[16][16];
int rch[16][16];
double pm0[16][16];
double pm1[16][16];
double pm2[16][16];
double pm3[16][16];
double ps0[16][16];
double ps1[16][16];
double ps2[16][16];
double ps3[16][16];

double Ex4He;
double K4He, th4He, ph4He;
double K12C, th12C, ph12C;
double K3a[3], th3a[3], ph3a[3];
int chf4He, chr4He, chf12C, chr12C, chf3a[3], chr3a[3];
double Amax12C, Amax3a[3];

int main(int argc, char *argv[]){
  
  TFile *fin =new TFile(Form("../rootfile/%s.root",argv[1]));
  TTree *tree = (TTree*)fin->Get("tree");
  
  
  tree->SetBranchAddress("mex4He",&Ex4He);
  tree->SetBranchAddress("mK4He",&K4He);
  tree->SetBranchAddress("mth4He",&th4He);
  tree->SetBranchAddress("mph4He",&ph4He);
  tree->SetBranchAddress("m4HechF",&chf4He);
  tree->SetBranchAddress("m4HechR",&chr4He);
  
  tree->SetBranchAddress("mK12C",&K12C);
  tree->SetBranchAddress("mth12C",&th12C);
  tree->SetBranchAddress("mph12C",&ph12C);
  tree->SetBranchAddress("m12CchF",&chf12C);
  tree->SetBranchAddress("m12CchR",&chr12C);
  tree->SetBranchAddress("mAmax12C",&Amax12C);

  tree->SetBranchAddress("mKf3a",K3a);
  tree->SetBranchAddress("mth3a",th3a);
  tree->SetBranchAddress("mph3a",ph3a);
  tree->SetBranchAddress("m3achF",chf3a);
  tree->SetBranchAddress("m3achR",chr3a);
  tree->SetBranchAddress("mAmax3a",Amax3a);
  
  
  // Amax generator                                                                                       
  ifstream ifs("/home/sakra/exp/Bucharest2022/simulation/Amax_gen/prm/fit_Amax_mean2305.txt");
  for(int i=0; i<16; i++){
    for(int j=0; j<16; j++){
      ifs >> fch[i][j] >> rch[i][j] >> pm0[i][j] >> pm1[i][j] >> pm2[i][j] >> pm3[i][j];
    }
  }
  ifstream ifs1("/home/sakra/exp/Bucharest2022/simulation/Amax_gen/prm/fit_Amax_sigma2305.txt");
  for(int i=0; i<16; i++){
    for(int j=0; j<16; j++){
      ifs1 >> fch[i][j] >> rch[i][j] >> ps0[i][j] >> ps1[i][j] >> ps2[i][j] >> ps3[i][j];
    }
  }



  TFile *fout =new TFile(Form("../rootfile/multi_%s.root",argv[1]),"recreate");
  TTree *multi = new TTree("multi","multi"); 

  double Ka, tha, pha;
  int chfa, chra;
  int C_cth, C_cph;
  int C_ath, C_aph;
  int flag_2a;
  double Amaxa;
  

  multi->Branch("Ex4He",&Ex4He,"Ex4He/D");
  multi->Branch("K4He",&K4He,"K4He/D");
  multi->Branch("th4He",&th4He,"th4He/D");
  multi->Branch("ph4He",&ph4He,"ph4He/D");
  multi->Branch("chf4He",&chf4He,"chf4He/I");
  multi->Branch("chr4He",&chr4He,"chr4He/I");

  multi->Branch("K12C",&K12C,"K12C/D");
  multi->Branch("th12C",&th12C,"th12C/D");
  multi->Branch("ph12C",&ph12C,"ph12C/D");
  multi->Branch("chf12C",&chf12C,"chf12C/I");
  multi->Branch("chr12C",&chr12C,"chr12C/I");
  multi->Branch("Amax12C",&Amax12C,"Amax12C/D");

  multi->Branch("Ka",&Ka,"Ka/D");
  multi->Branch("tha",&tha,"tha/D");
  multi->Branch("pha",&pha,"pha/D");
  multi->Branch("chfa",&chfa,"chfa/I");
  multi->Branch("Amaxa",&Amaxa,"Amaxa/D");

  multi->Branch("C_cth",&C_cth,"C_cth/I");
  multi->Branch("C_cph",&C_cph,"C_cph/I");
  multi->Branch("C_ath",&C_ath,"C_ath/I");
  multi->Branch("C_aph",&C_aph,"C_aph/I");
  multi->Branch("flag_2a",&flag_2a,"flag_2a/I");
  
  ULong64_t N=tree->GetEntries();
  cout << "Total entry: " << N << endl;
  for(ULong64_t evtn=0; evtn<N; evtn++){
    if(evtn%100000==0) cout << evtn << endl;
    tree->GetEntry(evtn);

    C_cth=0; C_cph=0;
    C_ath=0; C_aph=0;
    flag_2a=0;
    Ka=-1000; tha=-1000; pha=-1000; chfa=-1; chra=-1;
    
    //C_flag
    if(th12C>0 && th4He>0 && th4He>-4*th12C+200) C_cth=1;
    if(abs(ph4He-ph12C)>155 && abs(ph4He-ph12C)<205) C_cph=1;
    //	if(K4He+K12C>15.6 && K4He+K12C<16.7) C_E=1; 



    // 12C
    Amax12C = Amax12C - culc_Amax(K12C,chf12C,chr12C);

    
    // multihit alpha
    for(int i=0; i<3; i++){
      if(chf3a[i%3]==chf3a[(i+1)%3] && chr3a[i]==chr3a[(i+1)%3] && chf3a[i]>-1 && chr3a[i]>-1){ //a0 && a1
	Ka = K3a[i] + K3a[(i+1)%3];
	tha = th3a[i];
	pha = ph3a[i];
	chfa = chf3a[i];
	chra = chr3a[i];
	Amaxa = (Amax3a[i]*K3a[i]+Amax3a[(i+1)%3]*K3a[(i+1)%3])/(K3a[i]+K3a[(i+1)%3]) - culc_Amax(Ka,chfa,chra);
	flag_2a = 1;
	if(tha>0 && th4He>0 && th4He>-4*tha+200) C_ath=1;
	if(abs(ph4He-pha)>155 && abs(ph4He-pha)<205) C_aph=1;
      }
    }
    multi->Fill();
  }
  
  multi->AutoSave();
  fout->Close();
  
  return 0;
}



/******************************************************/

double culc_Amax(double ene, int fch, int rch){
  double tmp_Amax;
  
  if(fch>-1 && rch>-1){
    tmp_Amax = pm0[fch][rch] + pm1[fch][rch]*ene + pm2[fch][rch]/(ene-pm3[fch][rch]);
  }
  return tmp_Amax;
}
