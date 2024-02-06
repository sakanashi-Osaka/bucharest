#include <iostream>
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

using namespace std;

const float s1_r1 = 48/2.0;  // inner radius                                                  
const float s1_r2 = 96/2.0;  // outer radiu                                                   
const float s1_dist = 40;
//const float s1_dist = 35;

int entry_type(int board, int ch);
int get_front_ch_order(Int_t mod, Int_t ch);
int get_rear_ch_order(Int_t mod, Int_t ch);
float get_rear_phi(Int_t ch);
float get_front_r(int front_ch);
float get_front_theta(int front_ch);
float get_front_ex(float theta, float ene);

int ene_culc(int run){

  TFile *fin =new TFile(Form("../eventbuild/event%d.root",run));
  //  TFile *fin =new TFile(Form("../eventbuild/test%d.root",run));
  TTree *tree = (TTree*)fin->Get("event");

  float tmp_energy[N_BOARD][N_CH];
  float tmp_amax[N_BOARD][N_CH];
  float tmp_ex[N_BOARD][N_CH];

  tree->SetBranchAddress("Energy",tmp_energy);
  tree->SetBranchAddress("Amax",tmp_amax);
  tree->SetBranchAddress("front_ex",tmp_ex);
  
  
  int board=-1;
  int ch=-1;

  int hit_front;
  int hit_rear;
  int board_front;
  int board_rear;
  int ch_front;
  int ch_rear;
  double ene_front;
  double ene_rear;

  int fr;
  float rear_phi;
  int front_ch;
  float front_phi, front_r, front_theta;
  float front_ex[N_BOARD][N_CH];
  float cor_ex[N_BOARD][N_CH];

  TCanvas* c0 = new TCanvas( "name0" , "title", 800, 800);
  TCanvas* c1 = new TCanvas( "name1" , "title", 800, 800);
  TCanvas* c2 = new TCanvas( "name2" , "title", 800, 800);
  TCanvas* c3 = new TCanvas( "name3" , "title", 800, 800);
  
  TH1F *hist0 = new TH1F("name0","title",1000,-1,14);
  TH1F *hist1 = new TH1F("name1","title",1000,-1,14);
  TH1F *hist2 = new TH1F("name2","title",1000,-1,14);
  TH1F *hist3 = new TH1F("name3","title",1000,-1,14);

  TH1F *hist4 = new TH1F("name4","title",1000,-1,14);
  TH1F *hist5 = new TH1F("name5","title",1000,-1,14);
  TH1F *hist6 = new TH1F("name6","title",1000,-1,14);
  TH1F *hist7 = new TH1F("name7","title",1000,-1,14);

  TH1F *hist8 = new TH1F("name8","title",1000,-1,14);
  TH1F *hist9 = new TH1F("name9","title",1000,-1,14);
  TH1F *hist10 = new TH1F("name10","title",1000,-1,14);
  TH1F *hist11 = new TH1F("name11","title",1000,-1,14);

  TH1F *hist12 = new TH1F("name12","title",1000,-1,14);
  TH1F *hist13 = new TH1F("name13","title",1000,-1,14);
  TH1F *hist14 = new TH1F("name14","title",1000,-1,14);
  TH1F *hist15 = new TH1F("name15","title",1000,-1,14);

  ULong64_t N=tree->GetEntries();
  cout << "Total entry: " << N << endl;
  for(ULong64_t evtn=0; evtn<N; evtn++){
    
    if(evtn%10000==0) cout << "\rAnalyzed entry:" << evtn; std::cout << flush;
    
    tree->GetEntry(evtn);

    hit_front = 0; hit_rear = 0;
    
    for(int i=0; i<N_BOARD; i++){
      for(int j=0; j<N_CH; j++){

	if(i*16+j>65 && i*16+j<138 && j!=0 && j!=1){
	  if(tmp_energy[i][j]>0.5){
	    hit_front++; 
	    board_front=i;
	    ch_front=j;
	    ene_front=tmp_energy[i][j];				   
	  }
	}
	if(i*16+j>137 && i*16+j<156 && j!=0 && j!=1){
	  if(tmp_energy[i][j]>0.5){
	    hit_rear++;
	    board_rear=i;
	    ch_rear=j;
	    ene_rear=tmp_energy[i][j];				   
	  }
	}
	
      }
    }
   
    if(hit_front==1 && hit_rear==1 && abs(ene_front-ene_rear)<1){
      //      cout << board_front << " " << ch_front << " " <<  board_rear << " " << ch_rear << endl;

      int front_ch = get_front_ch_order(board_front, ch_front);
      float front_theta = get_front_theta(front_ch);
      int rear_ch = get_rear_ch_order(board_rear, ch_rear);
      float rear_phi = get_rear_phi(rear_ch);
      //      front_ex[board_front][ch_front] = get_front_ex(front_theta, ene_rear);
      front_ex[board_front][ch_front] = get_front_ex(front_theta, ene_front);

      float beam_x=0;
      //      float beam_x=2;
      float beam_y=2.0;
      //      float beam_y=-2;
      
      float tmp_x, tmp_y, tmp_z, tmp_theta;
      tmp_x = s1_dist * tan(front_theta)*cos(rear_phi) - beam_x;
      tmp_y = s1_dist * tan(front_theta)*sin(rear_phi) - beam_y;
      tmp_z = s1_dist;
      tmp_theta = acos(tmp_z/pow(pow(tmp_x,2)+pow(tmp_y,2)+pow(tmp_z,2),0.5));
      cor_ex[board_front][ch_front] = get_front_ex(tmp_theta, ene_front);

      if(board_front==4 && ch_front==2 && rear_ch==1){
	//      if(board_front==4 && ch_front==10 && rear_ch==0){
	hist0->Fill(front_ex[board_front][ch_front]);
	hist1->Fill(cor_ex[board_front][ch_front]);
      }
      if(board_front==4 && ch_front==2 && rear_ch==2){
      //      if(board_front==4 && ch_front==10 && rear_ch==3){
	hist2->Fill(front_ex[board_front][ch_front]);
	hist3->Fill(cor_ex[board_front][ch_front]);
      }
      if(board_front==5 && ch_front==6 && rear_ch==5){	
	//      if(board_front==5 && ch_front==12 && rear_ch==4){
	hist4->Fill(front_ex[board_front][ch_front]);
	hist5->Fill(cor_ex[board_front][ch_front]);
      }
      if(board_front==5 && ch_front==6 && rear_ch==6){
	//      if(board_front==5 && ch_front==12 && rear_ch==7){
	hist6->Fill(front_ex[board_front][ch_front]);
	hist7->Fill(cor_ex[board_front][ch_front]);
      }
      if(board_front==6 && ch_front==8 && rear_ch==9){
	//      if(board_front==6 && ch_front==14 && rear_ch==8){
	hist8->Fill(front_ex[board_front][ch_front]);
	hist9->Fill(cor_ex[board_front][ch_front]);
      }
      if(board_front==6 && ch_front==8 && rear_ch==10){
	//      if(board_front==6 && ch_front==14 && rear_ch==11){
	hist10->Fill(front_ex[board_front][ch_front]);
	hist11->Fill(cor_ex[board_front][ch_front]);
      }
      if(board_front==8 && ch_front==2 && rear_ch==13){
	//      if(board_front==8 && ch_front==2 && rear_ch==12){
	hist12->Fill(front_ex[board_front][ch_front]);
	hist13->Fill(cor_ex[board_front][ch_front]);
      }
      if(board_front==8 && ch_front==2 && rear_ch==14){
	//      if(board_front==8 && ch_front==2 && rear_ch==15){
	hist14->Fill(front_ex[board_front][ch_front]);
	hist15->Fill(cor_ex[board_front][ch_front]);
      }

    }
  }

  hist0->SetLineColor(1);
  hist1->SetLineColor(2);
  hist2->SetLineColor(3);
  hist3->SetLineColor(4);

  hist4->SetLineColor(1);
  hist5->SetLineColor(2);
  hist6->SetLineColor(3);
  hist7->SetLineColor(4);

  hist8->SetLineColor(1);
  hist9->SetLineColor(2);
  hist10->SetLineColor(3);
  hist11->SetLineColor(4);

  hist12->SetLineColor(1);
  hist13->SetLineColor(2);
  hist14->SetLineColor(3);
  hist15->SetLineColor(4);

  c0->cd();
  hist0->Draw();
  hist1->Draw("same");
  hist2->Draw("same");
  hist3->Draw("same");

  c1->cd();
  hist4->Draw();
  hist5->Draw("same");
  hist6->Draw("same");
  hist7->Draw("same");

  c2->cd();
  hist8->Draw();
  hist9->Draw("same");
  hist10->Draw("same");
  hist11->Draw("same");

  c3->cd();
  hist12->Draw();
  hist13->Draw("same");
  hist14->Draw("same");
  hist15->Draw("same");

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

float get_rear_phi(Int_t ch){
  float phi=0;
  if(ch%2==0 && ch>=0) phi = 90 + 22.5/2.0 + 22.5*(ch+1);
  if(ch%2==1 && ch>=0) phi = 90 +  22.5/2.0 + 22.5*(ch-1);  
  return phi*TMath::DegToRad();
  //  return phi;
}

int get_front_ch_order(Int_t mod, Int_t ch){
  if(mod==4 && ch>=2) return ch-2;
  if(mod==5 && ch>=2) return ch+12;
  if(mod==6 && ch>=2) return ch+26;    
  if(mod==7 && ch>=2) return ch+40;
  if(mod==8 && ch>=2 && ch<=9) return ch+54;      
  return 100;
}

int get_rear_ch_order(Int_t mod, Int_t ch){
  if(mod==8 && ch>=10 && ch<14) return ch+2;      
  if(mod==8 && ch>=14 && ch<16) return ch-14;      
  if(mod==9 && ch>=2 && ch<8) return ch;      
  if(mod==9 && ch>=8 && ch<12) return 19-ch;      
  return 100;
}

float get_front_phi(int front_ch){
  float phi=0;
  if(front_ch>=0 && front_ch<64){
    phi = 135 + 90*((int)(front_ch/16));
  }
  return phi*TMath::DegToRad();
}

float get_front_r(int front_ch){
  float strp_wid = (s1_r2 - s1_r1)/16.0;
  int section = (int)(front_ch/16);
  float r=0;
  
  if(section==0 || section==3){
    if(front_ch%2==0) r = s1_r1 + strp_wid*(1.5) + strp_wid*(front_ch%16);
    if(front_ch%2==1) r = s1_r1 + strp_wid*(0.5) + strp_wid*((front_ch%16) - 1);    
  }

  if(section==1 || section==2){
    if(front_ch%2==0) r = s1_r2 - strp_wid*(1.5) - strp_wid*(front_ch%16);
    if(front_ch%2==1) r = s1_r2 - strp_wid*(0.5) - strp_wid*((front_ch%16) - 1);    
  }
  
  return r;
}

float get_front_theta(int front_ch){
  
  float strp_wid = (s1_r2 - s1_r1)/16.0;
  int section = (int)(front_ch/16);
  float r=0;
  
  if(section==0 || section==3){
    if(front_ch%2==0) r = s1_r1 + strp_wid*(1.5) + strp_wid*(front_ch%16);
    if(front_ch%2==1) r = s1_r1 + strp_wid*(0.5) + strp_wid*((front_ch%16) - 1);
  }
  
  if(section==1 || section==2){
    if(front_ch%2==0) r = s1_r2 - strp_wid*(1.5) - strp_wid*(front_ch%16);
    if(front_ch%2==1) r = s1_r2 - strp_wid*(0.5) - strp_wid*((front_ch%16) - 1);
  }
  
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
