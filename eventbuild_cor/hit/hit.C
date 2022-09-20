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

#define PI 3.14159265
#define N_BOARD 10
#define N_CH 16
#define N_HIT_MAX 10

const float s1_r1 = 48/2.0;  // inner radius      
const float s1_r2 = 96/2.0;  // outer radius
const float s1_dist = 40;

int entry_type(int board, int ch);
int get_ch(int domain);
float get_front_theta(int front_ch);
std::pair<float,float> get_angle(float x, float y, int ch_f, int ch_r);
float get_front_ex(float theta, float ene);
void fill_data(int seg, int num_fi, int num_ff, int num_ri, int num_rf);
double detect3a(double *enea, int *chfa, int *chra);

//branch
int hit_n[4]={0,0,0,0};
int ch_f[4][N_HIT_MAX]={};
int ch_r[4][N_HIT_MAX]={};
int ch_g[N_HIT_MAX]={};
int ch_b[N_HIT_MAX]={};
float Energy_f[4][N_HIT_MAX]={};
float Energy_r[4][N_HIT_MAX]={};
float front_ex[4][N_HIT_MAX]={};
float cor_ex[4][N_HIT_MAX]={};
double ts_diff_f[4][N_HIT_MAX]={};
double ts_diff_r[4][N_HIT_MAX]={};
double ts_diff_g[N_HIT_MAX]={};
double ts_diff_b[N_HIT_MAX]={};
float Amax[4][N_HIT_MAX]={};
float Gamma[N_HIT_MAX]={};
float BGO[N_HIT_MAX]={};
int run_n=-1;
double mex=-1;

//analysis
int count_f[4]={0,0,0,0};
int count_r[4]={0,0,0,0};
int count_g=0;
int count_b=0;
float ene_ex[4][N_HIT_MAX]={};
float ene_f[4][N_HIT_MAX]={};
float ene_r[4][N_HIT_MAX]={};
float Amax_f[4][N_HIT_MAX]={};
float Amax_r[4][N_HIT_MAX]={};
double ts_f[4][N_HIT_MAX]={};
double ts_r[4][N_HIT_MAX]={};
int domain_f[4][N_HIT_MAX]={};
int domain_r[4][N_HIT_MAX]={};
int type_f[4][N_HIT_MAX]={};
int type_r[4][N_HIT_MAX]={};

double enea[3]={};
int chfa[3]={};
int chra[3]={};
double tha[3]={};
double pha[3]={};


int hit(int run){

  ifstream ifs("../timing.prm");
  int domain[10][16];
  double p0[10][16];
  double p1[10][16];
  double p2[10][16];
  for(int i=0; i<10; i++){
    for(int j=0; j<16; j++){
    ifs >> domain[i][j] >> p0[i][j]>> p1[i][j] >> p2[i][j];
    }
  }
  
  //  TFile *fin =new TFile(Form("../rootfile/event%d.root",run));
  TFile *fin =new TFile(Form("../rootfile/test%d.root",run));
  TTree *event = (TTree*)fin->Get("event");

  float tmp_energy[N_BOARD][N_CH]={};
  float tmp_amax[N_BOARD][N_CH]={};
  float tmp_front_ex[N_BOARD][N_CH]={};
  double tmp_ts_diff[N_BOARD][N_CH]={};
  
  event->SetBranchAddress("Energy",tmp_energy);
  event->SetBranchAddress("Amax_cor",tmp_amax);
  event->SetBranchAddress("front_ex",tmp_front_ex);
  event->SetBranchAddress("ts_diff",tmp_ts_diff);
  
  //  TFile *fout =new TFile(Form("rootfile/hit%d.root",run),"recreate");
  TFile *fout =new TFile(Form("rootfile/test%d.root",run),"recreate");
  TTree *hit = new TTree("hit","hit"); 

  
  hit->Branch("hit_n",hit_n,"hit_n[4]/I");
  hit->Branch("ch_f",ch_f,"ch_f[4][10]/I");
  hit->Branch("ch_r",ch_r,"ch_r[4][10]/I");
  hit->Branch("ch_g",ch_g,"ch_g[10]/I");
  hit->Branch("ch_b",ch_b,"ch_b[10]/I");
  hit->Branch("Energy_f",Energy_f,"Energy_f[4][10]/F");
  hit->Branch("Energy_r",Energy_r,"Energy_r[4][10]/F");
  hit->Branch("front_ex",front_ex,"front_ex[4][10]/F");
  hit->Branch("cor_ex",cor_ex,"cor_ex[4][10]/F");
  hit->Branch("ts_diff_f",ts_diff_f,"ts_diff_f[4][10]/D");
  hit->Branch("ts_diff_r",ts_diff_r,"ts_diff_r[4][10]/D");
  hit->Branch("ts_diff_g",ts_diff_g,"ts_diff_g[10]/D");
  hit->Branch("ts_diff_b",ts_diff_b,"ts_diff_b[10]/D");
  hit->Branch("Amax",Amax,"Amax[4][10]/F");
  hit->Branch("Gamma",Gamma,"Gamma[10]/F");
  hit->Branch("BGO",BGO,"BGO[10]/F");
  hit->Branch("run_n",&run_n,"run_n/I");
  hit->Branch("mex",&mex,"mex/D");

  
  
  ULong64_t N=event->GetEntries();
  cout << "Total entry: " << N << endl;
  for(ULong64_t evtn=0; evtn<N; evtn++){
    //  for(ULong64_t evtn=0; evtn<1000000; evtn++){

    if(evtn%10000==0) cout << "\rAnalyzed entry:" << evtn; std::cout << flush;

    event->GetEntry(evtn); 

    for(int i=0; i<4; i++){
      for(int j=0; j<N_HIT_MAX; j++){
	ch_f[i][j]=-1;
	ch_r[i][j]=-1;
	ch_g[j]=-1;
	ch_b[j]=-1;
	Energy_f[i][j]=-1;
	Energy_r[i][j]=-1;
	Amax[i][j]=-1;
	front_ex[i][j]=-1;
	cor_ex[i][j]=-1; 
	ts_diff_f[i][j]=-1e6;
	ts_diff_r[i][j]=-1e6;
	ts_diff_g[j]=-1e6;
	ts_diff_b[j]=-1e6;
	Gamma[i]=-1;
	BGO[i]=-1;
	run_n=-1;
	hit_n[i]=0;
	mex=-1;
	
	count_f[i]=0;
	count_r[i]=0;
	count_g=0;
	count_b=0;
	ene_ex[i][j]=0;
	ene_f[i][j]=0;
	ene_r[i][j]=0;
	Amax_f[i][j]=0;
	Amax_r[i][j]=0;
	ts_f[i][j]=0;
	ts_r[i][j]=0;
	domain_f[i][j]=0;
	domain_r[i][j]=0;
	type_f[i][j]=0;
	type_r[i][j]=0;
      }
    }


    run_n=run;
    
    for(int i=0; i<N_BOARD; i++){
      for(int j=0; j<N_CH; j++){
	int type = entry_type(i,j);

	if(tmp_energy[i][j]>0.3 && tmp_energy[i][j]<30){
	  
	  if(abs(type-22)<3){
	    ene_f[type-20][count_f[type-20]]=tmp_energy[i][j];
	    ene_ex[type-20][count_f[type-20]]=tmp_front_ex[i][j];
	    ts_f[type-20][count_f[type-20]]=tmp_ts_diff[i][j]-p1[i][j];
	    Amax_f[type-20][count_f[type-20]]=tmp_amax[i][j];
	    domain_f[type-20][count_f[type-20]]=i*16+j;
	    type_f[type-20][count_f[type-20]]=type;
	    count_f[type-20]++;
	  }
	  if(abs(type-32)<3){
	    ene_r[type-30][count_r[type-30]]=tmp_energy[i][j];
	    ts_r[type-30][count_r[type-30]]=tmp_ts_diff[i][j]-p1[i][j];
	    Amax_r[type-30][count_r[type-30]]=tmp_amax[i][j];
	    domain_r[type-30][count_r[type-30]]=i*16+j;
	    type_r[type-30][count_r[type-30]]=type;
	    count_r[type-30]++;
	  }

	  if(type==10){
	    Gamma[count_g]=tmp_energy[i][j];
	    ts_diff_g[count_g]=tmp_ts_diff[i][j]-p1[i][j];
	    ch_g[count_g]=i*16+j;
	    count_g++;
	  }
	  if(type==11){
	    BGO[count_b]=tmp_energy[i][j];
	    ts_diff_b[count_b]=tmp_ts_diff[i][j]-p1[i][j];
	    ch_b[count_b]=i*16+j;
	    count_b++;
	  }
	  
	}
      }
    }

    for(int seg=0; seg<4; seg++){
      if(count_f[seg]==0 || count_r[seg]==0) continue;

      if(count_f[seg]==1 && count_r[seg]==1){ //1vs1 event
	if(abs(ene_f[seg][0]-ene_r[seg][0])<0.5){;
	  fill_data(seg,0,0,0,0);	  
	  hit_n[seg]=11;
	}
      }
      if(count_f[seg]==2 && count_r[seg]==1){ //2vs1 event
	if(abs(ene_f[seg][0]+ene_f[seg][1]-ene_r[seg][0])<1){;
	  fill_data(seg,0,0,0,0);	  
	  fill_data(seg,1,1,0,1);	  
	  hit_n[seg]=21;
	}
      }
      if(count_f[seg]==2 && count_r[seg]==2){ //2vs2 event
	if(abs(ene_f[seg][0]-ene_r[seg][0])<0.5 && abs(ene_f[seg][1]-ene_r[seg][1])<0.5){;
	  fill_data(seg,0,0,0,0);	  
	  fill_data(seg,1,1,1,1);	  
	  hit_n[seg]=22;
	}else if(abs(ene_f[seg][0]-ene_r[seg][1])<0.5 && abs(ene_f[seg][1]-ene_r[seg][0])<0.5){;
	  fill_data(seg,0,0,1,0);	  
	  fill_data(seg,1,1,0,1);	  
	  hit_n[seg]=22;
	}
      }
      if(count_f[seg]==3 && count_r[seg]==1){ //3vs1 event
	if(abs(ene_f[seg][0]+ene_f[seg][1]+ene_f[seg][2]-ene_r[seg][0])<1){;
	  fill_data(seg,0,0,0,0);	  
	  fill_data(seg,1,1,0,1);	  
	  fill_data(seg,2,2,0,2);	  
	  hit_n[seg]=31;
	}
      }
      if(count_f[seg]==3 && count_r[seg]==2){ //3vs2 event
	if(abs(ene_f[seg][0]+ene_f[seg][1]+ene_f[seg][2]-ene_r[seg][0]-ene_r[seg][1])<1){
	  if(abs(ene_f[seg][0]+ene_f[seg][1]-ene_r[seg][0])<0.5 && abs(ene_f[seg][2]-ene_r[seg][1])<0.5){
	    fill_data(seg,0,0,0,0);	  
	    fill_data(seg,1,1,0,1);	  
	    fill_data(seg,2,2,1,2);	  
	    hit_n[seg]=32;
	  }else if(abs(ene_f[seg][0]+ene_f[seg][1]-ene_r[seg][1])<0.5 && abs(ene_f[seg][2]-ene_r[seg][0])<0.5){
	    fill_data(seg,0,0,1,0);	  
	    fill_data(seg,1,1,1,1);	  
	    fill_data(seg,2,2,0,2);	  
	    hit_n[seg]=32;
	  }else if(abs(ene_f[seg][0]+ene_f[seg][2]-ene_r[seg][0])<0.5 && abs(ene_f[seg][1]-ene_r[seg][1])<0.5){
	    fill_data(seg,0,0,0,0);	  
	    fill_data(seg,2,2,0,2);	  
	    fill_data(seg,1,1,1,1);	  
	    hit_n[seg]=32;
	  }else if(abs(ene_f[seg][0]+ene_f[seg][2]-ene_r[seg][1])<0.5 && abs(ene_f[seg][1]-ene_r[seg][0])<0.5){
	    fill_data(seg,0,0,1,0);	  
	    fill_data(seg,2,2,1,2);	  
	    fill_data(seg,1,1,0,1);	  
	    hit_n[seg]=32;
	  }else if(abs(ene_f[seg][1]+ene_f[seg][2]-ene_r[seg][0])<0.5 && abs(ene_f[seg][0]-ene_r[seg][1])<0.5){
	    fill_data(seg,1,1,0,1);	  
	    fill_data(seg,2,2,0,2);	  
	    fill_data(seg,0,0,1,0);	  
	    hit_n[seg]=32;
	  }else if(abs(ene_f[seg][1]+ene_f[seg][2]-ene_r[seg][1])<0.5 && abs(ene_f[seg][0]-ene_r[seg][0])<0.5){
	    fill_data(seg,1,1,1,1);	  
	    fill_data(seg,2,2,1,2);	  
	    fill_data(seg,0,0,0,0);	  
	    hit_n[seg]=32;
	  }else{
	    continue;
	  }
	  double tmp_cor_ex =cor_ex[(seg+2)%4][0];
	  for(int tmp=0; tmp<3; tmp++){
	    chfa[tmp]=ch_f[seg][tmp];
	    chra[tmp]=ch_r[seg][tmp];
	    enea[tmp]=Energy_f[seg][tmp];
	  }      
	  double ex12C = detect3a(enea,chfa,chra);
	  cout << ex12C << " " << tmp_cor_ex  << endl;
	  cout << endl;
	  mex=ex12C;
	}
      }
      if(count_f[seg]==3 && count_r[seg]==3){ //3vs3 event
	if(abs(ene_f[seg][0]+ene_f[seg][1]+ene_f[seg][2]-ene_r[seg][0]-ene_r[seg][1]-ene_r[seg][2])<1){
	  if(abs(ene_f[seg][0]-ene_r[seg][0])<0.5 && abs(ene_f[seg][1]-ene_r[seg][1])<0.5 && abs(ene_f[seg][2]-ene_r[seg][2])<0.5){
	    fill_data(seg,0,0,0,0);	  
	    fill_data(seg,1,1,1,1);	  
	    fill_data(seg,2,2,2,2);	  
	    hit_n[seg]=33;
	  }else if(abs(ene_f[seg][0]-ene_r[seg][0])<0.5 && abs(ene_f[seg][1]-ene_r[seg][2])<0.5 && abs(ene_f[seg][2]-ene_r[seg][1])<0.5){
	    fill_data(seg,0,0,0,0);	  
	    fill_data(seg,1,1,2,1);	  
	    fill_data(seg,2,2,1,2);	  
	    hit_n[seg]=33;
	  }else if(abs(ene_f[seg][0]-ene_r[seg][1])<0.5 && abs(ene_f[seg][1]-ene_r[seg][0])<0.5 && abs(ene_f[seg][2]-ene_r[seg][2])<0.5){
	    fill_data(seg,0,0,1,0);	  
	    fill_data(seg,1,1,0,1);	  
	    fill_data(seg,2,2,2,2);	  
	    hit_n[seg]=33;
	  }else if(abs(ene_f[seg][0]-ene_r[seg][1])<0.5 && abs(ene_f[seg][1]-ene_r[seg][2])<0.5 && abs(ene_f[seg][2]-ene_r[seg][0])<0.5){
	    fill_data(seg,0,0,1,0);	  
	    fill_data(seg,1,1,2,1);	  
	    fill_data(seg,2,2,0,2);	  
	    hit_n[seg]=33;
	  }else if(abs(ene_f[seg][0]-ene_r[seg][2])<0.5 && abs(ene_f[seg][1]-ene_r[seg][0])<0.5 && abs(ene_f[seg][2]-ene_r[seg][1])<0.5){
	    fill_data(seg,0,0,1,2);	  
	    fill_data(seg,1,1,2,0);	  
	    fill_data(seg,2,2,0,1);	  
	    hit_n[seg]=33;
	  }else if(abs(ene_f[seg][0]-ene_r[seg][2])<0.5 && abs(ene_f[seg][1]-ene_r[seg][1])<0.5 && abs(ene_f[seg][2]-ene_r[seg][0])<0.5){
	    fill_data(seg,0,0,1,2);	  
	    fill_data(seg,1,1,2,1);	  
	    fill_data(seg,2,2,0,0);
	    hit_n[seg]=33;
	  }
	}
      }
      
    }
    //    cout << hit_n[0] << hit_n[1] << hit_n[2] << hit_n[3] << endl;
    
    hit->Fill();
    
  }
  
  hit->AutoSave();
  fout->Close();
  
  return 0;
}


int get_ch(int domain){

  //f_seg0
  if(domain>=66 && domain<80){
    if(domain%2==0) return 80-domain;
    if(domain%2==1) return 82-domain;
  }
  if(domain>=82 && domain<84){
    if(domain%2==0) return 82-domain;
    if(domain%2==1) return 84-domain;
  }
  //f_seg1
  if(domain>=84 && domain<96){
    if(domain%2==0) return domain-83;
    if(domain%2==1) return domain-85;
  }
  if(domain>=98 && domain<102){
    if(domain%2==0) return domain-85;
    if(domain%2==1) return domain-87;
  }
  //f_seg2
  if(domain>=102 && domain<112){
    if(domain%2==0) return domain-101;
    if(domain%2==1) return domain-103;
  }
  if(domain>=114 && domain<120){
    if(domain%2==0) return domain-103;
    if(domain%2==1) return domain-105;
  }
  //f_seg3
  if(domain>=120 && domain<128){
    if(domain%2==0) return 134-domain;
    if(domain%2==1) return 136-domain;
  }
  if(domain>=130 && domain<138){
    if(domain%2==0) return 136-domain;
    if(domain%2==1) return 138-domain;
  }
  //r
  if(domain>=138 && domain<142){
    if(domain%2==0) return domain-125;
    if(domain%2==1) return domain-127;
  }  
  if(domain>=142 && domain<144){
    return 143-domain;
  }  
  if(domain>=146 && domain<152){
    if(domain%2==0) return domain-143;
    if(domain%2==1) return domain-145;
  }
  //test
  if(domain==152)  return 9;
  if(domain==153)  return 11;
  if(domain==154)  return 8;
  if(domain==155)  return 10;

  return -1;
}

int entry_type(int board, int ch){
  int type_n=-1; // 0:trigger, 10:Gamma_signal, 11:BGO_signal, 20+n:Si_f_signal(n=seg), 30+n:Si_r_signal(n=seg), 40:scaler
  
  if(board==0 && ch==0) type_n = 0;
  if(board==0 && ch>1){
    if(ch%2==0) type_n = 10;
    if(ch%2==1) type_n = 11;
  }
  if(board==1 || board==2 || board==3){
    if(ch%2==0) type_n = 10;
    if(ch%2==1) type_n = 11;
    if(ch>11) type_n = 40;
  }
  //f_seg0
  if(board==4 && ch>1) type_n = 20;
  if(board==5 && ch>1 && ch<4) type_n = 20;
  //f_seg1
  if(board==5 && ch>3) type_n = 21;
  if(board==6 && ch>1 && ch<6) type_n = 21;
  //f_seg2
  if(board==6 && ch>5) type_n = 22;
  if(board==7 && ch>1 && ch<8) type_n = 22;
  //f_seg3
  if(board==7 && ch>7) type_n = 23;
  if(board==8 && ch>1 && ch<10) type_n = 23;
  //r_seg0
  if(board==8 && ch>13) type_n = 30;
  if(board==9 && ch>1 && ch<4) type_n = 30;
  //r_seg1
  if(board==9 && ch>3 && ch<8) type_n = 31;
  //r_seg2
  if(board==9 && ch>7 && ch<12) type_n = 32;
  //r_seg3
  if(board==8 && ch>9 && ch<14) type_n = 33;
  
  return type_n;
}; 

float get_front_theta(int front_ch){
  
  float strp_wid = (s1_r2 - s1_r1)/16.0;
  float r=0;  
  r = s1_r1 + strp_wid*(0.5) + strp_wid*(15-front_ch);
  
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



std::pair<float,float> get_angle(float x, float y, int ch_f, int ch_r){

  float tmp_x, tmp_y, tmp_z, tmp_theta, tmp_phi;
  float front_theta = get_front_theta(ch_f);
  float rear_phi = (90 + 22.5/2.0 + 22.5*ch_r)*PI/180;

  tmp_x = s1_dist * tan(front_theta)*cos(rear_phi) - x;
  tmp_y = s1_dist * tan(front_theta)*sin(rear_phi) - y;
  tmp_z = s1_dist;
  tmp_theta = acos(tmp_z/pow(pow(tmp_x,2)+pow(tmp_y,2)+pow(tmp_z,2),0.5));
  if(tmp_y>=0) tmp_phi = acos(tmp_x/pow(pow(tmp_x,2)+pow(tmp_y,2),0.5));
  if(tmp_y<0) tmp_phi = -1 * acos(tmp_x/pow(pow(tmp_x,2)+pow(tmp_y,2),0.5));
  
  return std::make_pair(tmp_theta,tmp_phi);
}




void fill_data(int seg, int num_fi, int num_ff, int num_ri, int num_rf){
  Energy_f[seg][num_ff]=ene_f[seg][num_fi];
  Energy_r[seg][num_rf]=ene_r[seg][num_ri];
  Amax[seg][num_ff]=Amax_f[seg][num_fi];
  ts_diff_f[seg][num_ff]=ts_f[seg][num_fi];
  ts_diff_r[seg][num_rf]=ts_r[seg][num_ri];
  ch_f[seg][num_ff]=get_ch(domain_f[seg][num_fi]);
  ch_r[seg][num_rf]=get_ch(domain_r[seg][num_ri]);

  double tmp_theta1 = get_angle(0, 2, ch_f[seg][num_ff], ch_r[seg][num_rf]).first;
  cor_ex[seg][num_ff] = get_front_ex(tmp_theta1, Energy_f[seg][num_ff]);
}


double detect3a(double *enea, int *chfa, int *chra){

  double m4He=931.494*4+2.4249;
  double m12C=931.494*12;
  double tha[3]={};
  double pha[3]={};
  double moma[3][3]={};
  double mom12C[3]={0,0,0};
  double total12C=-1;
  double ex12C=-1;
  
  for(int n=0; n<3; n++){
    tha[n] = get_angle(0,2,chfa[n],chra[n]).first;
    pha[n] = get_angle(0,2,chfa[n],chra[n]).second;
    //    cout << tha[n] << " " << pha[n] << endl;
    //    cout << tha[n]*180/3.14 << " " << pha[n]*180/3.14 << endl;
    
    moma[n][0] = pow(2*m4He*enea[n]+enea[n]*enea[n],0.5) * sin(tha[n]) * cos(pha[n]);
    moma[n][1] = pow(2*m4He*enea[n]+enea[n]*enea[n],0.5) * sin(tha[n]) * sin(pha[n]);
    moma[n][2] = pow(2*m4He*enea[n]+enea[n]*enea[n],0.5) * cos(tha[n]);
    //    cout << enea[n] << " " << moma[n][0]  << " " << moma[n][1] << " "  << moma[n][2] << endl;
    //    cout << m4He+enea[n] << " " << moma[n][0]  << " " << moma[n][1] << " "  << moma[n][2] << endl;

    mom12C[0] += moma[n][0];
    mom12C[1] += moma[n][1];
    mom12C[2] += moma[n][2];
  }
  //  cout << mom12C[0]  << " " << mom12C[1] << " "  << mom12C[2] << endl;

  total12C = pow(mom12C[0]*mom12C[0]+mom12C[1]*mom12C[1]+mom12C[2]*mom12C[2]+m12C*m12C,0.5);
  ex12C = enea[0]+enea[1]+enea[2]+m4He*3 - total12C;
  
  return ex12C;
}
