#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <TCanvas.h>
#include <TEventList.h>
#include <cstdint>
#include <vector>
#include <signal.h>
#include <stdlib.h>
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

double cor_ts(double chf4He);
double get_front_theta(double front_ch);
std::pair<double,double> get_angle(double ch_f, double ch_r);
double get_front_ex(double theta, double ene);
int get_g_slot(int ch_g);
std::pair<double,double> get_g_angle(int g_slot);

int getStrFromText(string filename, vector<string> &vstr){
  ifstream ifs(filename);
  string tmp;
  while (getline(ifs, tmp))
    vstr.push_back(tmp);
  return 0;
}

const double s1_r1 = 48/2.0;  // inner radius      
const double s1_r2 = 96/2.0;  // outer radius
double s1_dist = 40;

//for beam position correction
int run_tmp[300]={};
double px[300]={};
double py[300]={};
double beam_x=0;
double beam_y=0;


int main(int argc, char *argv[]){
  
  TFile *fin =new TFile(Form("run%d.root",atoi(argv[1])));
  TTree *hit = (TTree*)fin->Get("hit");

  if(atoi(argv[1])>=2305 && atoi(argv[1])<=2318) s1_dist=41.0;
  if(atoi(argv[1])>=2319 && atoi(argv[1])<=2327) s1_dist=40.5;
  if(atoi(argv[1])>=2328 && atoi(argv[1])<=2337) s1_dist=40.5;
  if(atoi(argv[1])>=2338 && atoi(argv[1])<=2345) s1_dist=41.0;
  if(atoi(argv[1])>=2346 && atoi(argv[1])<=2362) s1_dist=41.5;
  if(atoi(argv[1])>=2363 && atoi(argv[1])<=2376) s1_dist=40.5;
  if(atoi(argv[1])>=2377 && atoi(argv[1])<=2394) s1_dist=41.0;
  if(atoi(argv[1])>=2394 && atoi(argv[1])<=2413) s1_dist=40.5;

  double chr[16][16]={};
  double chf[16][16]={};
  double pr[16][16]={};
  double pf[16][16]={};
  ifstream ifs("calib_ts/ts.prm");
  for(int i=0; i<16; i++){
    for(int j=0; j<16; j++){
      ifs >> chf[i][j] >> chr[i][j] >> pr[i][j] >> pf[i][j];
      if(pr[i][j]<0||pr[i][j]>100) pr[i][j] = 50.0;
      if(pf[i][j]<0||pf[i][j]>100) pf[i][j] = 50.0;
    }
  }

  double pHe0[16][16]={};
  double pHe1[16][16]={};
  double pHe2[16][16]={};
  double pHe3[16][16]={};
  double pHe4[16][16]={};
  ifstream ifs1(Form("calib_Amax/4He/test%d.txt",atoi(argv[1])));
  for(int i=0; i<16; i++){
    for(int j=0; j<16; j++){
      ifs1 >> pHe0[i][j] >> pHe1[i][j] >> pHe2[i][j] >> pHe3[i][j] >> pHe4[i][j];
    }
  }

  double pa0[16][16]={};
  double pa1[16][16]={};
  double pa2[16][16]={};
  double pa3[16][16]={};
  double pa4[16][16]={};
  ifstream ifs2(Form("calib_Amax/12C/gomi%d.txt",atoi(argv[1])));
  for(int i=0; i<16; i++){
    for(int j=0; j<16; j++){
      ifs2 >> pa0[i][j] >> pa1[i][j] >> pa2[i][j] >> pa3[i][j];
    }
  }

  vector<string> vstr;
  getStrFromText("beam_pos.prm",vstr);
  ifstream ifs3("beam_pos.prm");
  for(int i=0; i<(int)vstr.size(); i++){
    ifs3 >> run_tmp[i] >> px[i] >> py[i];
  }
  for(int i=0; i<(int)vstr.size(); i++){
    if(run_tmp[i]==atoi(argv[1]) && abs(px[i])<1 && abs(py[i]-2)<1){
      beam_x = px[i];
      beam_y = py[i];
      //      beam_x = -0.5;
      //      beam_y = 2.4;
      break;
    }else{
      beam_x = -0.5;
      beam_y = 2.4;
    }
  }
  
  int hit_n[4]={0,0,0,0};
  int ch_f[4][N_HIT_MAX]={};
  int ch_r[4][N_HIT_MAX]={};
  int ch_g[N_HIT_MAX]={};
  int ch_b[N_HIT_MAX]={};
  double Energy_f[4][N_HIT_MAX]={};
  double Energy_r[4][N_HIT_MAX]={};
  double cor_ex[4][N_HIT_MAX]={};
  double ts_diff_f[4][N_HIT_MAX]={};
  double ts_diff_r[4][N_HIT_MAX]={};
  double ts_diff_g[N_HIT_MAX]={};
  double ts_diff_b[N_HIT_MAX]={};
  //  double TDC_f[4][N_HIT_MAX]={};
  //  double TDC_r[4][N_HIT_MAX]={};
  double Amax[4][N_HIT_MAX]={};
  double Gamma[N_HIT_MAX]={};
  double BGO[N_HIT_MAX]={};

  
  hit->SetBranchAddress("hit_n",hit_n);
  hit->SetBranchAddress("ch_f",ch_f);
  hit->SetBranchAddress("ch_r",ch_r);
  hit->SetBranchAddress("ch_g",ch_g);
  hit->SetBranchAddress("ch_b",ch_b);
  hit->SetBranchAddress("Energy_f",Energy_f);
  hit->SetBranchAddress("Energy_r",Energy_r);
  hit->SetBranchAddress("cor_ex",cor_ex);
  hit->SetBranchAddress("ts_diff_f",ts_diff_f);
  hit->SetBranchAddress("ts_diff_r",ts_diff_r);
  hit->SetBranchAddress("ts_diff_g",ts_diff_g);
  hit->SetBranchAddress("ts_diff_b",ts_diff_b);
  //  hit->SetBranchAddress("TDC_f",TDC_f);
  //  hit->SetBranchAddress("TDC_r",TDC_r);
  hit->SetBranchAddress("Amax",Amax);
  hit->SetBranchAddress("Gamma",Gamma);
  hit->SetBranchAddress("BGO",BGO);
  
  
  TFile *fout =new TFile(Form("gamma_angle%d.root",atoi(argv[1])),"recreate");
  TTree *tree = new TTree("tree","tree"); 

  //branch
  double K4He;
  double Ex4He;
  double chf4He;
  double chr4He;
  double tsf4He;
  double tsr4He;
  double Amax4He;
  double th4He;
  double ph4He;
  
  double K12C;
  double chf12C;
  double chr12C;
  double tsf12C;
  double tsr12C;
  double Amax12C;
  double th12C;
  double ph12C;

  double Kg[2];
  double Kg_cor;
  int chg[2];
  double thg[2];
  double phg[2];
  double tsg[2];
  double Kb;
  int chb;
  double thb;
  double phb;
  double tsb;
  double g_angle;
  
  int C_E=0;
  int C_th=0;
  int C_ph=0;
  int C_ts=0;

  int run_n;
  
  tree->Branch("Ex4He",&Ex4He,"Ex4He/D");
  tree->Branch("K4He",&K4He,"K4He/D");
  tree->Branch("chf4He",&chf4He,"chf4He/D");
  tree->Branch("chr4He",&chr4He,"chr4He/D");
  tree->Branch("tsf4He",&tsf4He,"tsf4He/D");
  tree->Branch("tsr4He",&tsr4He,"tsr4He/D");
  tree->Branch("Amax4He",&Amax4He,"Amax4He/D");
  tree->Branch("th4He",&th4He,"th4He/D");
  tree->Branch("ph4He",&ph4He,"ph4He/D");
  
  tree->Branch("K12C",&K12C,"K12C/D");
  tree->Branch("chf12C",&chf12C,"chf12C/D");
  tree->Branch("chr12C",&chr12C,"chr12C/D");
  tree->Branch("tsf12C",&tsf12C,"tsf12C/D");
  tree->Branch("tsr12C",&tsr12C,"tsr12C/D");
  tree->Branch("Amax12C",&Amax12C,"Amax12C/D");
  tree->Branch("th12C",&th12C,"th12C/D");
  tree->Branch("ph12C",&ph12C,"ph12C/D");
  tree->Branch("g_angle",&g_angle,"g_angle/D");

  tree->Branch("C_E",&C_E,"C_E/I");
  tree->Branch("C_th",&C_th,"C_th/I");
  tree->Branch("C_ph",&C_ph,"C_ph/I");
  tree->Branch("C_ts",&C_ts,"C_ts/I");
  
  tree->Branch("Kg",Kg,"Kg[2]/D");
  tree->Branch("Kg_cor",&Kg_cor,"Kg_cor/D");
  tree->Branch("chg",chg,"chg[2]/I");
  tree->Branch("thg",thg,"thg[2]/D");
  tree->Branch("phg",phg,"phg[2]/D");
  tree->Branch("tsg",tsg,"tsg[2]/D");
  tree->Branch("Kb",&Kb,"Kb/D");
  tree->Branch("chb",&chb,"chb/I");
  tree->Branch("thb",&thb,"thb/D");
  tree->Branch("phb",&phb,"phb/D");
  tree->Branch("tsb",&tsb,"tsb/D");

  tree->Branch("run",&run_n,"run/I");
  
  ULong64_t N=hit->GetEntries();
  cout << "Total entry: " << N << endl;
  for(ULong64_t evtn=0; evtn<N; evtn++){
  //  for(ULong64_t evtn=0; evtn<1e6; evtn++){

    if(evtn%100000==0) cout << "\rAnalyzed entry:" << evtn; std::cout << flush;


    
    K4He = -100; chf4He = -100; chr4He = -100; tsf4He = -100; tsr4He = -100;
    Amax4He = -100; th4He = -1000; ph4He = -1000; Ex4He = -100; 
    K12C = -100; chf12C = -100; chr12C = -100; tsf12C = -100; tsr12C = -100;
    Amax12C = -100; th12C = -1000; ph12C = -1000;
    g_angle = -100;
    
    for(int i=0; i<2; i++){    
      Kg[i] = -100; chg[i] = -100; thg[i] = -1000; phg[i] = -1000; tsg[i] = -100; Kg_cor = -100;
    }
    Kb = -100; chb = -100; thb = -1000; phb = -1000; tsb = -100;
    C_E=0; C_th=0; C_ph=0; C_ts=0;

    hit->GetEntry(evtn); 
    run_n=atoi(argv[1]);
    
    int seg_i=-1; // selected high energy 4He 
    double ene_max=0;
    for(int i=0; i<4; i++){
      for(int j=0; j<N_HIT_MAX; j++){
	//	cout << Energy_f[i][j] << endl;
	if(ene_max < Energy_f[i][j]){
	  ene_max = Energy_f[i][j];
	  seg_i = i;
	}
      }
    }
    //    cout << evtn << " " << ene_max << " " << seg_i << endl;
    if(seg_i<0) continue;
      
    double tmp_theta = -1000;
    double tmp_phi = -1000;
    double tmp_ene = -100;
    double tmp_ex = -100;
    double tmp_chf = -100;
    double tmp_chr = -100;
    
    if(hit_n[seg_i]==11){
      tmp_chf = (double)ch_f[seg_i][0];
      tmp_chr = (double)ch_r[seg_i][0];
      tmp_theta = get_angle(tmp_chf,tmp_chr).first;
      tmp_phi = get_angle(tmp_chf,tmp_chr).second;
      tmp_ene = Energy_f[seg_i][0];
      tmp_ex = get_front_ex(tmp_theta,tmp_ene);
      //      tmp_ex = cor_ex[seg_i][0];

    }else if(hit_n[seg_i]==21 && abs(ch_f[seg_i][0]-ch_f[seg_i][1])<=1){
      continue;
      tmp_chf = (double)(ch_f[seg_i][0]+ch_f[seg_i][1])/2;
      tmp_chr = (double)ch_r[seg_i][0];
      tmp_theta = get_angle(tmp_chf,tmp_chr).first;
      tmp_phi = get_angle(tmp_chf,tmp_chr).second;
      tmp_ene = Energy_f[seg_i][0] + Energy_f[seg_i][1];
      tmp_ex = get_front_ex(tmp_theta, tmp_ene);	
      //      tmp_ex = cor_ex[seg_i][0];
    }else{
      continue;
    }

    if(4.1 < tmp_ex && tmp_ex<4.7){  // selected leading alpha //gomi*.root
      Ex4He = tmp_ex;
      K4He = tmp_ene; 
      chf4He = tmp_chf;
      chr4He = tmp_chr;
      tsf4He = ts_diff_f[seg_i][0]/1e3;
      tsr4He = ts_diff_r[seg_i][0]/1e3;
      Amax4He = Amax[seg_i][0] - pHe3[(int)chf4He][(int)chr4He]; 
      th4He = tmp_theta *180/PI;
      ph4He = tmp_phi *180/PI;

      for(int i=0; i<2; i++){
	Kg[i] = Gamma[i];
	chg[i] = get_g_slot(ch_g[i]);
	thg[i] = get_g_angle(chg[i]).first;
	phg[i] = get_g_angle(chg[i]).second;
	tsg[i] = ts_diff_g[i]/1e3;
      }
      Kb = BGO[0];
      chb = get_g_slot(ch_b[0]-1);
      thb = get_g_angle(chb).first;
      phb = get_g_angle(chb).second;
      tsb = ts_diff_b[0]/1e3;
	


      int seg_j = (seg_i+2)%4; // opposite segment

      //      if(hit_n[seg_j]==11){ 
      tmp_chf = (double)ch_f[seg_j][0];
      tmp_chr = (double)ch_r[seg_j][0];
      tmp_theta = get_angle(tmp_chf,tmp_chr).first;
      tmp_phi = get_angle(tmp_chf,tmp_chr).second;
      
      K12C = Energy_f[seg_j][0]; 
      chf12C = ch_f[seg_j][0];
      chr12C = ch_r[seg_j][0];
      tsf12C = ts_diff_f[seg_j][0]/1e3 + K4He*40/18 -pf[(int)chf12C][(int)chr12C];
      tsr12C = ts_diff_r[seg_j][0]/1e3 + K4He*40/18 -pr[(int)chf12C][(int)chr12C];
      Amax12C = Amax[seg_j][0] - (pa0[(int)chf12C][(int)chr12C] + pa1[(int)chf12C][(int)chr12C]*K12C + pa2[(int)chf12C][(int)chr12C]/(K12C-pa3[(int)chf12C][(int)chr12C]));
      //      th12C = tmp_theta *180/PI;
      //      ph12C = tmp_phi *180/PI;


      if(Kg[0]<0.5) continue;
      th12C = 75. - 0.5*th4He;
      ph12C = 180. + ph4He;      
      if(ph12C>360.) ph12C = ph12C-360.; 
      //      cout << "4He " << th4He << " " << ph4He << " " << K4He << endl;
      //      cout << "12C " << th12C << " " << ph12C << " " << K12C << endl;
      //      cout << "gamma " << thg[0] << " " << phg[0] << " " << Kg[0] << endl;
      //      cout << endl;
      // angle between a-> and b-> (12C && gamma)
      // ang_cos = a*b / |a|*|b| = (ax*bx + ay*by + az*bz) / |a|*|b|
      double tmp_xx = sin(th12C*PI/180)*cos(ph12C*PI/180) * sin(thg[0]*PI/180)*cos(phg[0]*PI/180);
      double tmp_yy = sin(th12C*PI/180)*sin(ph12C*PI/180) * sin(thg[0]*PI/180)*sin(phg[0]*PI/180);
      double tmp_zz = cos(th12C*PI/180) * cos(thg[0]*PI/180);
      double ang_cos = tmp_xx+tmp_yy+tmp_zz;
      double tmp_angle = acos(ang_cos)*180/PI;

      g_angle = tmp_angle + 10*(0.5-(double)rand()/RAND_MAX);
      
      /*
      double totalE_12C = K12C+AMU*12+7.65;
      double tmpmom = sqrt(totalE_12C*totalE_12C-(AMU*12+7.65)*(AMU*12+7.65));
      double tmpgam = pow(1+pow(tmpmom/(AMU*12),2),0.5);
      double tmpbeta = pow(1-1/pow(tmpgam,2),0.5);
      Kg_cor = Kg[0]*(1-ang_cos*tmpbeta)/pow(1-pow(tmpbeta,2),0.5);
      */
      
	//      }
      tree->Fill(); //gomi*.root 
    }
  }
  
  tree->AutoSave();
  fout->Close();
  
  return 0;
}



int get_g_slot(int ch_g){
  int g_slot=-1; //{A,B,C,-B,-A}={0,10,20,30,40}
  if(ch_g==14) g_slot = 1;
  if(ch_g==2) g_slot = 2;
  if(ch_g==4) g_slot = 3;
  if(ch_g==6) g_slot = 4;
  if(ch_g==8) g_slot = 5;

  if(ch_g==10) g_slot = 11;
  if(ch_g==12) g_slot = 12;
  if(ch_g==16) g_slot = 13;
  if(ch_g==18) g_slot = 14;
  if(ch_g==20) g_slot = 15;

  if(ch_g==22) g_slot = 21;
  if(ch_g==24) g_slot = 25;
  if(ch_g==26) g_slot = 24;
  if(ch_g==32) g_slot = 23;
  if(ch_g==34) g_slot = 22;

  if(ch_g==36) g_slot = 31;
  if(ch_g==38) g_slot = 32;
  if(ch_g==40) g_slot = 33;
  if(ch_g==42) g_slot = 34;
  if(ch_g==48) g_slot = 35;

  if(ch_g==50) g_slot = 41;
  if(ch_g==52) g_slot = 45;
  if(ch_g==54) g_slot = 44;
  if(ch_g==56) g_slot = 43;
  if(ch_g==58) g_slot = 42;
  
  return g_slot;
}

std::pair<double,double> get_g_angle(int g_slot){
  double tmp_theta=-1000, tmp_phi=-1000;
  if(00<g_slot && g_slot<10) tmp_theta=37.0;
  if(10<g_slot && g_slot<20) tmp_theta=70.0;
  if(20<g_slot && g_slot<30) tmp_theta=90.0;
  if(30<g_slot && g_slot<40) tmp_theta=110.0;
  if(40<g_slot && g_slot<50) tmp_theta=143.0;

  if((0<g_slot&&g_slot<10) || (20<g_slot&&g_slot<30) || (40<g_slot&&g_slot<50)){
    tmp_phi = -90.+72*(g_slot%10-1);
  }
  if((10<g_slot&&g_slot<20) || (30<g_slot&&g_slot<40)){
    tmp_phi = -54.+72*(g_slot%10-1);
  }
  if(tmp_phi>180) tmp_phi=tmp_phi-360.;
    
  return std::make_pair(tmp_theta,tmp_phi);
}


double cor_ts(double chf4He){
  double p0 = 16.973;
  double p1 = 0.117552;
  
  double output = p0 - p1*pow(chf4He-13,2);
  return output;
}

double get_front_theta(double front_ch){
  
  double strp_wid = (s1_r2 - s1_r1)/16.0;
  double r=0;  
  r = s1_r1 + strp_wid*(0.5) + strp_wid*(15-front_ch);
  
  return atan(r/s1_dist);
}

std::pair<double,double> get_angle(double ch_f, double ch_r){
  
  double tmp_x, tmp_y, tmp_z, tmp_theta, tmp_phi;
  double front_theta = get_front_theta(ch_f);
  double rear_phi = (90 + 22.5/2.0 + 22.5*ch_r)*PI/180;

  tmp_x = s1_dist * tan(front_theta)*cos(rear_phi) - beam_x;
  tmp_y = s1_dist * tan(front_theta)*sin(rear_phi) - beam_y;
  tmp_z = s1_dist;
  tmp_theta = acos(tmp_z/pow(pow(tmp_x,2)+pow(tmp_y,2)+pow(tmp_z,2),0.5));
  if(tmp_y>=0) tmp_phi = acos(tmp_x/pow(pow(tmp_x,2)+pow(tmp_y,2),0.5));
  if(tmp_y<0) tmp_phi = -1 * acos(tmp_x/pow(pow(tmp_x,2)+pow(tmp_y,2),0.5));
  
  return std::make_pair(tmp_theta,tmp_phi);
}

double get_front_ex(double theta, double ene){
  const double beam_ene = 25.0;
  const double mass_ex_12c = 0;
  const double mass_ex_4he = 2.425;
  const double mass_ex_13c = 3.125;    

  double mass_12c = AMU*12 + mass_ex_12c;
  double mass_4he = AMU*4  + mass_ex_4he;
  double mass_13c = AMU*13 + mass_ex_13c;  
  
  double m1, m2, m3, m4;
  double E1, E3;
  double p1, p3;

  double s,t,u;
  double total_m4;
  double ex4;
  
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
