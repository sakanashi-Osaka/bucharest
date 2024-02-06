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
#include <TF1.h>
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

using namespace std;

const double s1_r1 = 48/2.0;  // inner radius      
const double s1_r2 = 96/2.0;  // outer radius
const double s1_dist = 40;

double get_front_theta(int front_ch);
std::pair<double,double> get_angle(double x, double y, int ch_f, double ch_r);
double get_front_ex(double theta, double ene);
double detect3a(double *enea, int *chfa, int *chra);
std::pair<double,double> thetaCM_to_phiLAB(double thetaCM);
int getStrFromText(string filename, vector<string> &vstr){
  ifstream ifs(filename);
  string tmp;
  while (getline(ifs, tmp))
    vstr.push_back(tmp);
  return 0;
}


int hit_n[4]={0,0,0,0};
int ch_f[4][N_HIT_MAX]={};
int ch_r[4][N_HIT_MAX]={};
double cor_ex[4][N_HIT_MAX]={};


int count(int run){
//int count(){
  
  //  ofstream ofs(Form("log%d.txt",run),std::ios::out);
  
  
  int run_tmp[300];
  double beam_x, beam_y, px[300], py[300];
  vector<string> vstr;
  getStrFromText("../beam_pos.prm",vstr);
  ifstream ifs3("../beam_pos.prm");
  for(int i=0; i<(int)vstr.size(); i++){
    ifs3 >> run_tmp[i] >> px[i] >> py[i];
  }
  for(int i=0; i<(int)vstr.size(); i++){
    if(run_tmp[i]==run && abs(px[i]-20)<10 && abs(py[i]-30)<10){ 
      beam_x = (px[i]-20.0) * 0.2;
      beam_y = (py[i]-20.0) * 0.2;
      break;
    }else{
      beam_x = 0.;
      beam_y = 2.;
    }
  }
  cout <<"beam pos: " << beam_x << " " << beam_y << endl;

  
  TFile *fin =new TFile(Form("../rootfile/hit%d.root",run));
  TTree *hit = (TTree*)fin->Get("hit");
  
  
  TH1F *h[4][16];
  for(int i=0;i<4;i++){
    for(int j=0;j<16;j++){
      h[i][j] = new TH1F(Form("h_%d_%d",i,j),Form("h_%d_%d",i,j),100,7,8);
    }
  }
  TF1 *f[4][16];
  for(int i=0;i<4;i++){
    for(int j=0;j<16;j++){
      f[i][j] = new TF1(Form("f_%d_%d",i,j),"[0]+[1]*x+[2]*exp(-pow(x-[3],2)/(2*pow([4],2)))+[5]*exp(-pow(x-[6],2)/(2*pow([7],2)))",7,8);
    }
  }
  
  
  hit->SetBranchAddress("hit_n",hit_n);
  hit->SetBranchAddress("ch_f",ch_f);
  hit->SetBranchAddress("ch_r",ch_r);
  hit->SetBranchAddress("cor_ex",cor_ex);
  
  
  ULong64_t N=hit->GetEntries();
  cout << "Total entry: " << N << endl;
  for(ULong64_t evtn=0; evtn<N; evtn++){
  //  for(ULong64_t evtn=0; evtn<1000000; evtn++){
      if(evtn%10000==0) cout << "\rAnalyzed entry:" << evtn; std::cout << flush;
      
    hit->GetEntry(evtn); 
    
    for(int i=0; i<4; i++){
      
      if(hit_n[i]==11 && cor_ex[i][0]<8 && cor_ex[i][0]>7){
	h[i][ch_f[i][0]]->Fill(cor_ex[i][0]);
      }
    } 
  }
  
  for(int i=0;i<4;i++){
    for(int j=0;j<16;j++){
      double tmp0 = h[i][j]->GetBinContent(1);
      double tmp1 = h[i][j]->GetBinContent(99);
      double tmp2 = h[i][j]->GetMaximum();
      double tmp3 = h[i][j]->GetMaximumBin();
      double tmp4 = h[i][j]->GetXaxis()->GetXmin();
      double tmp5 = h[i][j]->GetXaxis()->GetXmax();
      double tmp6 = h[i][j]->GetXaxis()->GetNbins();
      double tmp7 = tmp4 + (tmp5-tmp4)/tmp6*tmp3;
      
      /*
      f[i][j]->SetParameter(0,tmp0+(tmp1-tmp0)*4);
      f[i][j]->SetParameter(1,(tmp1-tmp0)/2);
      f[i][j]->SetParameter(2,tmp2-tmp0);
      f[i][j]->SetParLimits(2,(tmp2-tmp0)/2,(tmp2-tmp0)*1.5);
      f[i][j]->SetParameter(3,7.2); ///
      f[i][j]->SetParLimits(3,6.9,7.7); ///
      f[i][j]->SetParameter(4,0.1);
      f[i][j]->SetParLimits(4,0.05,0.2);
      */

      f[i][j]->SetParameter(0,tmp0+(tmp1-tmp0)*4);
      f[i][j]->SetParameter(1,(tmp1-tmp0)/2);
      f[i][j]->SetParameter(2,tmp2-tmp0);
      f[i][j]->SetParLimits(2,(tmp2-tmp0)/2,(tmp2-tmp0)*1.5);
      f[i][j]->SetParameter(3,tmp7); ///
      f[i][j]->SetParLimits(3,tmp7-0.05,tmp7+0.05); ///
      f[i][j]->SetParameter(4,0.1);
      f[i][j]->SetParLimits(4,0.05,0.2);
      f[i][j]->SetParameter(5,2000);
      f[i][j]->SetParLimits(5,0,5000);
      f[i][j]->SetParameter(6,tmp7+0.3); ///
      f[i][j]->SetParLimits(6,tmp7+0.1,tmp7+0.8); ///
      f[i][j]->SetParameter(7,0.1);
      f[i][j]->SetParLimits(7,0.05,0.2);
      h[i][j]->Fit(f[i][j]);

      double p0 = f[i][j]->GetParameter(0);
      double p1 = f[i][j]->GetParameter(1);
      double p2 = f[i][j]->GetParameter(2);
      double p3 = f[i][j]->GetParameter(3);
      double p4 = f[i][j]->GetParameter(4);
      double p5 = f[i][j]->GetParameter(5);
      double p6 = f[i][j]->GetParameter(6);
      double p7 = f[i][j]->GetParameter(7);
      double nev0 = p2 * pow(2*3.1416*p4*p4,0.5) * 50;
      double nev1 = p5 * pow(2*3.1416*p7*p7,0.5) * 50;
      //      ofs << i <<" "<< j <<" "<< nev0 << " " << nev1 <<" "<< 180.0/3.141592*get_angle(beam_x, beam_y, j, (double)i*4+1.5).first  << " " << p0 <<" "<< p1 <<" "<< p2 <<" "<< p3 <<" "<< p4 << " " << p5 << " " << p6 << " " << p7 << endl;
    }
  }

  return 0;
}




double get_front_theta(int front_ch){
  
  double strp_wid = (s1_r2 - s1_r1)/16.0;
  double r=0;  
  r = s1_r1 + strp_wid*(0.5) + strp_wid*(15-front_ch);
  
  return atan(r/s1_dist);
}

double get_front_ex(double theta, double ene){
  const double beam_ene = 25.0;
  const double AMU = 931.4943;
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

std::pair<double,double> get_angle(double x, double y, int ch_f, double ch_r){

  double tmp_x, tmp_y, tmp_z, tmp_theta, tmp_phi;
  double front_theta = get_front_theta(ch_f);
  double rear_phi = (90 + 22.5/2.0 + 22.5*ch_r)*PI/180;

  tmp_x = s1_dist * tan(front_theta)*cos(rear_phi) - x;
  tmp_y = s1_dist * tan(front_theta)*sin(rear_phi) - y;
  tmp_z = s1_dist;
  tmp_theta = acos(tmp_z/pow(pow(tmp_x,2)+pow(tmp_y,2)+pow(tmp_z,2),0.5));
  if(tmp_y>=0) tmp_phi = acos(tmp_x/pow(pow(tmp_x,2)+pow(tmp_y,2),0.5));
  if(tmp_y<0) tmp_phi = -1 * acos(tmp_x/pow(pow(tmp_x,2)+pow(tmp_y,2),0.5));
  
  return std::make_pair(tmp_theta,tmp_phi);
}

std::pair<double,double> thetaCM_to_phiLAB(double thetaCM){
  double R_TO_D = 180.0/PI;
  double D_TO_R = PI/180.0;

  double thetarCM = thetaCM * D_TO_R;
  
  double m1 = 4.;
  double m2 = 12.;
  double k1 = 25.;
  double q_value = 7.65;

  double tmp_A = pow( 1+(1+m1/m2)*q_value/k1 ,0.5);
  double theta = atan(sin(thetarCM)/(cos(thetarCM)+m1/m2/tmp_A)) * R_TO_D;
  double phi = atan(sin(thetarCM)/(-cos(thetarCM)+1/tmp_A)) * R_TO_D;

  return std::make_pair(theta,phi);
}


