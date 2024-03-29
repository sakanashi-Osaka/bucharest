double connect(double input, double *x_array, double *y_array, int num);
double culc_omega(int channel);


int hist13C(){
  
  int pa[64], pb[64];
  double pc[64], pd[64], pe[64], pf[64], pg[64], ph[64], pi[64], pj[64], pk[64], pl[64], pm[64];  
  //j:theta3lab, k:thetacm, l:k3, m:g-factor, n:not used, o:theta4lab, p:k4
  
  ifstream ifs0("./log_13C_740.txt");
  TH2D *h0= new TH2D("h0","h0",100,20,60,100,0,20000);
  TH2D *h1= new TH2D("h1","h1",100,20,60,100,0,20000);
  for(int ch=0;ch<64;ch++){
    ifs0 >> pa[ch] >> pb[ch] >> pc[ch] >> pd[ch] >> pe[ch] >> pf[ch] >> pg[ch] >> ph[ch] >> pi[ch] >> pj[ch] >> pk[ch] >> pl[ch] >> pm[ch];
    h0->Fill(pe[ch],pc[ch]);
    h1->Fill(pe[ch],pd[ch]);
  }
  h0->Draw("col");
  //  h1->Draw("col same");
  
  
  TF1 *f2;
  
  //gaussian
  f2 = new TF1("f2","[0] + [1]*x + [2]*exp(-pow(x-[3],2)/(2*pow([4],2)))",45,70);
  f2->SetParameter(0,0);
  f2->SetParameter(1,0);
  f2->SetParameter(2,5000);
  f2->SetParLimits(2,2000,10000);
  f2->SetParameter(3,37);
  f2->SetParLimits(3,35,40);
  f2->SetParameter(4,5);
  f2->SetParLimits(4,3,10);
  
  h0->Fit(f2);

  //  double x[9]={28.5,31.7,35.0,38.4,42.3,48.8,55.4,62.3,68.8};
  //  double y[9]={49.84,23.62,8.831,4.429,10.77,22.28,11.82,9.388,13.42};
  //  for(int i=0; i<9; i++) y[i]=y[i]*1.5;
//  TGraph *g2 = new TGraph(9,x,y);
  //  g2->Draw("same");

  return 0;
}


double connect(double input, double *x_array, double *y_array, int num){
  int tmp_x = -1;
  double diff0, diff1;
  
  for(int i=0; i<num-1; i++){
    diff0 = x_array[i] - input;
    diff1 = x_array[i+1] - input;

    if(diff0==0) tmp_x = i;
    if(diff0*diff1 < 0) tmp_x = i;
  }

  if(tmp_x < 0){
    std::cout << "input value ls wrong!!" << std::endl;
    return -1;
  }

  double a, b;
  double output;
  a = y_array[tmp_x+1] - y_array[tmp_x];
  b = y_array[tmp_x] - a*x_array[tmp_x];
  output = a*input + b;

  return output;
}

double culc_omega(int channel){
  double PI = 3.141592;
  double Si_r1 = 24.;
  double Si_r2 = 48.;
  double d_Si = 40.;
  
  double omega = -1;
  omega = (pow(Si_r2-1.5*channel,2) - pow(Si_r2-1.5*channel-1.5,2))*PI/4/(pow(d_Si,2)+pow(Si_r2-1.5*channel-1.5,2));
    
  return omega;
}

