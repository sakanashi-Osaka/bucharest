#define CH 92
#define PEAK_N 4

int fit_zero()
{
  double ch[CH][PEAK_N];
  double x[CH][PEAK_N];
  double y[CH][PEAK_N];
  double sigma[CH][PEAK_N];

  double ene_source[PEAK_N]={0,5.139,5.474,5.789};
  
  ifstream ifs("log.txt");
  ofstream ofs("tmp_zero.prm",std::ios::app); 
  
  for(int i=0; i<CH; i++){
    for(int j=1; j<PEAK_N; j++){
      ifs >> ch[i][j] >> x[i][j] >> y[i][j] >> sigma[i][j];
      x[i][0]=0;
    }
  }

  TF1 *f = new TF1("fit","[0]+[1]*x",0,4000);
  TGraph *g1;                                                               

  for(int i=0; i<CH; i++){
    g1 = new TGraph(PEAK_N, y[i], ene_source);
    g1->Draw("AP");                                               
    g1->Fit("fit");             
    double p0 = f->GetParameter(0);
    double p1 = f->GetParameter(1);

    if(i%16==0 || i%16==1){
      ofs << i+64 << " " << 0 << " " << 1 << " " << endl;
    }else{
      ofs << i+64 << " " << p0 << " " << p1 << " " << endl;
    }
  }
  return 0;
}
