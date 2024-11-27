int culc_Dron(){
  int i,j;
  double k,l;
  std::ifstream fin("log_Dron.txt");
  ofstream ofs("ene_Dron.prm",std::ios::app);

  TH2D *h[256];
  TF1 *f[256];
  
  double mean[3];
  double ene[3]={0.0,4.44,9.64};
  
  for(int m=0; m<256; m++){

    h[m] = new TH2D(Form("h%d",m),"title",500,-1,14,500,-1,14);
    f[m] = new TF1(Form("f%d",m),"[0] + [1]*x",-1,14);
    
    for(int n=0; n<3; n++){
      fin >> i >> j >> k >> l;
      cout << i << " " << j << " " << k << " " << l << endl;
      mean[n]=k;
      h[m]->Fill(mean[n],ene[n]);
    }
    h[m]->Draw("col");
    h[m]->Fit(f[m]);

    double a = f[m]->GetParameter(0);
    double b = f[m]->GetParameter(1);
    ofs << a << " " << b << endl;

  }
  
  return 0;
}
