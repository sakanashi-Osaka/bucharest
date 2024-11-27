{
  double x[16];
  double num[16];
  double num1[16];

  for(int i=0; i<4; i++){
    num[i] = hit->Draw("",Form("ch_r[0][0]==%d&&ch_f[0][0]<14&&Amax[0][0]<0.035&&abs(cor_ex[0][0]-7.65)<0.2&&hit_n[0]==11",i),"col",2e6);
  }
  for(int i=4; i<8; i++){
    num[i] = hit->Draw("",Form("ch_r[1][0]==%d&&ch_f[1][0]<14&&Amax[1][0]<0.05&&abs(cor_ex[1][0]-7.65)<0.2&&hit_n[1]==11",i),"col",2e6);
  }
  for(int i=8; i<12; i++){
    num[i] = hit->Draw("",Form("ch_r[2][0]==%d&&ch_f[2][0]<14&&Amax[2][0]<0.055&&abs(cor_ex[2][0]-7.65)<0.2&&hit_n[2]==11",i),"col",2e6);
  }
  for(int i=12; i<16; i++){
    num[i] = hit->Draw("",Form("ch_r[3][0]==%d&&ch_f[3][0]<14&&Amax[3][0]<0.04&&abs(cor_ex[3][0]-7.65)<0.2&&hit_n[3]==11",i),"col",2e6);
  }

  for(int i=0; i<4; i++){
    num1[i] = hit->Draw("",Form("ch_r[0][0]==%d&&ch_f[0][0]<14&&Amax[0][0]<0.035&&abs(cor_ex[0][0]-7.65)<0.2&&hit_n[0]==11&&Energy_f[2][0]>1",i),"col",2e6);
  }
  for(int i=4; i<8; i++){
    num1[i] = hit->Draw("",Form("ch_r[1][0]==%d&&ch_f[1][0]<14&&Amax[1][0]<0.05&&abs(cor_ex[1][0]-7.65)<0.2&&hit_n[1]==11&&Energy_f[3][0]>1",i),"col",2e6);
  }
  for(int i=8; i<12; i++){
    num1[i] = hit->Draw("",Form("ch_r[2][0]==%d&&ch_f[2][0]<14&&Amax[2][0]<0.055&&abs(cor_ex[2][0]-7.65)<0.2&&hit_n[2]==11&&Energy_f[0][0]>1",i),"col",2e6);
  }
  for(int i=12; i<16; i++){
    num1[i] = hit->Draw("",Form("ch_r[3][0]==%d&&ch_f[3][0]<14&&Amax[3][0]<0.04&&abs(cor_ex[3][0]-7.65)<0.2&&hit_n[3]==11&&Energy_f[1][0]>1",i),"col",2e6);
  }

  for(int i=0; i<16; i++) x[i]=i;
  auto g = new TGraph(16,x,num);
  auto g1 = new TGraph(16,x,num1);
  g->Draw("");
  g1->Draw("same");
}
