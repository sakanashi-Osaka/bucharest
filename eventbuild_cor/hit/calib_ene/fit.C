#define ch 256
#define peak_n 3

int fit()
{
  TH1F *h[ch];
  TF1 *f[ch][peak_n];
  double mean[ch][peak_n];
  double sigma[ch][peak_n];
  int n;
  
  ofstream ofs("log.txt",std::ios::app);
  
  //  for(Int_t i=0;i<ch;i++){ 
  for(Int_t i=0;i<256;i++){ 
    h[i] = new TH1F(Form("h_%d",i),Form("h_%d",i),500,-1,14);

    hit->Draw(Form("cor_ex[][0]>>h_%d",i),Form("ch_f[][0]==%d&&ch_r[][0]==%d&&hit_n[]==11",i%16,(i-i%16)/16));


    f[i][0]=new TF1(Form("f%d_%d",i,0),"gaus",-0.5,0.5);
    f[i][1]=new TF1(Form("f%d_%d",i,1),"gaus",4,5);
    f[i][2]=new TF1(Form("f%d_%d",i,2),"gaus",9,10);
    
    h[i]->Fit(Form("f%d_%d",i,0),"+","",-0.5,0.5);
    h[i]->Fit(Form("f%d_%d",i,1),"+","",4,5);
    h[i]->Fit(Form("f%d_%d",i,2),"+","",9,10);
    for(Int_t j=0;j<3;j++){ 
      mean[i][j]=f[i][j]->GetParameter(1);
      sigma[i][j]=f[i][j]->GetParameter(2);
      cout << mean[i][j] << " " << sigma[i][j] << endl;
      ofs << i << " " << j << " " << mean[i][j] << " " << sigma[i][j] << endl;
    }
  }

  ofs.close();
  
  return 0;
}
