#define ch 32
#define peak_n 1

int count()
{
  
  TSpectrum *s[peak_n];
  TH1F *h[ch];
  TF1 *f[ch][peak_n];
  double *xpeaks[ch];
  double peak[peak_n];
  double height[ch][peak_n];
  double mean[ch][peak_n];
  double sigma[ch][peak_n];
  int n;

  ofstream ofs("log.txt",std::ios::app);

  for(Int_t i=0;i<ch;i++){
    s[i] = new TSpectrum(peak_n);
    //    h[i] = new TH1F(Form("h_%d",i),Form("h_%d",i),100,4.2,4.6); // for run2060
    h[i] = new TH1F(Form("h_%d",i),Form("h_%d",i),40,3.0,3.4);
    tree->Draw(Form("Energy>>h_%d",i),Form("domain==%d",2*i),"",2e8);
    s[i]->Search(h[i],1,"new",0.1);
    xpeaks[i]=s[i]->GetPositionX();

    for(Int_t j=0;j<peak_n;j++){
      peak[j]=xpeaks[i][j];
    }
    n = sizeof(peak)/sizeof(peak[0]);
    std::sort(peak,peak+n);


    double gomi0 = h[i]->GetBinContent(1);
    //    double gomi1 = h[i]->GetBinContent(100-1); // for run2060
    double gomi1 = h[i]->GetBinContent(40-1);

    //    double tmp1 = (gomi1-gomi0)/(4.6-4.2); //slope // for run2060
    double tmp1 = (gomi1-gomi0)/(3.4-3.0); //slope
    double tmp0 = gomi0-tmp1*3.0; //interseption

    //    double tmp_hight = h[i]->GetBinContent(50); // for run2060
    double tmp_hight = h[i]->GetBinContent(30);
    
    for(Int_t j=0;j<peak_n;j++){
      //      f[i][j]=new TF1(Form("f%d_%d",i,j),"[0]+[1]*x+[2]*exp(-pow(x-[3],2)/(2*pow([4],2)))",peak[j]-40,peak[j]+40);
      f[i][j]=new TF1(Form("f%d_%d",i,j),"[0]+[1]*x+[2]*exp(-pow(x-[3],2)/(2*pow([4],2)))",peak[j]-5,peak[j]+5);
      f[i][j]->SetParameter(0,tmp0);
      f[i][j]->SetParameter(1,tmp1);
      f[i][j]->SetParameter(2,tmp_hight);
      //      f[i][j]->SetParameter(3,(4.2+4.6)/2);
      f[i][j]->SetParameter(3,(3.0+3.4)/2);
      f[i][j]->SetParameter(4,0.03);
    }

    for(Int_t j=0;j<peak_n;j++){
      //      h[i]->Fit(Form("f%d_%d",i,j),"+","",peak[j]-40,peak[j]+40);
      h[i]->Fit(Form("f%d_%d",i,j),"+","",peak[j]-5,peak[j]+5);
      height[i][j]=f[i][j]->GetParameter(2);
      mean[i][j]=f[i][j]->GetParameter(3);
      sigma[i][j]=f[i][j]->GetParameter(4);
      ofs << 2*i << " " << height[i][j] << " " << mean[i][j] << " " << sigma[i][j] << endl;
    }

    h[i]->Draw();
    sleep(5);
  }

  ofs.close();

  return 0;
}

