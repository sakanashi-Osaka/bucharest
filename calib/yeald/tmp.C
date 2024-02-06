#define ch 2
#define peak_n 1

int tmp()
{


  TSpectrum *s[peak_n];
  TH1D *h[ch];
  TF1 *f[ch][peak_n];
  double *xpeaks[ch];
  double peak[peak_n];
  double height[ch][peak_n];
  double mean[ch][peak_n];
  double sigma[ch][peak_n];
  int n;

  ofstream ofs("log.txt",std::ios::app);

  for(Int_t i=ch;i<ch+1;i++){
    s[i] = new TSpectrum(peak_n);
    h[i] = new TH1D(Form("h_%d",i),Form("h_%d",i),40,4.1,4.5);
    tree->Draw(Form("Energy>>h_%d",i),Form("domain==%d",2*i),"",5e7);
    s[i]->Search(h[i],1,"new",0.1);
    xpeaks[i]=s[i]->GetPositionX();

    for(Int_t j=0;j<peak_n;j++){
      peak[j]=xpeaks[i][j];
    }
    n = sizeof(peak)/sizeof(peak[0]);
    std::sort(peak,peak+n);


    double gomi0 = h[i]->GetBinContent(1);
    double gomi1 = h[i]->GetBinContent(20-1);

    double tmp1 = (gomi1-gomi0)/(4.5-4.1); //slope
    double tmp0 = gomi0-tmp1*2.0; //interseption
    
    for(Int_t j=0;j<peak_n;j++){
      f[i][j]=new TF1(Form("f%d_%d",i,j),"[0]+[1]*x+[2]*exp(-pow(x-[3],2)/(2*pow([4],2)))",peak[j]-4,peak[j]+4);
      f[i][j]->SetParameter(0,tmp0);
      f[i][j]->SetParameter(1,tmp1);
      f[i][j]->SetParameter(2,1500);
      f[i][j]->SetParameter(3,(4.1+4.5)/2);
      f[i][j]->SetParameter(4,0.03);
    }

    for(Int_t j=0;j<peak_n;j++){
      h[i]->Fit(Form("f%d_%d",i,j),"+","",peak[j]-8,peak[j]+8);
      height[i][j]=f[i][j]->GetParameter(2);
      mean[i][j]=f[i][j]->GetParameter(3);
      sigma[i][j]=f[i][j]->GetParameter(4);
      ofs << 2*i << " " << height[i][j] << " " << mean[i][j] << " " << sigma[i][j] << endl;
    }

    h[i]->Draw();
  }

  ofs.close();

  return 0;
}

