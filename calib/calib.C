#define ch 96                                                                                                               
#define peak_n 3

int calib()                                                                                                                                                
{                                                                                                                                                         
 TSpectrum *s[peak_n];                                                                                                                                        
  TH1F *h[ch];                                                                                                                                            
  TF1 *f[ch][peak_n];                                                                                                                                         
  double *xpeaks[ch];                                                                                                                                     
  double peak[peak_n];                                                                                                                                        
  double mean[ch][peak_n];                                                                                                                                    
  double sigma[ch][peak_n];                                                                                                                                   
  int n;                                                                                                                                                  
                                                                                                                                                          
  ofstream ofs("log.txt",std::ios::app);                                                                                                                  
                                                                                                                                                          
  for(Int_t i=0;i<ch;i++){                                                                                 
    s[i] = new TSpectrum(peak_n);                                                                                                                             
    h[i] = new TH1F(Form("h_%d",i),Form("h_%d",i),500,1000,4000);                                                                                                                                     
    tree->Draw(Form("ADC>>h_%d",i),Form("domain==%d",i+64),"");                                                              
    s[i]->Search(h[i],1,"new",0.1);                                                                                                                       
    xpeaks[i]=s[i]->GetPositionX();                                                                                                                       
                                                                                                                                                          
    for(Int_t j=0;j<peak_n;j++){                                                                                                                              
      peak[j]=xpeaks[i][j];                                                                                                                               
    }                                                                                                                                                     
    n = sizeof(peak)/sizeof(peak[0]);                                                                                                                     
    std::sort(peak,peak+n);                                                                                                                               
                                                                                                                                                          
    cout << i << " ";                                                                                                                                     
    for(Int_t j=0;j<peak_n;j++){                                                                                                                              
      cout << peak[j] << " ";                                                                                                                             
    }                                                                                                                                                     
    cout << endl;                                                                                                                                         
                                                                                                                                                          
    for(Int_t j=0;j<peak_n;j++){                                                                                                                              
      f[i][j]=new TF1(Form("f%d_%d",i,j),"gaus",peak[j]-20,peak[j]+20);                                                                                   
      f[i][j]->SetParameters(100,peak[j],5);                                                                                                              
    }                                                                                                                                                     
                                                                                                                                                          
    for(Int_t j=0;j<peak_n;j++){                                                                                                                              
      h[i]->Fit(Form("f%d_%d",i,j),"+","",peak[j]-20,peak[j]+20);                                                                                         
      mean[i][j]=f[i][j]->GetParameter(1);                                                                                                                
      sigma[i][j]=f[i][j]->GetParameter(2);                                                                                                               
      cout << mean[i][j] << " " << sigma[i][j] << endl;                                                                                                   
      ofs << i << " " << j << " " << mean[i][j] << " " << sigma[i][j] << endl;                                                                            
    }                                                                                                                                                     

    h[i]->Draw();
  }                                                                                                                                                       
                                                                                                                                                          
  ofs.close();                                                                                                                                            
                                                                                                                                                          
  return 0;                                                                                                                                               
}                                                                                                                                                         
                       
