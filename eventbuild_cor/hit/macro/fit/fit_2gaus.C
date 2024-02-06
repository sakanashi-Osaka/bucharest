{
  TF1 *f2 = new TF1("f2","[0] + [1]*x + [2]*exp(-pow(x-[3],2)/(2*pow([4],2))) + [5]*exp(-pow(x-[6],2)/(2*pow([7],2)))",7,8);

  TH1D *h = new TH1D("h","h",200,6,8);
  tree->Draw("Ex4He>>h","chf4He<5 && chr4He<4 && chr4He>-1","");
  

  //for 13C
  double pi = h->GetBinContent(0);
  double pj = h->GetBinContent(200);
  double pmax = h->GetMaximum();
  double pmin = h->GetMinimum();
  double pex = 6.0 + (8.0-6.0)/200*h->GetMaximumBin();
  
  
  f2->SetParameter(0,pi);
  f2->SetParameter(1,0);
  
  f2->SetParameter(2,(pmax+pmin)/2);
  f2->SetParLimits(2,pmin,pmax);
  f2->SetParameter(3,pex-0.1);
  f2->SetParLimits(3,pex-0.1,pex);
  f2->SetParameter(4,0.1);
  f2->SetParLimits(4,0.05,0.15);
  
  f2->SetParameter(5,(pmax+pmin)/2);
  f2->SetParLimits(5,pmin,pmax);
  f2->SetParameter(6,pex+0.1);
  f2->SetParLimits(6,pex+0.01,pex+0.15);
  f2->SetParameter(7,0.1);
  f2->SetParLimits(7,0.05,0.15);
   
  
  h->Fit(f2,"","",6.8,7.6);

  double p0=f2->GetParameter(0);
  double p1=f2->GetParameter(1);
  double p2=f2->GetParameter(2);
  double p3=f2->GetParameter(3);
  double p4=f2->GetParameter(4);
  double p5=f2->GetParameter(5);
  double p6=f2->GetParameter(6);
  double p7=f2->GetParameter(7);

  double x[150]={};
  double y0[150]={};
  double y1[150]={};
  for(int i=0; i<150; i++){
    x[i]=6.0+(double)i*3/150;
    y0[i]=p2*exp(-pow(x[i]-p3,2)/(2*pow(p4,2)));
    y1[i]=p5*exp(-pow(x[i]-p6,2)/(2*pow(p7,2)));
  }
  
  
  TGraph *g0 = new TGraph(150,x,y0);
  TGraph *g1 = new TGraph(150,x,y1);
  g0->Draw("same");
  g1->Draw("same");

  cout << "yield: " << p2*p4*pow(3.1415,0.5)*150/3 << endl;
  cout << "yield: " << p5*p7*pow(3.1415,0.5)*150/3 << endl;
}
