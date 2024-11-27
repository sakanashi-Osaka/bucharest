{
  //  TF1 *f2 = new TF1("f2","[0] + [1]*x + [2]*(1858.18*exp(-pow(x-7.12369,2)/(2*pow(0.0712369,2))) + 1277.87*exp(-pow(x-7.22852,2)/(2*pow(0.111216,2))))",6.9,7.4);

  //  TF1 *f2 = new TF1("f2","[0] + [1]*x + [2]*(714.865*exp(-pow(x-7.19236,2)/(2*pow(0.0736675,2))) + 278.834*exp(-pow(x-7.32913,2)/(2*pow(0.0963947,2))))",6.9,7.4); // seg0 raw
  TF1 *f2 = new TF1("f2","[0] + [1]*x + [2]*(714.865*exp(-pow(x-7.19236+[3],2)/(2*pow(0.0736675,2))) + 278.834*exp(-pow(x-7.32913+[3],2)/(2*pow(0.0963947,2))))",6.95,7.35); //seg0 sift

  TH1D *h = new TH1D("h","h",200,6,8);
  tree->Draw("Ex4He>>h","chf4He<5 && chr4He<4 && chr4He>-1","");  

  //for 13C
  f2->SetParameter(0,300);
  f2->SetParameter(1,0);
  
  f2->SetParameter(2,0.2);
  f2->SetParLimits(2,0.05,1);

  f2->SetParameter(3,0);
  f2->SetParLimits(3,-0.05,0.05);
   
  
  h->Fit(f2,"","",6.95,7.35);

  /*
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
  */
}
