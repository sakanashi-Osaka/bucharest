int circle(int ev){
  double Max = TMath::Pi()*2;
  double theta0[500]={};
  double x0[1000]={};
  double y0[1000]={};
  double theta1[16]={};
  double x1[1600]={};
  double y1[1600]={};

  for (int i=0; i<500; i++) {
    theta0[i] = Max/500.*i;
    x0[i] =24*cos(theta0[i]);
    y0[i] =24*sin(theta0[i]);
    x0[i+500] =48*cos(theta0[i]);
    y0[i+500] =48*sin(theta0[i]);
  }

  for (int i=0; i<16; i++) {
    theta1[i] = Max/16.*i;
    for (int j=0; j<100; j++) {
      x1[i*100+j]=48*std::cos(theta1[i])*j/100;
      y1[i*100+j]=48*std::sin(theta1[i])*j/100;
    }
  }

  //  TEventList *gate = new TEventList("gate","gate");
  //  tree->Draw(">>gate","idp4He==1&&m4HechR<4","col");
  //  tree->Draw("s3apos[][1]:s3apos[][0]>>h(50,-100,50,50,-100,50)","idp4He==1&&m4HechR<4","colz",1,gate->GetEntry(ev));

  //    TGraph *g0 = new TGraph(1000, x0, y0);
  //    TGraph *g1 = new TGraph(1600, x1, y1);
  //    g0->Draw("same");
  //    g1->Draw("same");

    tree->Scan("s3apos[][1]:s3apos[][0]","idp4He==1&&m4HechR<4","",1,gate->GetEntry(ev));
  
  return 0;
}

