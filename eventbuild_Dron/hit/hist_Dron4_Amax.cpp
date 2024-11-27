{
  TCanvas *c0 = new TCanvas("c0","c0",700,600);
  gStyle->SetOptStat(0);

  TFile *f1 = TFile::Open("rootfile/Dron2236.root");
  hit->Draw("Amax_origin[1][0]:Energy_f[1][0]>>hh(400,1,25,400,0.3,0.7)","hit_n[1]==11&&ch_f[1][0]==8&&ch_r[1][0]==5","",1e2);
  hh->SetMarkerStyle(8);
  hh->SetMarkerSize(0.3);
  hh->SetMarkerColor(1);
  hh->Draw("");

  TFile *f2 = TFile::Open("rootfile/Dron2287.root");
  hit->Draw("Amax_origin[1][0]:Energy_f[1][0]*1.00>>hhh(400,1,25,400,0.3,0.7)","hit_n[1]==11&&ch_f[1][0]==8&&ch_r[1][0]==5","same",2e6);
  hhh->SetMarkerStyle(8);
  hhh->SetMarkerSize(0.3);
  hhh->SetMarkerColor(1);
  hhh->Draw("same");

  TFile *f0 = TFile::Open("rootfile/Dron2305.root");
  hit->Draw("Amax_origin[1][0]:Energy_f[1][0]>>h(400,1,25,400,0.3,0.7)","hit_n[1]==11&&ch_f[1][0]==8&&ch_r[1][0]==5","same",2e6);
  h->SetMarkerStyle(8);
  h->SetMarkerSize(0.3);
  h->SetMarkerColor(2);
  h->Draw("same");

}
