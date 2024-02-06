{
  TCanvas *c0 = new TCanvas("c0","c0",600,600);
  c0->Divide(2,2);
  c0->cd(1);
  int run = 2317;
  TFile *fin =new TFile(Form("../../eventbuild_cor/hit/rootfile/run%d.root",run));
  TTree *hit = (TTree*)fin->Get("hit");
  hit->Draw("ch_f[][0]*16+ch_r[][0]:ts_diff_f[][0]/1e3>>h(100,-100,100,256,0,256)","","col",1e6);
}
