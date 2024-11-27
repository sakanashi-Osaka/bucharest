{
  TFile f12C("../rootfile/gomi12C_68.root");
  TFile f13C("../rootfile/gomi13C_11.root");
  TFile fb("../rootfile/gomib_11.root");

  TTree* tree12C = (TTree*)f12C.Get("tree");
  TTree* tree13C = (TTree*)f13C.Get("tree");
  TTree* treeb = (TTree*)fb.Get("tree");

  TH1D *h12C_Kg = new TH1D("h12C_Kg","h12C_Kg",200,1,5);
  TH1D *h13C_Kg = new TH1D("h13C_Kg","h13C_Kg",200,1,5);
  TH1D *hb_Kg = new TH1D("hb_Kg","hb_Kg",200,1,5);
  
  tree12C->Draw("Kg[0]>>h12C_Kg","abs(Ex4He-3.5)<0.35","");
  tree13C->Draw("Kg[0]>>h13C_Kg","abs(Ex4He-3.5)<0.35","");
  treeb->Draw("Kg[0]>>hb_Kg","abs(Ex4He-3.5)<0.35","");

  TFile *fout12C =new TFile("hist12C.root","recreate");
  h12C_Kg->Write();
  TFile *fout13C =new TFile("hist13C.root","recreate");
  h13C_Kg->Write();
  TFile *foutb =new TFile("histb.root","recreate");
  hb_Kg->Write();
}
