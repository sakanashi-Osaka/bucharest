{
  TFile f12C("../rootfile/gomi12C_68.root");
  TFile f13C("../rootfile/gomi13C_11.root");
  TFile fb("../rootfile/gomib_11.root");

  TTree* tree12C = (TTree*)f12C.Get("tree");
  TTree* tree13C = (TTree*)f13C.Get("tree");
  TTree* treeb = (TTree*)fb.Get("tree");

  TH1D *h12C_Ex = new TH1D("h12C_Ex","h12C_Ex",500,-1,9);
  TH1D *h13C_Ex = new TH1D("h13C_Ex","h13C_Ex",500,-1,9);
  TH1D *hb_Ex = new TH1D("hb_Ex","hb_Ex",500,-1,9);
  
  tree12C->Draw("Ex4He>>h12C_Ex","","");
  tree13C->Draw("Ex4He>>h13C_Ex","","");
  treeb->Draw("Ex4He>>hb_Ex","","");

  TFile *fout12C =new TFile("hist12C.root","recreate");
  h12C_Ex->Write();
  TFile *fout13C =new TFile("hist13C.root","recreate");
  h13C_Ex->Write();
  TFile *foutb =new TFile("histb.root","recreate");
  hb_Ex->Write();
}
