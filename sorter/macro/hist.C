{
  TFile *f0 =new TFile("../run2230_-1_ssgant1.root");
  TTree *tree0 = (TTree*)f0->Get("tree");  
  TFile *f1 =new TFile("../run2269_-1_ssgant1.root");
  TTree *tree1 = (TTree*)f1->Get("tree");  
  TFile *f2 =new TFile("../run2287_-1_ssgant1.root");
  TTree *tree2 = (TTree*)f2->Get("tree");  

  tree0->Draw("Energy>>h0(400,9,21)","domain==86","",5e7,0);
  tree1->Draw("Energy>>h1(400,9,21)","domain==86","",5e7,0);
  tree2->Draw("Energy>>h2(400,9,21)","domain==86","",5e7,0);

  h0->SetLineColor(1);
  h1->SetLineColor(4);
  h2->SetLineColor(2);

  h0->Draw();
  h1->Draw("same");
  h2->Draw("same");

}
 
