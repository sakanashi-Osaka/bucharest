{
  double factor_13C = 0.20;
  double Eshift_13C = -0.0;
  double Eshift_12C = 0.;
  
  TFile *f0 =new TFile("../rootfile/test12C_68.root");
  TTree *tree0 = (TTree*)f0->Get("tree");
  TFile *f1 =new TFile("../rootfile/test13C_11.root");
  TTree *tree1 = (TTree*)f1->Get("tree");
  
  TH1D *hadd = new TH1D("hadd","hadd",100,6,8);
  TH1D *h = new TH1D("h","h",100,6,8);
  TH1D *hexp = new TH1D("hexp","hexp",100,6,8);

  const char *cut_name = "C_ts&&C_ph&&C_th&&K12C>9.3-K4He*0.5&&K12C<16.8-K4He&&K12C>2.7&&Amax12C<0.018";
  const char *cut_name2 = "chr4He!=3&&chr4He!=4&&chr4He!=11&&chr4He!=12";
  tree0->Draw(Form("Ex4He+%f>>h12C(100,6,8)",Eshift_12C),Form("%s&&%s",cut_name,cut_name2),"");
  tree1->Draw(Form("Ex4He+%f>>h13C(100,6,8)",Eshift_13C),Form("%s&&%s",cut_name,cut_name2),"");

  //  /*
  h12C->SetLineColor(2);
  h12C->Draw();
  h->Add(h13C,factor_13C);
  h->SetLineColor(1);
  h->Draw("same");
  //  */
  /*  
  hadd->Add(h12C,1);
  hadd->Add(h13C,-1*factor_13C);
  hadd->Draw();

  TF1 *fit = new TF1("fit","[0]*exp(-pow(x-[1],2)/(2*pow([2],2))) + [3]*exp(-pow(x-[4],2)/(2*pow([5],2))) + [6]*(x-[7])",7.2,8.2);
  fit->SetParameter(0,2000);
  fit->SetParameter(1,9);
  fit->SetParameter(2,0.2);

  fit->SetParameter(3,100);
  fit->SetParLimits(3,50,200);
  fit->SetParameter(4,7.6);
  fit->SetParameter(5,0.08);

  fit->SetParameter(6,10);
  fit->SetParameter(7,7);
  hadd->Fit("fit","","",7.2,8.2);


  double p0, p1, p2, p3, p4, p5, p6, p7;
  p0 = fit->GetParameter(0);
  p1 = fit->GetParameter(1);
  p2 = fit->GetParameter(2);
  p6 = fit->GetParameter(6);
  p7 = fit->GetParameter(7);
  */
  /*
  double tmp[100];
  for(int i=0; i<100; i++){
    double x = (8.2-6.2)/100*i + 6.2;
    tmp[i] = p0*exp(-pow(x-p1,2)/(2*pow(p2,2)))+p6*(x-p7);
    hexp->SetBinContent(i+1,tmp[i]);
  }
  hadd->Add(hexp,-1);
  hadd->Draw("");
  */
}
