{
  TCanvas *c = new TCanvas("c","c",700,700);
  bool flag_bg = 1;
  bool flag_12C = 1;

  double pb0, pb1, pb2, pb3, pb4, pb5;
  double p12C0, p12C1, p12C2, p12C3, p12C4, p12C5;

  TFile fb("histb.root");    
  TF1 *fib = new TF1("fib","exp([0]+[1]*x+[2]*x*x) + [3]*exp(-pow(x-[4],2)/(2*pow([5],2)))",1,5);
  if(flag_bg){
    //*** bg function ***// 
    fib->SetParameter(0,4);
    fib->SetParameter(1,-1);
    fib->SetParameter(2,0);
    fib->SetParameter(3,100);
    fib->SetParameter(4,1.42);
    fib->SetParameter(5,0.1);

    hb_Kg->Fit("fib","","",1,5);
    
    pb0=fib->GetParameter(0);
    pb1=fib->GetParameter(1);
    pb2=fib->GetParameter(2);
    pb3=fib->GetParameter(3);
    pb4=fib->GetParameter(4);
    pb5=fib->GetParameter(5);
  }

  TFile f12C("hist12C.root");  
  TF1 *fi12C = new TF1("fi12C","(exp([0]+[1]*x+[2]*x*x) + [3]*exp(-pow(x-[4],2)/(2*pow([5],2)))) *[6]",1,5);
  if(flag_12C){
    //*** 12C+bg function ***// 
    fi12C->SetParameter(0,pb0);
    fi12C->SetParameter(1,pb1);
    fi12C->SetParameter(2,pb2);
    fi12C->SetParameter(3,pb3);
    fi12C->SetParameter(4,pb4);
    fi12C->SetParameter(5,pb5);
    fi12C->SetParLimits(0,pb0,pb0);
    fi12C->SetParLimits(1,pb1,pb1);
    fi12C->SetParLimits(2,pb2,pb2);
    fi12C->SetParLimits(3,pb3,pb3);
    fi12C->SetParLimits(4,pb4,pb4);
    fi12C->SetParLimits(5,pb5,pb5);
    //factor of bg
    fi12C->SetParameter(6,2.35);
    fi12C->SetParLimits(6,2.35,2.35);
    
    h12C_Kg->Fit("fi12C","","",1,5);
    
    p12C0=fi12C->GetParameter(0);
    p12C1=fi12C->GetParameter(1);
    p12C2=fi12C->GetParameter(2);
  }
  h12C_Kg->SetLineColor(1);
  h12C_Kg->Draw("");

  
  TFile f13C("hist13C.root");
  TH1D *h = new TH1D("h","h",200,1,5);                                                        
  h->Add(h13C_Kg,0.18);
  h->Add(fi12C,1.2);
  h->SetLineColor(2);
  h->Draw("same")
  
}
