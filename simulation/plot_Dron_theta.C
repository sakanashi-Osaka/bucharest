void plot_Dron_theta() {
 
    gStyle->SetOptStat(0);
    TCanvas *c1 = new TCanvas("c1", "Canvas 1", 600, 600);

    // 最初のROOTファイルを開く                                                                                                            
    TFile *file1 = TFile::Open("rootfile/Dron.root");

    // Treeを取得                                                                                                                          
    TTree *tree1 = (TTree*)file1->Get("tree");

    TH2D *h = new TH2D("h","h",50,28,52,50,28,52);
    TH2D *hh = new TH2D("hh","hh",50,28,52,50,28,52);
    TH2D *hhh = new TH2D("hhh","hhh",50,28,52,50,28,52);

    // ヒストグラムhとhhを作成して散布図として表示                                                                                         
    tree1->Draw("mth3a[]:mth4He>>h", "m12CchF>-1&&m12CchR>-1&&m4HechR>-1&&m4HechF>-1&&Entry$%50==0", "");
    tree1->Draw("mth12C:mth4He>>hh", "m12CchF>-1&&m12CchR>-1&&m4HechR>-1&&m4HechF>-1&&Entry$%20==0", "same");

    h->SetMarkerStyle(8);
    h->SetMarkerColor(1);
    h->SetMarkerSize(0.2);

    hh->SetMarkerStyle(8);
    hh->SetMarkerColor(2);
    hh->SetMarkerSize(0.4);

    TF1 *line1 = new TF1("line1", "52 - 0.2*x", 28, 52);
    line1->SetLineColor(kRed);
    line1->SetLineWidth(2);
    line1->Draw("same");
    TF1 *line2 = new TF1("line1", "68 - 0.4*x", 28, 52);
    line2->SetLineColor(kRed);
    line2->SetLineWidth(2);
    line2->Draw("same");

    // キャンバス2を作成                                                                                                                   
    TCanvas *c2 = new TCanvas("c2", "Canvas 2", 600, 600);

    // 次のROOTファイルを開く                                                                                                              
    TFile *file2 = TFile::Open("rootfile/Dron12C_all.root");

    // Treeを取得                                                                                                                          
    TTree *tree2 = (TTree*)file2->Get("tree");

    hhh->SetMarkerStyle(8);
    hhh->SetMarkerColor(2);
    hhh->SetMarkerSize(0.4);

    // 新たなヒストグラムを作成                                                                                                         
    tree2->Draw("th12C:th4He>>hhh", "Amax4He>-0.015&&Amax4He<0.02&&fabs(Ex4He-7.65)<0.25&&chf4He<14&&fabs(Amax12C)<1.2", "col", 2e7);

    // y = 9.3 - x/2 の直線を追加
    line1->Draw("same");
    line2->Draw("same");
}
