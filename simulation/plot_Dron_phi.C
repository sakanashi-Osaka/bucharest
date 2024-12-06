void plot_Dron_phi() {
 
    gStyle->SetOptStat(0);
    TCanvas *c1 = new TCanvas("c1", "Canvas 1", 600, 600);

    // 最初のROOTファイルを開く                                                                                                            
    TFile *file1 = TFile::Open("rootfile/Dron.root");

    // Treeを取得                                                                                                                          
    TTree *tree1 = (TTree*)file1->Get("tree");

    TH1D *h = new TH1D("h","h",100,130,230);
    TH1D *hh = new TH1D("hh","hh",100,130,230);

    h->SetLineColor(1);
    hh->SetLineColor(2);
    hh->SetLineWidth(2);

    // ヒストグラムhとhhを作成して散布図として表示                                                                                         
    tree1->Draw("(mph3a[]-mph4He+360)%360>>h", "m12CchF>-1&&m12CchR>-1&&m4HechR>-1&&m4HechF>-1&&Entry$%1==0", "");
    tree1->Draw("(mph12C-mph4He+360)%360>>hh", "m12CchF>-1&&m12CchR>-1&&m4HechR>-1&&m4HechF>-1&&Entry$%1==0", "same");



}
