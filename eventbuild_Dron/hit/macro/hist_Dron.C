void hist_Dron() {
    // ROOTファイルを開く
  TFile *file1 = TFile::Open("rootfile/Dron_2236_2239.root");
  //  TFile *file1 = TFile::Open("rootfile/Dron_2243_2247.root");
  //  TFile *file1 = TFile::Open("rootfile/Dron_2250_2269.root");
  //  TFile *file1 = TFile::Open("rootfile/Dron_2287_2294.root");
  //  TFile *file1 = TFile::Open("rootfile/Dron_2297_2306.root");
  //  TFile *file1 = TFile::Open("rootfile/Dron_2310_2318.root");
  //  TFile *file1 = TFile::Open("rootfile/Dron_2322_2327.root");
  //  TFile *file1 = TFile::Open("rootfile/Dron_2330_2337.root");
  //  TFile *file1 = TFile::Open("rootfile/Dron_2340_2345.root");
  //  TFile *file1 = TFile::Open("rootfile/Dron_2349_2354.root");
  //  TFile *file1 = TFile::Open("rootfile/Dron_2357_2362.root");
  //  TFile *file1 = TFile::Open("rootfile/Dron_2371_2376.root");
  //  TFile *file1 = TFile::Open("rootfile/Dron_2379_2386.root");
  //  TFile *file1 = TFile::Open("rootfile/Dron_2389_2394.root");
  //  TFile *file1 = TFile::Open("rootfile/Dron_2404_2411.root");


  
  TFile *file2 = TFile::Open("rootfile/Dron13C_all.root");
  TFile *file3 = TFile::Open("rootfile/Dronb_all.root");

    if (!file1 || !file2 || !file3) {
        std::cerr << "Error: One or more files could not be opened!" << std::endl;
        return;
    }

    // ツリーを取得
    TTree *tree1 = (TTree*)file1->Get("tree");
    TTree *tree2 = (TTree*)file2->Get("tree");
    TTree *tree3 = (TTree*)file3->Get("tree");

    if (!tree1 || !tree2 || !tree3) {
        std::cerr << "Error: One or more trees could not be found!" << std::endl;
        return;
    }

    // ヒストグラムを生成
    TH1F *h12C = new TH1F("h12C", "h12C", 200, 6.5, 8.5);
    TH1F *h13C = new TH1F("h13C", "h13C", 200, 6.5, 8.5);
    TH1F *hb = new TH1F("hb", "hb", 200, 6.5, 8.5);

    // ヒストグラムにデータを描画
    tree1->Draw("Ex4He>>h12C", "abs(Amax4He+0.0025)<0.0175 && chf4He<14", "goff");
    tree2->Draw("Ex4He>>h13C", "abs(Amax4He+0.0025)<0.0175 && chf4He<14", "goff");
    //    tree2->Draw("Ex4He>>h13C", "abs(Amax4He+0.0025)<0.0175 && chf4He<14", "goff");
    //    tree2->Draw("Ex4He+0.01>>h13C", "abs(Amax4He+0.0025)<0.0175 && chf4He<14", "goff");
    //    tree2->Draw("Ex4He+0.01>>h13C", "abs(Amax4He+0.0025)<0.0175 && chf4He<14", "goff");
    //  tree2->Draw("Ex4He+0.01>>h13C", "abs(Amax4He+0.0025)<0.0175 && chf4He<14", "goff");
    //    tree2->Draw("Ex4He+0.05>>h13C", "abs(Amax4He+0.0025)<0.0175 && chf4He<14", "goff");
    //    tree2->Draw("Ex4He>>h13C", "abs(Amax4He+0.0025)<0.0175 && chf4He<14", "goff");
    //    tree2->Draw("Ex4He>>h13C", "abs(Amax4He+0.0025)<0.0175 && chf4He<14", "goff");
    //  tree2->Draw("Ex4He+0.01>>h13C", "abs(Amax4He+0.0025)<0.0175 && chf4He<14", "goff");
    //    tree2->Draw("Ex4He+0.01>>h13C", "abs(Amax4He+0.0025)<0.0175 && chf4He<14", "goff");
    //    tree2->Draw("Ex4He>>h13C", "abs(Amax4He+0.0025)<0.0175 && chf4He<14", "goff");
    //    tree2->Draw("Ex4He-0.01>>h13C", "abs(Amax4He+0.0025)<0.0175 && chf4He<14", "goff");
    //    tree2->Draw("Ex4He-0.005>>h13C", "abs(Amax4He+0.0025)<0.0175 && chf4He<14", "goff");
    //    tree2->Draw("Ex4He>>h13C", "abs(Amax4He+0.0025)<0.0175 && chf4He<14", "goff");
    //    tree2->Draw("Ex4He-0.01>>h13C", "abs(Amax4He+0.0025)<0.0175 && chf4He<14", "goff");
    tree3->Draw("Ex4He>>hb", "abs(Amax4He)<0.015 && chf4He<14", "goff");

    // h13C と hb をキャプチャしたスケーリング関数にガウス項を追加
    auto scalingFunction = [h13C, hb](double *x, double *par) {
        double xValue = x[0];
        int bin13C = h13C->FindBin(xValue);
        int binb = hb->FindBin(xValue);

        double content13C = h13C->GetBinContent(bin13C);
        double contentb = hb->GetBinContent(binb);

        // ガウス項
        double gaussian = par[2] * exp(-0.5 * pow((xValue - par[3]) / par[4], 2));
        // ガウス項2
        double gaussian_high = par[5] * exp(-0.5 * pow((xValue - par[6]) / par[7], 2));

        return par[0] * content13C + par[1] * contentb + gaussian + gaussian_high;
    };

    // フィット関数を作成
    TF1 *fitFunc = new TF1("fitFunc", scalingFunction, 6.5, 8.5, 8);
    //    TF1 *fitFunc = new TF1("fitFunc", scalingFunction, 6.7, 8.2, 8);
    fitFunc->SetParNames("Scale_h13C", "Scale_hb", "Gauss_Amplitude", "Gauss_Mean", "Gauss_Sigma", "Gauss2_Amplitude", "Gauss2_Mean", "Gauss2_Sigma");

    // パラメータの初期値と範囲を設定
    //1
    fitFunc->SetParameters(0.00825, 0.2326, 1000, 7.6, 0.1, 100, 9, 1);
    fitFunc->SetParLimits(0, 0.00825, 0.00825);  // Scale_h13C
    fitFunc->SetParLimits(1, 0.2326, 0.2326);    // Scale_hb
    //2
    //    fitFunc->SetParameters(0.00681, 0.192, 1000, 7.6, 0.1, 100, 9, 1);
    //    fitFunc->SetParLimits(0, 0.00681, 0.00681);  // Scale_h13C
    //    fitFunc->SetParLimits(1, 0.192, 0.192);    // Scale_hb
    //3
    //    fitFunc->SetParameters(0.01157, 0.326, 1000, 7.6, 0.1, 100, 9, 1);
    //    fitFunc->SetParLimits(0, 0.01157, 0.01157);  // Scale_h13C
    //    fitFunc->SetParLimits(1, 0.326, 0.326);    // Scale_hb
    //4
    //    fitFunc->SetParameters(0.01379, 0.388, 1000, 7.6, 0.1, 100, 9, 1);
    //    fitFunc->SetParLimits(0, 0.01379, 0.01379);  // Scale_h13C
    //    fitFunc->SetParLimits(1, 0.388, 0.388);    // Scale_hb
    //5
    //    fitFunc->SetParameters(0.0167, 0.47, 1000, 7.6, 0.1, 100, 9, 1);
    //    fitFunc->SetParLimits(0, 0.0167, 0.0167);  // Scale_h13C
    //    fitFunc->SetParLimits(1, 0.47, 0.47);    // Scale_hb
    //6
    //    fitFunc->SetParameters(0.012, 0.34, 1000, 7.6, 0.1, 100, 9, 1);
    //    fitFunc->SetParLimits(0, 0.012, 0.012);  // Scale_h13C
    //    fitFunc->SetParLimits(1, 0.34, 0.34);    // Scale_hb

    
    //    fitFunc->SetParameters(0.01, 0.2, 1000, 7.6, 0.1, 100, 9, 1);
    //    fitFunc->SetParLimits(0, 0.001, 0.2);  // Scale_h13C
    //    fitFunc->SetParLimits(1, 0.001, 10.0);    // Scale_hb

    fitFunc->SetParLimits(2, 10, 1000000);  // Gauss_Amplitude
    fitFunc->SetParLimits(3, 7.55, 7.75);  // Gauss_Mean
    fitFunc->SetParLimits(4, 0.01, 0.1);   // Gauss_Sigma
    fitFunc->SetParLimits(5, 1, 10000);    // Gauss2_Amplitude
    fitFunc->SetParLimits(6, 8.0, 10.0);    // Gauss2_Mean
    fitFunc->SetParLimits(7, 0.1, 2);   // Gauss2_Sigma

    // フィットを実行
    h12C->Fit(fitFunc, "R");

    
      
        // 各項の関数を個別に定義
    // h13C のスケールした内容
    TH1F *scaled_h13C = (TH1F*)h13C->Clone("scaled_h13C");
    scaled_h13C->Scale(fitFunc->GetParameter(0));

    // hb のスケールした内容
    TH1F *scaled_hb = (TH1F*)hb->Clone("scaled_hb");
    scaled_hb->Scale(fitFunc->GetParameter(1));

    // ガウス項1の関数
    TF1 *fGaussian1 = new TF1("fGaussian1", "[0]*exp(-0.5*((x-[1])/[2])^2)", 6.5, 8.5);
    fGaussian1->SetParameters(fitFunc->GetParameter(2), fitFunc->GetParameter(3), fitFunc->GetParameter(4));

    // ガウス項2の関数
    TF1 *fGaussian2 = new TF1("fGaussian2", "[0]*exp(-0.5*((x-[1])/[2])^2)", 6.5, 8.5);
    fGaussian2->SetParameters(fitFunc->GetParameter(5), fitFunc->GetParameter(6), fitFunc->GetParameter(7));

    // 結果を描画
    TCanvas *c2 = new TCanvas("c2", "Fit Results", 800, 600);
    h12C->SetLineColor(1);
    h12C->GetYaxis()->SetRangeUser(0, 5000);
    h12C->Draw();
    
    fitFunc->SetLineColor(kRed);
    fitFunc->Draw("same");

    scaled_h13C->SetLineColor(kBlue);
    scaled_h13C->Draw("same hist");

    scaled_hb->SetLineColor(kGreen);
    scaled_hb->Draw("same hist");

    fGaussian1->SetLineColor(kMagenta);
    fGaussian1->Draw("same");

    fGaussian2->SetLineColor(kOrange);
    fGaussian2->Draw("same");


    double areaGaussian1 = sqrt(2 * M_PI) * fitFunc->GetParameter(2) * fitFunc->GetParameter(4);
    std::cout << "Gaussian 1 Area: " << areaGaussian1 / (2./200.) << std::endl;

    double chi2 = fitFunc->GetChisquare(); // カイ二乗
    double ndf = fitFunc->GetNDF();        // 自由度
    double chi2_ndf = chi2 / ndf;          // カイ二乗/自由度

    // 結果を出力
    std::cout << "Chi2 = " << chi2 << std::endl;
    std::cout << "NDF = " << ndf << std::endl;
    std::cout << "Chi2/NDF = " << chi2_ndf << std::endl;
}
