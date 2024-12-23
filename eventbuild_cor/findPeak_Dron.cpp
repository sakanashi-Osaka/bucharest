#include <vector>

void findPeaks(int run, float &peak1, float &peak2) {
    // ファイル名を構築
    TString filename = TString::Format("rootfile/Dron%d.root", run);
    
    // ファイルを開く
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file '" << filename << "'" << std::endl;
        peak1 = -1; // エラー値を設定
        peak2 = -1; // エラー値を設定
        return;
    }
    cout << "Analyzing run" << run << endl;
    
    // Treeを取得
    TTree *tree = (TTree*)file->Get("event");
    if (!tree) {
        std::cerr << "Error: Cannot find tree in the file." << std::endl;
        file->Close();
        peak1 = -1; // エラー値を設定
        peak2 = -1; // エラー値を設定
        return;
    }

    // ブランチをリンク
    float Energy[10][16];
    //    tree->SetBranchAddress("Energy", Energy);
    tree->SetBranchAddress("Energy2423", Energy);

    // ヒストグラムを作成
    TH1D *hist1 = new TH1D("hist1", "Energy[5][10] (16 < E < 18)", 100, 16, 18);
    TH1D *hist2 = new TH1D("hist2", "Energy[5][10] (20 < E < 22)", 100, 20, 22);

    // ヒストグラムにデータをフィル
    Long64_t nEntries = tree->GetEntries();
    //    for (Long64_t i = 0; i < nEntries; ++i) {
    for (Long64_t i = 0; i < 1e6; ++i) {
      tree->GetEntry(i);
        float value = Energy[5][10];
        if (16 < value && value < 18) {
            hist1->Fill(value);
        } else if (20 < value && value < 22) {
            hist2->Fill(value);
        }
    }
    // ピーク位置を探す
    TF1 *fit1 = new TF1("fit1", "gaus", 16, 18);
    hist1->Fit(fit1, "Q");
    peak1 = fit1->GetParameter(1);

    TF1 *fit2 = new TF1("fit2", "gaus", 20, 22);
    hist2->Fit(fit2, "Q");
    peak2 = fit2->GetParameter(1);

    // 後始末
    file->Close();
    //    delete hist1;
    //    delete hist2;
    //    delete fit1;
    //    delete fit2;
}

void processRuns(int startRun, int endRun) {
    // 配列を定義
    std::vector<int> runNumbers;
    std::vector<float> peak1Values;
    std::vector<float> peak2Values;

    // 指定範囲でfindPeaksを実行
    for (int run = startRun; run <= endRun; ++run) {
        float peak1, peak2;
        findPeaks(run, peak1, peak2);
        if (peak1 >= 0 && peak2 >= 0) { // 有効な結果のみ保存
            runNumbers.push_back(run);
            peak1Values.push_back(peak1);
            peak2Values.push_back(peak2);
        }
    }

    // グラフを作成
    TGraph *graph1 = new TGraph(runNumbers.size());
    TGraph *graph2 = new TGraph(runNumbers.size());
    
    for (size_t i = 0; i < runNumbers.size(); ++i) {
      graph1->SetPoint(i, runNumbers[i], peak1Values[i]);
      graph2->SetPoint(i, runNumbers[i], peak2Values[i]);
    }
    
    graph1->SetTitle("Peak1 vs RunNumbers;Run Number;Peak1 Value");
    graph2->SetTitle("Peak2 vs RunNumbers;Run Number;Peak2 Value");
    
    // 新しいROOTファイルに保存
    TFile *outputFile = new TFile("Peaks_energy2423_5_10_Dron.root", "RECREATE");
    graph1->Write("Peak1");
    graph2->Write("Peak2");
    outputFile->Close();


    
    // 結果を出力
    std::cout << "Run\tPeak1\tPeak2" << std::endl;
    for (size_t i = 0; i < runNumbers.size(); ++i) {
        std::cout << runNumbers[i] << "\t" << peak1Values[i] << "\t" << peak2Values[i] << std::endl;
    }
}
