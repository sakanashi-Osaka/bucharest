#include <fstream>
#include <iostream>

void calib() {
    // 出力ファイルを開く
    std::ofstream outputFile("ex_cor.prm");
    if (!outputFile.is_open()) {
        std::cerr << "Error opening output file!" << std::endl;
        return;
    }

    // ループでrun番号を処理
    for (int run = 2236; run <= 2413; ++run) {
      cout << run << endl;
      // ファイル名を作成
        TString fileName = TString::Format("../rootfile/Dron%04d.root", run);

        // ファイルを開く
        TFile *file = TFile::Open(fileName);
        if (!file || file->IsZombie()) {
            std::cerr << "Error opening file: " << fileName << std::endl;
            continue;
        }

        // ツリーを取得
        TTree *tree = (TTree*)file->Get("tree");
        if (!tree) {
            std::cerr << "Tree not found in file: " << fileName << std::endl;
            file->Close();
            continue;
        }

        // 1次元ヒストグラムを作成
        TH1F *h = new TH1F("h", "Ex4He Distribution", 100, 7, 8);
        tree->Draw("Ex4He>>h", "Amax4He<0.015&&Amax4He>-0.02", "col");

        // 最大値のbinを取得
        Int_t maxBin = h->GetMaximumBin();
        double centerX = h->GetXaxis()->GetBinCenter(maxBin);

        // フィット
        TF1 *gausFit = new TF1("gausFit", "gaus", centerX - 0.1, centerX + 0.1);
        h->Fit(gausFit, "Q"); // 静かにフィット

        // フィットパラメータを取得
        double mean = gausFit->GetParameter(1);
        double sigma = gausFit->GetParameter(2);

        // 結果を出力
        outputFile << run << " " << mean << " " << sigma << std::endl;

        // メモリ解放
        delete h;
        delete gausFit;
        file->Close();
        delete file;
    }

    // 出力ファイルを閉じる
    outputFile.close();
    std::cout << "Analysis completed. Results saved to ex_cor.prm." << std::endl;
}
