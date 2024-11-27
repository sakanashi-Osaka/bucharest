void plot_hist_for_Dron() {
    // ファイルを開く
    std::ifstream infile("hist_for_Dron.txt");
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open data file." << std::endl;
        return;
    }

    // データを格納する変数
    std::vector<double> run, x, y;

    // ファイルからデータを読み込む
    double run_val, x_val, y_val;
    while (infile >> run_val >> x_val >> y_val) {
        run.push_back(run_val);
        x.push_back(x_val);
        y.push_back(y_val);
    }

    // データポイント数
    int n = run.size();

    // グラフを作成
    TGraph *graph = new TGraph(n, x.data(), y.data());
    graph->SetTitle("Graph from data.txt;X-axis;Y-axis"); // タイトルと軸ラベルを設定

    // キャンバスを作成して描画
    TCanvas *c1 = new TCanvas("c1", "Graph", 800, 600);
    graph->SetMarkerStyle(20); // マーカーのスタイル
    graph->SetMarkerSize(0.8); // マーカーのサイズ
    graph->SetMarkerColor(kBlue); // マーカーの色
    graph->Draw("AP"); // A: Axes, P: Points



    // 平均と分散を計算する関数
    auto calculate_mean_and_variance = [](const std::vector<double>& data, double& mean, double& variance) {
        mean = 0.0;
        variance = 0.0;
        int n = data.size();
        for (const auto& val : data) {
            mean += val;
        }
        mean /= n;
        for (const auto& val : data) {
            variance += (val - mean) * (val - mean);
        }
        variance /= n;
    };

    // x_val と y_val の平均と分散を計算
    double mean_x, variance_x, mean_y, variance_y;
    calculate_mean_and_variance(x, mean_x, variance_x);
    calculate_mean_and_variance(y, mean_y, variance_y);

    // 結果を表示
    std::cout << "X: Mean = " << mean_x << ", Variance = " << variance_x 
              << ", Sigma = " << std::sqrt(variance_x) << std::endl;
    std::cout << "Y: Mean = " << mean_y << ", Variance = " << variance_y 
              << ", Sigma = " << std::sqrt(variance_y) << std::endl;

}
