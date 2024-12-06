#include <cmath>
#include <vector>

const double s1_r1 = 48.0 / 2.0;  // 内側の半径
const double s1_r2 = 96.0 / 2.0;  // 外側の半径
const double s1_dist = 41.0;      // 検出器までの距離
const double beam_x = 0.0, beam_y = 2.4;
const double PI = TMath::Pi();


// 先ほど計算した定数
const double n_effective = 5.683677395054122e18;  // 標的の有効原子密度 (atoms/cm^2)
const double N_beam = 2.069128914783562e15;      // ビーム粒子数

double get_front_theta(double front_ch) {
    double strp_wid = (s1_r2 - s1_r1) / 16.0;
    double r = s1_r1 + strp_wid * (0.5) + strp_wid * (15 - front_ch);

    return atan(r / s1_dist);
}

std::pair<double,double> get_angle(double ch_f, double ch_r){

  double tmp_x, tmp_y, tmp_z, tmp_theta, tmp_phi;
  double front_theta = get_front_theta(ch_f);
  double rear_phi = (90 + 22.5/2.0 + 22.5*ch_r)*PI/180;

  tmp_x = s1_dist * tan(front_theta)*cos(rear_phi) - beam_x;
  tmp_y = s1_dist * tan(front_theta)*sin(rear_phi) - beam_y;
  tmp_z = s1_dist;
  tmp_theta = acos(tmp_z/pow(pow(tmp_x,2)+pow(tmp_y,2)+pow(tmp_z,2),0.5));
  if(tmp_y>=0) tmp_phi = acos(tmp_x/pow(pow(tmp_x,2)+pow(tmp_y,2),0.5));
  if(tmp_y<0) tmp_phi = -1 * acos(tmp_x/pow(pow(tmp_x,2)+pow(tmp_y,2),0.5));

  return std::make_pair(tmp_theta,tmp_phi);
}

double get_distance(double ch_f, double ch_r) {

    double front_theta = get_front_theta(ch_f);
    double rear_phi = (90.0 + 22.5 / 2.0 + 22.5 * ch_r) * PI / 180.0;

    double tmp_x = s1_dist * tan(front_theta) * cos(rear_phi) - beam_x;
    double tmp_y = s1_dist * tan(front_theta) * sin(rear_phi) - beam_y;
    double tmp_z = s1_dist;

    double tmp_d = sqrt(tmp_x * tmp_x + tmp_y * tmp_y + tmp_z * tmp_z);
    return tmp_d;
}

double get_area(double ch_f, double ch_r) {
    const double PI = TMath::Pi();
    double tmp_r_outer = 48.0 - 1.5 * ch_f;
    double tmp_r_inner = 46.5 - 1.5 * ch_f;
    double tmp_area = PI * (pow(tmp_r_outer, 2) - pow(tmp_r_inner, 2)) / 16.0;

    return tmp_area;
}

void hist_Dron_dist() {
    // ROOTファイルを開く
  //    TFile *file = new TFile("rootfile/Dron_2297_2306.root");
    TFile *file = new TFile("rootfile/Dron12C_all.root");
    TTree *tree = (TTree*)file->Get("tree");

    const int n_bins = 16; // chf4He, chr4He は 0-15 の値を取る

    // "gomi.root" に出力するためのファイルを開く
    TFile *outputFile = new TFile("gomi.root", "RECREATE");

    // 二次元グラフ用の TGraph を作成
    TGraph *graph = new TGraph();
    int point_index = 0;

    // 微分断面積計算
    printf("chf\tchr\tYield\tSolid_Angle\tDiff_Cross_Section (mb/sr)\n");
    for (int chf = 0; chf < n_bins; ++chf) {
      for (int chr = 0; chr < n_bins; ++chr) {
      //        for (int chr = 0; chr < 1; ++chr) {
            // chf と chr に基づくヒストグラムを作成
            TH1D *h1 = new TH1D(TString::Format("h1_%d_%d", chf, chr), "Ex4He", 100, 6.5, 8.5);
            TString condition = TString::Format("chf4He==%d && chr4He==%d && Amax4He>-0.2 && Amax4He<0.15", chf, chr);
            tree->Draw("Ex4He>>" + TString(h1->GetName()), condition);

            // バックグラウンドの初期値をヒストグラムの両端から推定
            double bg_left = h1->GetBinContent(1);
            double bg_right = h1->GetBinContent(h1->GetNbinsX());
            double bg_slope = (bg_right - bg_left) / (8.5 - 6.5);
            double bg_intercept = bg_left - bg_slope * 6.5;

            // フィット関数を作成 (ガウシアン + 線形バックグラウンド)
            TF1 *fitFunc = new TF1("fitFunc", "[0]*exp(-0.5*((x-[1])/[2])**2) + [3]*x + [4]", 6.5, 8.5);
            fitFunc->SetParameters(200, 7.65, 0.05, bg_slope, bg_intercept); // 初期値の設定

            // フィット実行
            h1->Fit(fitFunc, "Q");

            // ガウシアンパラメータから収量を計算
            double amplitude = fitFunc->GetParameter(0);
            double sigma = fitFunc->GetParameter(2);
            double yield = amplitude * sqrt(2 * TMath::Pi()) * sigma / (2.0 / 100);

            // 立体角計算
            double area = get_area(chf, chr);
            double distance = get_distance(chf, chr);
            double solid_angle = area / (distance * distance);

            // 微分断面積の計算（修正、単位を mb に変換）
            double diff_cross_section = (solid_angle > 0) ? yield / (N_beam * n_effective * solid_angle) * 1e27 : 0.0;

            // front_theta を計算
            double front_theta = get_angle(chf,chr).first;

            // TGraph にデータポイントを追加
            if (diff_cross_section > 0 && diff_cross_section < 1e5) {
                graph->SetPoint(point_index, front_theta * 180.0 / TMath::Pi(), diff_cross_section); // θを度に変換
                ++point_index;
            }

            // 結果を表示
            printf("%d\t%d\t%.5f\t%.5e\t%.5f\n", chf, chr, yield, solid_angle, diff_cross_section);

            // フィッティング結果を "gomi.root" に保存
            h1->Write();

            // メモリ解放
            delete h1;
            delete fitFunc;
        }
    }

    // 二次元グラフを保存
    graph->Write("diff_cross_section_graph");

    // ファイルを閉じる
    outputFile->Close();

    // 二次元グラフを描画
    TCanvas *c = new TCanvas("c", "Differential Cross Section vs Front Theta", 800, 600);
    graph->SetTitle("Differential Cross Section vs Front Theta;Front Theta (deg);Differential Cross Section (mb/sr)");
    graph->SetMarkerStyle(20);
    graph->Draw("AP");
}
