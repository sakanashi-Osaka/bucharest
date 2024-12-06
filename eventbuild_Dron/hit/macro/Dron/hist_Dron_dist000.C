#include <cmath>
#include <vector>


const double s1_r1 = 48.0 / 2.0;  // 内側の半径
const double s1_r2 = 96.0 / 2.0;  // 外側の半径
const double s1_dist = 41.0;      // 検出器までの距離
const double beam_x = 0.0, beam_y = 2.4;
const double PI = TMath::Pi();


// データ格納用の構造体
struct AngleData {
    double labAngle;
    double cmAngle;
};
// 一次補間関数
double Interpolate(double x1, double y1, double x2, double y2, double x) {
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
}

// 実験室系から重心系の角度を計算する関数
double Lab_to_CM(double labAngle) {
    // データ読み込み
    ifstream infile("log/theta_000.txt");
    if (!infile) {
        cerr << "Error: File could not be opened." << endl;
        return -1; // エラー値
    }

    vector<AngleData> data;
    double lab, cm;
    while (infile >> lab >> cm) {
        data.push_back({lab, cm});
    }
    infile.close();

    // データが降順の場合は逆順に並び替える
    if (!data.empty() && data.front().labAngle > data.back().labAngle) {
        reverse(data.begin(), data.end());
    }

    // 範囲外のチェック
    if (labAngle < data.front().labAngle || labAngle > data.back().labAngle) {
        cerr << "Error: Lab angle is out of data range." << endl;
        return -1; // エラー値
    }

    // データに直接一致する場合
    for (size_t i = 0; i < data.size(); ++i) {
        if (fabs(data[i].labAngle - labAngle) < 1e-6) {
            return data[i].cmAngle;
        }
    }
    // 線形補間
    for (size_t i = 0; i < data.size() - 1; ++i) {
        if (labAngle > data[i].labAngle && labAngle < data[i + 1].labAngle) {
            return Interpolate(data[i].labAngle, data[i].cmAngle,
                               data[i + 1].labAngle, data[i + 1].cmAngle,
                               labAngle);
        }
    }
    cerr << "Error: Unexpected behavior in LabToCM function." << endl;
    return -1; // エラー値
}


double gfactor(double theta_lab_deg) {
  //  double gfac = theta_lab_deg*0.0049 + 0.3828; // 7.65
  //  double gfac = theta_lab_deg*0.0046 + 0.4276; // 4.44
  double gfac = theta_lab_deg*0.0043 + 0.4728; // g.s.
  return gfac;
}

double rho_target = 2.267;  // 炭素の密度 (g/cm^3)
double M_target = 12.01;  // 炭素のモル質量 (g/mol)
double thickness_target = 50e-6;  // 標的の厚さ (50 μg/cm^2 -> g/cm^2)
double N_A = 6.022e23;  // アボガドロ定数
double n_target = (rho_target / M_target) * N_A;  // 原子密度 (atoms/cm^3)
double n_effective = n_target * thickness_target;  // 有効原子密度 (atoms/cm^2)

  //double beam_current_nAs = 663022;  // ビーム量 (nA*s) //12C_all
//double beam_current_nAs = 8800;  // ビーム量 (nA*s) //2305
double beam_current_nAs = 59643;  // ビーム量 (nA*s) //2297-2306
double beam_charge = 2;  //4He2+ は2価
double beam_current_C = beam_current_nAs*1e-9;  // nA*s を C に変換
double e_charge = 1.60217e-19; //電子1個あたりの電荷
double N_beam = beam_current_C / (beam_charge * e_charge); //ビーム粒子数 



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


#include <TGraphErrors.h>

void hist_Dron_dist000() {
    // ROOTファイルを開く
  //    TFile *file = new TFile("../rootfile/gamma2305.root");
  TFile *file = new TFile("../rootfile/gamma_2297_2306.root");
    TTree *tree = (TTree*)file->Get("tree");

    const int n_bins = 16; // chf4He, chr4He は 0-15 の値を取る

    // "gomi.root" に出力するためのファイルを開く
    TFile *outputFile = new TFile("gomi_with_errors000.root", "RECREATE");

    // 二次元グラフ用の TGraphErrors を作成
    TGraphErrors *graph = new TGraphErrors();
    TGraphErrors *graph2 = new TGraphErrors();
    int point_index = 0;

    // 微分断面積計算
    printf("chf\tchr\tYield\tSolid_Angle\tDiff_Cross_Section (mb/sr)\n");
    for (int chf = 0; chf < 16; ++chf) {
        for (int chr = 0; chr < 16; ++chr) {
            // chf と chr に基づくヒストグラムを作成
            TH1D *h1 = new TH1D(TString::Format("h1_%d_%d", chf, chr), "Ex4He", 100, -1, 1);
            TString condition = TString::Format("chf4He==%d && chr4He==%d && Amax4He>-0.2 && Amax4He<0.2", chf, chr);
            tree->Draw("Ex4He>>" + TString(h1->GetName()), condition);

            // バックグラウンドの初期値をヒストグラムの両端から推定
            double bg_left = h1->GetBinContent(1);
            double bg_right = h1->GetBinContent(h1->GetNbinsX());
            double bg_slope = (bg_right - bg_left) / (1. + 1.);
            double bg_intercept = bg_left - bg_slope * (-1);

            int maxBin = h1->GetMaximumBin();
            double maxBinHeight = h1->GetBinContent(maxBin);
            double val_maxBin = -1. + 0.02 * (double)maxBin;

            // フィット関数を作成 (ガウシアン + 線形バックグラウンド)
            TF1 *fitFunc = new TF1("fitFunc", "[0]*exp(-0.5*((x-[1])/[2])**2) + [3]*x + [4]", -1, 1);
            fitFunc->SetParameters(maxBinHeight, val_maxBin, 0.1, bg_slope, bg_intercept); // 初期値の設定

            // フィット実行
            h1->Fit(fitFunc, "Q");

            // ガウシアンパラメータから収量を計算
            double amplitude = fitFunc->GetParameter(0);
            double sigma = fitFunc->GetParameter(2);
            double yield = amplitude * sqrt(2 * TMath::Pi()) * sigma / (2.0 / 100);

            // 収量の統計誤差を計算
            double yield_error = amplitude > 0 ? sqrt(amplitude) * sqrt(2 * TMath::Pi()) * sigma / (2.0 / 100) /2: 0;
	    double theta_error = (get_angle(chf,chr).first-get_angle(chf+1,chr).first)/4 * 180/PI;
	    
            // 立体角計算
            double area = get_area(chf, chr);
            double distance = get_distance(chf, chr);
            double solid_angle = area / (distance * distance);

            // 微分断面積とその誤差を計算（単位を mb に変換）
            double diff_cross_section = (solid_angle > 0) ? yield / (N_beam * n_effective * solid_angle) * 1e27 : 0.0;
            double diff_cross_section_error = (solid_angle > 0 && yield > 0) ? yield_error / (N_beam * n_effective * solid_angle) * 1e27 : 0.0;

            // front_theta を計算
            double front_theta = get_angle(chf, chr).first;

            // TGraphErrors にデータポイントを追加
            double tmp_theta = front_theta * 180.0 / TMath::Pi();
            if (diff_cross_section > 0.05 && diff_cross_section < 100) {
                if (chr == 3 || chr == 4 || chr == 11 || chr == 12) {
                    graph->SetPoint(point_index, Lab_to_CM(tmp_theta), diff_cross_section * gfactor(tmp_theta)); // θを度に変換
                    graph->SetPointError(point_index, theta_error, diff_cross_section_error * gfactor(tmp_theta)); // 誤差設定
                } else {
                    graph2->SetPoint(point_index, Lab_to_CM(tmp_theta), diff_cross_section * gfactor(tmp_theta)); // θを度に変換
                    graph2->SetPointError(point_index, theta_error, diff_cross_section_error * gfactor(tmp_theta)); // 誤差設定
                }
                ++point_index;
            }

            // 結果を表示
            printf("%d\t%d\t%.5f\t%.5e\t%.5f\t%.5e\t%.5f\n", chf, chr, yield, front_theta * 180.0 / TMath::Pi(), diff_cross_section, Lab_to_CM(front_theta * 180.0 / TMath::Pi()), diff_cross_section * gfactor(front_theta * 180.0 / TMath::Pi()));

            // フィッティング結果を "gomi_with_errors.root" に保存
            h1->Write();

            // メモリ解放
            delete h1;
            delete fitFunc;
        }
    }

    // 二次元グラフを保存
    graph->Write("tmp_with_errors");
    graph2->Write("cross_section_cm_with_errors");

    // ファイルを閉じる
    outputFile->Close();

    // 二次元グラフを描画
    TCanvas *c = new TCanvas("c", "Differential Cross Section with Errors", 800, 600);
    graph2->SetTitle("Differential Cross Section with Errors;Front Theta (deg);Differential Cross Section (mb/sr)");
    graph2->SetMarkerStyle(4);
    graph2->SetMarkerSize(0.5);
    graph2->SetMarkerColor(1);
    graph2->Draw("AP");
    graph->SetMarkerStyle(4);
    graph->SetMarkerSize(0.5);
    graph->SetMarkerColor(2);
    graph->Draw("P SAME");
}





