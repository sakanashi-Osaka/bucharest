#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

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
    ifstream infile("~/exp/Bucharest2022/eventbuild_Dron/hit/macro/Dron/log/theta_444.txt");
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
