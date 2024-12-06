void plot_dist_444(){
  TString fileName = "gomi_with_errors444.root";
  TFile *file = TFile::Open(fileName);
  TString histName = "cross_section_cm_with_errors";
  TGraphErrors *g = (TGraphErrors*)file->Get(histName);


  TCanvas *c0 = new TCanvas("c0","c0",600,600);

  double xMin = 40.0, xMax = 75.0; // xの範囲
  double yMin = 0.7, yMax = 30.0; // yの範囲

  // 新しいTGraphを作成
  TGraphErrors *g0 = new TGraphErrors();
  int nPoints = g->GetN(); // 元のTGraphのデータ点数
  int newIndex = 0;
  
  for (int i = 0; i < nPoints; ++i) {
    double x, y, x_error, y_error;
    g->GetPoint(i, x, y);
    x_error = g->GetErrorX(i);
    y_error = g->GetErrorY(i);
    if (x >= xMin && x <= xMax && y >= yMin && y <= yMax) {
      g0->SetPoint(newIndex, x, y);
      g0->SetPointError(newIndex, x_error, y_error);
      newIndex++;
    }
  }
  
  g0->GetXaxis()->SetLimits(0,100);
  g0->GetYaxis()->SetLimits(0,100);
  g0->GetXaxis()->SetRangeUser(35,75);
  g0->GetYaxis()->SetRangeUser(0,8);
  g0->SetMarkerStyle(8);
  g0->SetMarkerSize(0.5);
  g0->Draw("AP");
}

