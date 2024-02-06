TGraph* vlocate_stop(){
  if(!gPad){
    return NULL;
  }
  else{
    TGraph *g = (TGraph*)gPad->GetListOfPrimitives()->Last();
    int n = g->GetN();
    double *x = g->GetX();
    double *y = g->GetY();
    for(int i = 0; i < n; ++i){
      std::cout << x[i] << "," << y[i] << std::endl;
    }
    return g;
  }
}
