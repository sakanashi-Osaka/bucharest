double tmp(double thetacm_i=40., double thetacm_f=75., int n_div=350, int i_div=0);

double thetacm(){
  double thetacm_i = 40.;
  double thetacm_f = 75.;
  const int n_div = 350;
  double theta_array[n_div];

  double total=0.;
  for(int i=0; i<n_div; i++){
    theta_array[i] = tmp(thetacm_i,thetacm_f,n_div,i);
    total += theta_array[i];
  }
  cout << total;
  
  return 0;
}


double tmp(double thetacm_i=40., double thetacm_f=75., int n_div=350, int i_div=0){
  double p0 = -241.512; //fitted using hist.C
  double p1 = 3.3118;
  double p2 = 112.729;
  double p3 = 54.2054;
  double p4 = 8.41897;
  double p5 = 130.893;
  double p6 = 36.5291;
  double p7 = 5.08719;
  
  double x = thetacm_i + (thetacm_f-thetacm_i)/n_div * (double)i_div;
  double num;
  num = p0 + p1*x + p2*exp(-pow(x-p3,2)/(2*pow(p4,2))) + p5*exp(-pow(x-p6,2)/(2*pow(p7,2)));
    
  return num;
}
