#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <TCanvas.h>
#include <TEventList.h>
#include <cstdint>
#include <vector>
#include <signal.h>
#include <map>
#include <algorithm>
#include <regex>
#include <TH1F.h>
#include <TH2.h>
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>

int time_cor(std::string fin);

std::vector< std::vector<double> > tmp_array(101);
std::vector< std::vector<double> > tmp_time(101);

int sorter_cor(int run){
    
  //  TFile *fin =new TFile(Form("../sorter/rootfile/run%d_-1_ssgant1.root",run));
  //    TFile *fin =new TFile(Form("../sorter_cor/rootfile/run%d_0_ssgant.root",run));
  TFile *fin =new TFile(Form("../sorter/rootfile/marged%d.root",run));
  //    TFile *fin =new TFile(Form("../data/marged%d.root",run));
  //  TFile *fin =new TFile(Form("./rootfile/run%d_-1.root",run));
  TTree *tree = (TTree*)fin->Get("tree");


  double array_array[20]={};
  double array_time[101][20]={};
  //  time_cor(Form("../tmp_time/log/cor/log%d.txt",run));
  time_cor(Form("../tmp_time/log/cor/cor2/log%d.txt",run)); //for run2236- rebuild
  for(int j=1; j<tmp_array.at(0).size(); j++){
    //    std::cout << tmp_array.at(0).at(j) << " ";                         
    array_array[j-1] = tmp_array.at(0).at(j);
  }
  
  for(int i=0; i<tmp_array.size(); i++){
    //    std::cout << (i+15)%16 << " ";                              
    for(int j=1; j<tmp_array.at(i).size(); j++){
      if(tmp_array.at(i).size()!=tmp_array.at(0).size()) continue;
      tmp_time.at(i).push_back(tmp_array.at(i).at(j)-tmp_array.at(0).at(j));

      bool flag = false;
      for(int k=-50; k<50; k++){
        if(abs(tmp_time.at(i).at(j-1)-2.147*k)<1){
	  //	  std::cout << k << " ";                                      
	  array_time[i][j-1] = k;
          flag = true;
        }
      }
      //      if(flag==false) std::cout << "error" << " ";                                      
    }
    //    std::cout << std::endl;
  }

  for(int j=0; j<20; j++){
    cout << array_array[j] << " ";
  }
  cout << endl;

  for(int i=0; i<100; i++){
    for(int j=0; j<20; j++){
      cout << array_time[i][j] << " ";
    }
    cout << endl;
  }

  
  int tmp_domain;
  int tmp_adc;
  double tmp_ts;
  //  double tmp_tdc;
  float tmp_energy;
  float tmp_amax;
  
  tree->SetBranchAddress("domain", &tmp_domain);
  tree->SetBranchAddress("ADC", &tmp_adc);
  tree->SetBranchAddress("FineTS", &tmp_ts);
  tree->SetBranchAddress("Energy", &tmp_energy);
  tree->SetBranchAddress("Amax", &tmp_amax);
  //  tree->SetBranchAddress("TDC", &tmp_tdc);
  
  
  //  TFile *fout =new TFile(Form("rootfile/run%d_0_ssgant.root",run),"recreate");
  TFile *fout =new TFile(Form("rootfile/run%d_-2.root",run),"recreate"); //after fixing ch48 build miss 
  TTree *tout = new TTree("tree","tree");  
  
  int domain;
  int ADC;
  double FineTS;
  //  double TDC;
  float Energy;
  float Amax;
  
  tout->Branch("domain",&domain,"domain/I");
  tout->Branch("FineTS",&FineTS,"FineTS/D");
  tout->Branch("ADC",&ADC,"ADC/I");
  tout->Branch("Energy",&Energy,"Energy/F");
  tout->Branch("Amax",&Amax,"Amax/F");
  //  tout->Branch("TDC",&TDC,"TDC/D");


  ULong64_t N=tree->GetEntries();
  cout << "Total entry: " << N << endl;
    for(ULong64_t evtn=0; evtn<N; evtn++){
  //  for(ULong64_t evtn=20e7; evtn<60e7; evtn++){

    int loop = (int)((evtn-evtn%(N/200))/(N/200));
      
    domain=-1;
    ADC=-1;
    FineTS=-1;
    Energy=-1;
    Amax=-1;
    //    TDC=-1;

    if(evtn%10000==0) cout << "\rAnalyzed entry:" << evtn << " loop:" << loop; std::cout << flush;

    tree->GetEntry(evtn);

    domain = tmp_domain;
    ADC = tmp_adc;
    FineTS = tmp_ts;
    //    TDC = tmp_tdc;

    
    if(run<2305){
      if(domain>62){
	for(int i=0; i<20; i++){
	  if(array_array[i]<0.1) continue;
	  double tmp_num = (double)array_time[domain-62][i]*pow(2,31)*1e3;
	  if(array_array[i]*1e12 + tmp_num-1e12 < tmp_ts){
	    FineTS = tmp_ts - tmp_num;
	  }
	}
      }
      if(domain==0 || domain==28 || domain==29 || domain==31 || domain==45 ||domain==46 ||  domain==61 || domain==63){
	for(int i=0; i<20; i++){
	  if(array_array[i]<0.1) continue;
	  double tmp_num = (double)array_time[1][i]*pow(2,31)*1e3;
	  if(array_array[i]*1e12 + tmp_num-1e12 < tmp_ts){
	    FineTS = tmp_ts - tmp_num;
	  }
	}
      }
    }

    if(domain==75){
      if(run==2238){
	if(FineTS > 2470*1e12 && FineTS<3205*1e12){
	  FineTS = FineTS - pow(2,31)*1e3 * 2;
	}
      }
    }
    if(domain==8 || domain==9){
      if(run==2244){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;      
      }
    }
    if(domain==149){
      if(run==2264){
	if(FineTS > 126*1e12 && FineTS<448*1e12){
	  FineTS = FineTS - pow(2,31)*1e3 * 2;
	}
      }
    }
    if(domain==2 || domain==3 || domain==6 || domain==7){
      if(run==2267){
	FineTS = FineTS - pow(2,31)*1e3 * 2;
      }
    }
    
    
    
    //time correction for run~2304
    if(domain==30 || domain==31){
      if(run==2300){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;      
      }
    }

    if(domain==2 || domain==3){
      if(run==2299){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;      
      }
    }

    if(domain==14 || domain==15){
      if(run==2294){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;      
      }
    }

    if(domain==0 || domain==4 || domain==5 || domain==8 || domain==9 || domain==10 || domain==11){
      if(run==2304){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;      
      }
    }

    if(domain==0 || domain==28 || domain==29 || domain==31 || domain==45 ||domain==46 ||  domain==61 || domain==63){      
      if(run==2287){
	if(tmp_ts > 266*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 684*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 2395*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 6;
      }
      if(run==2289){
	if(tmp_ts > 1892*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 12;
	if(tmp_ts > 2892*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 16;
      }
      if(run==2292){
	if(tmp_ts > 928*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 1192*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 8;
	if(tmp_ts > 1562*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 10;
	if(tmp_ts > 2560*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 12;
	if(tmp_ts > 2788*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
      }
      if(run==2293){
	if(tmp_ts > 472*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 2671*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 10;
      }
      if(run==2294){
	if(tmp_ts > 1473*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
      }
      if(run==2295){
	if(tmp_ts > 1021*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 2436*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
	if(tmp_ts > 2734*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 20;
	if(tmp_ts > 3119*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 24;
      } 
      if(run==2296){
	if(tmp_ts > 352*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 945*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 1379*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 6;
      }
     if(run==2297){
	if(tmp_ts > 584*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 1528*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 6;
      }
      if(run==2299){
	if(tmp_ts > 317*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 352*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 1105*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 16;
	if(tmp_ts > 2228*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 28;
	if(tmp_ts > 2294*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 32;
      }
      if(run==2300){
	if(tmp_ts > 1530*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
      }
      if(run==2301){
	if(tmp_ts > 1326*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 12;
	if(tmp_ts > 1662*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
      }
      if(run==2304){
	if(tmp_ts > 175*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 3092*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
      }
    }
    

    
    //time correction for run2305~
    if(domain==0 || domain==28 || domain==29 || domain==31 || domain==45 ||domain==46 ||  domain==61 || domain==63){
      
      if(run==2305){
	if(tmp_ts > 1903*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 2845*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
	if(tmp_ts > 2896*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 16;
      }
      if(run==2306){
	if(tmp_ts > 297*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 627*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
	if(tmp_ts > 687*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 16;
	if(tmp_ts > 1585*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 28;
	if(tmp_ts > 2612*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 32;
      }
      if(run==2307){
	if(tmp_ts > 3316*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 3321*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
      }
      if(run==2309){
	if(tmp_ts > 359*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 12;
      }
      if(run==2311){
	if(tmp_ts > 1977*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
      }
      if(run==2313){
	if(tmp_ts > 661*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 667*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 8;
	if(tmp_ts > 2248*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 10;
	if(tmp_ts > 4301*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 12;
      }
      if(run==2315){
	if(tmp_ts > 574*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 2732*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 6;
	if(tmp_ts > 3215*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 8;
	if(tmp_ts > 3444*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 20;
      }
      if(run==2316){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;
	if(tmp_ts > 2103*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 2109*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 8;
      }
      if(run==2317){
	if(tmp_ts > 118*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 12;
      }      
      if(run==2318){
	if(tmp_ts > 40*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 325*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 8;
	if(tmp_ts > 494*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 10;
	if(tmp_ts > 501*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 12;
	if(tmp_ts > 730*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 16;
	if(tmp_ts > 2343*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 18;
      }      
      if(run==2320){
	if(tmp_ts > 2143*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
      }      
      if(run==2321){
	if(tmp_ts > 637*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 12;
      }      
      if(run==2322){
	if(tmp_ts > 730*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 767*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
	if(tmp_ts > 2038*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 16;
      }      
      if(run==2323){
	if(tmp_ts > 264*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 12;
	if(tmp_ts > 475*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
	if(tmp_ts > 2311*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 16;
	if(tmp_ts > 3033*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 28;
      }      
      if(run==2324){
	if(tmp_ts > 1085*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
      }      
      if(run==2325){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;
      }      
      if(run==2327){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;
	if(tmp_ts > 703*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 0;
	if(tmp_ts > 1839*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 12;
	if(tmp_ts > 2239*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
      }      
      if(run==2328){
	if(tmp_ts > 495*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 1090*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 1192*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 8;
	if(tmp_ts > 1830*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 12;
	if(tmp_ts > 2354*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
	if(tmp_ts > 2558*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 18;
	if(tmp_ts > 2807*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 20;
      }      
      if(run==2329){
	if(tmp_ts > 1486*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
      }
      if(run==2330){
	if(tmp_ts > 554*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
      }
      if(run==2331){
	if(tmp_ts > 1616*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 1796*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 6;
	if(tmp_ts > 2035*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 18;
	if(tmp_ts > 2288*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 24;
	if(tmp_ts > 2312*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 36;
	if(tmp_ts > 3200*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 38;
	if(tmp_ts > 3380*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 40;
      }      
      if(run==2332){
	if(tmp_ts > 554*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 12;
      } 
      if(run==2333){
	if(tmp_ts > 1745*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 12;
	if(tmp_ts > 2320*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
      }
      if(run==2335){
	if(tmp_ts > 65*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 1383*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
      }
      if(run==2337){
	if(tmp_ts > 300*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 1660*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 2500*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 6;
	if(tmp_ts > 3194*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 8;
      }
      if(run==2338){
	if(tmp_ts > 298*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 1410*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 6;
	if(tmp_ts > 2308*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 10;
	if(tmp_ts > 3051*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
      }
      if(run==2339){
	if(tmp_ts > 213*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 273*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
      }
      if(run==2340){
	if(tmp_ts > 1920*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 1925*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 6;
	if(tmp_ts > 2380*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 20;
      }
      if(run==2341){
	if(tmp_ts > 469*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
      }
      if(run==2342){
	if(tmp_ts > 645*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 651*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 6;
	if(tmp_ts > 2175*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 10;
      }
      if(run==2343){
	if(tmp_ts > 3274*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
      }	
      if(run==2344){
	if(tmp_ts > 458*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 12;
	if(tmp_ts > 891*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
	if(tmp_ts > 1502*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 16;
	if(tmp_ts > 2012*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 18;
	if(tmp_ts > 2504*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 20;
	if(tmp_ts > 2510*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 24;
	if(tmp_ts > 3366*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 26;
      }
      if(run==2345){
	if(tmp_ts > 1005*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 2261*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 2536*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 6;
	if(tmp_ts > 3450*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 10;
      }
      if(run==2347){
	if(tmp_ts > 357*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 12;
	if(tmp_ts > 2652*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
      }
      if(run==2348){
	if(tmp_ts > 1001.5*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
      }
      if(run==2349){
	if(tmp_ts > 680*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 1804*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 78;
      }
      if(run==2350){
	if(tmp_ts > 253*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 854*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
	if(tmp_ts > 1482*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 16;
	if(tmp_ts > 1735*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 18;
	if(tmp_ts > 2026*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 20;
	if(tmp_ts > 2980*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 24;
      }
      if(run==2351){
	if(tmp_ts > 594*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 8;
	if(tmp_ts > 1994*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 12;
	if(tmp_ts > 2330*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 16;
	if(tmp_ts > 2577*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 20;
	if(tmp_ts > 3068*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 22;
      }
      if(run==2352){
	if(tmp_ts > 1155*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 1410*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
      }
      if(run==2353){
	if(tmp_ts > 84*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 364*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 1477*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 6;
	if(tmp_ts > 3143*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 10;
      }
      if(run==2354){
	if(tmp_ts > 992*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 6;
      }
      if(run==2355){
	if(tmp_ts > 2125*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
      }
      if(run==2356){
	if(tmp_ts > 806*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
      }
      if(run==2357){
	if(tmp_ts > 189*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 8;
	if(tmp_ts > 717*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 10;
	if(tmp_ts > 1387*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 12;
	if(tmp_ts > 2975*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
      }
      if(run==2358){
	if(tmp_ts > 390*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 1954*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 12;
	if(tmp_ts > 2271*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
	if(tmp_ts > 3013*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 16;
      }
      if(run==2359){
	if(tmp_ts > 292*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 2560*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 2874*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 8;
	if(tmp_ts > 3421*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 10;
      }
      if(run==2360){
	if(tmp_ts > 72*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 1104*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 2419*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 6;
	if(tmp_ts > 2794*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 8;
	if(tmp_ts > 2894*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 20;
	if(tmp_ts > 3275*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 22;
      }
      if(run==2361){
	if(tmp_ts > 1865*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 2347*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 2562*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 16;
	if(tmp_ts > 2711*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 18;
      }
      if(run==2362){
	if(tmp_ts > 150*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 898*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 12;
	if(tmp_ts > 2054*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
	if(tmp_ts > 2908*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 26;
	if(tmp_ts > 3319*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 28;
      }
      if(run==2363){
	if(tmp_ts > 865*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 1179*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
	if(tmp_ts > 2325*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 16;
	if(tmp_ts > 2688*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 20; 
	if(tmp_ts > 2773*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 22;
	if(tmp_ts > 3360*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 32;
      }
      if(run==2369){
	if(tmp_ts > 25*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
      }
      if(run==2371){
	if(tmp_ts > 340*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 1111*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 8;
	if(tmp_ts > 1510*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 10;
	if(tmp_ts > 2300*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 12; 
	if(tmp_ts > 2544*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
	if(tmp_ts > 2775*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 26;
	if(tmp_ts > 3233*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 28;
      }
      if(run==2372){
	if(tmp_ts > 42*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 495*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 6;
	if(tmp_ts > 3085*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 8;
	if(tmp_ts > 3210*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 10; 
	if(tmp_ts > 3510*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 20;
      }
      if(run==2373){
	if(tmp_ts > 235*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 8;
      }
      if(run==2374){
	if(tmp_ts > 3507*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 4106*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 10;
      }
      if(run==2375){
	if(tmp_ts > 947*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 2034*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 8;
	if(tmp_ts > 3121*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 10;
	if(tmp_ts > 3231*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 18;
      }
      if(run==2376){
	if(tmp_ts > 2182*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 3361*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
      }
      if(run==2377){
	if(tmp_ts > 1874*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 12;
	if(tmp_ts > 2798*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 16;
      }
      if(run==2378){
	if(tmp_ts > 384*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
      }
      if(run==2379){
	if(tmp_ts > 1427*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 1984*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 6;
      }
      if(run==2381){
	if(tmp_ts > 1158*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 1164*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 6;
      }
      if(run==2384){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;
	if(tmp_ts > 76*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 0;
	if(tmp_ts > 2368*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
      }
      if(run==2385){
	if(tmp_ts > 378*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 382*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 8;
	if(tmp_ts > 3010*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 10;
      }
      if(run==2387){
	if(tmp_ts > 3371*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 3376*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 10;
      }
      if(run==2389){
	if(tmp_ts > 2*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 2384*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
	if(tmp_ts > 2905*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 16;
      }
      if(run==2390){
	if(tmp_ts > 3066*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 3074*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 8;
      }
      if(run==2391){
	if(tmp_ts > 601*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 2387*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
      }
      if(run==2394){
	if(tmp_ts > 132*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 211*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 365*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 8;
	if(tmp_ts > 370*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 10;
	if(tmp_ts > 787*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
	if(tmp_ts > 792*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 20;
	if(tmp_ts > 1334*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 22;
	if(tmp_ts > 1684*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 24;
      }
      if(run==2398){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 10;
      }
      if(run==2403){
	if(tmp_ts > 240*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 249*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 6;
      }
      if(run==2404){
	if(tmp_ts > 1582*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 2707*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 8;
      }
      if(run==2405){
	if(tmp_ts > 435*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 2640*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 8;
	if(tmp_ts > 2648*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
      }
      if(run==2406){
	if(tmp_ts > 1127*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;
	if(tmp_ts > 1481*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
	if(tmp_ts > 2750*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 6;
	if(tmp_ts > 3106*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 8;
	if(tmp_ts > 3114*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 12;
	if(tmp_ts > 3436*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
	if(tmp_ts > 3444*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 18;
      }
      if(run==2407){
	if(tmp_ts > 843*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 12;
      }
      if(run==2410){
	if(tmp_ts > 614*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 4;
      }
      if(run==2411){
	if(tmp_ts > 170*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 6;
	if(tmp_ts > 871*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 14;
	if(tmp_ts > 1496*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 26;
      }
      if(run==2413){
	if(tmp_ts > 1105*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 6;
	if(tmp_ts > 1265*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 8;
	if(tmp_ts > 1287*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 16;
      }      
    }
    
    if( (domain>1 && domain<28) || (domain>31 && domain<44) || (domain>47 &&  domain<60) ){
      if(run==2316){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;
      }
      if(run==2325){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;
      }
      if(run==2327){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;
      }
      if(run==2384){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;
      }
      if(run==2392){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;
      }
      if(run==2411){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;
      }
    }


    if(domain==20||domain==34||domain==42||domain==52||domain==58){
      if(run==2411){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 0;  
      }
    }	
    
    if(domain==0){
      if(run==2392){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;      
	if(tmp_ts > 1923*1e12) FineTS = tmp_ts + pow(2,31)*1e3 * 2;      
      }
    }

    if(domain==2 || domain==3){
      if(run==2385){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;      
      }
    }
    
    if(domain==8 || domain==9){
      if(run==2313){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;      
      }
      if(run==2316){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 0;  
      }
      if(run==2325){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 0;  
      }
      if(run==2327){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 0;  
      }
      if(run==2335){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;  
      }
      if(run==2348){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;      
      }
      if(run==2353){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;      
      }
      if(run==2392){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 0;      
      }
    }
    
    if(domain==14 || domain==15){
      if(run==2321){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;      
      }
    }
    
    if(domain==30 || domain==31){
      if(run==2348){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;      
      }
      if(run==2360){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;      
      }
      if(run==2377){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;      
      }
      if(run==2392){
	if(tmp_ts > 0*1e12) FineTS = tmp_ts - pow(2,31)*1e3 * 2;      
      }
    }

    Energy = tmp_energy;
    Amax = tmp_amax;

    tout->Fill();
  }
  
  tout->AutoSave();
  fout->Close();
  
  return 0;
}

int time_cor(std::string fin){

  std::ifstream ifs(fin);

  std::string str = "";
  int i = 0;
  
  while(getline(ifs, str)){
    std::string tmp = "";
    std::istringstream stream(str);
    while(getline(stream, tmp, ' ')){
      tmp_array.at(i).push_back(std::stod(tmp));
    }
    i++;
  }

  for(int j=1; j<tmp_array.at(0).size(); j++){
    //    std::cout << tmp_array.at(0).at(j) << " ";
  }
  //  std::cout << std::endl;
  
  for(int i=0; i<tmp_array.size(); i++){
    //    std::cout << (i+15)%16 << " ";
    for(int j=1; j<tmp_array.at(i).size(); j++){
      if(tmp_array.at(i).size()!=tmp_array.at(0).size()) continue;
      tmp_time.at(i).push_back(tmp_array.at(i).at(j)-tmp_array.at(0).at(j));
         
      bool flag = false;
      for(int k=-50; k<50; k++){
        if(abs(tmp_time.at(i).at(j-1)-2.147*k)<1){
	  //          std::cout << k << " ";
          flag = true;
        }
      }
      //      if(flag==false) std::cout << "error" << " ";

    }
    //    std::cout << std::endl;
  }

  return 0;
}
