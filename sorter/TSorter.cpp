#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <climits>

#include <TFile.h>
#include <TTree.h>

#include "TSorter.hpp"

using namespace std;

TSorter::TSorter() {}

TSorter::~TSorter()
{
  for(auto &&pointer: fDataVec) delete pointer;
}

void TSorter::LoadFiles(Int_t runNumber, Int_t start, Int_t stop)
{
  fRunNumber = runNumber;
  
  //Check file exist or not, then load
  if(stop < 0) stop = INT_MAX; 
  //  for(auto version = start; version <= stop; version++) {
  for(auto version = 0; version < 1; version++) {
    //    std::string fileName = fInputDir.Data() + std:: string("/run") + std::to_string(runNumber) + std::string("_") + std::to_string(version) + std::string("_ssgant.root");
    std::string fileName = fInputDir.Data() + std:: string("/marged") + std::to_string(runNumber) + std::string(".root");
    
    std::ifstream fin(fileName);
    if(!fin) {
      cout << "filename was wrong!! " << fileName << endl;
      fin.close();
      break;
    } else {
      std::cout << "Loading " << fileName << std::endl;
      fin.close();
      auto input = new TFile(fileName.c_str(), "READ");
      auto tree = (TTree *)input->Get("tree");
      if(tree) {
	const auto nEvents = tree->GetEntries();
	std::cout << nEvents << " events" << std::endl;

	TreeData treeData;
	tree->SetBranchAddress("domain", &(treeData.domain));
	tree->SetBranchAddress("FineTS", &(treeData.FineTS));
	tree->SetBranchAddress("ADC", &(treeData.ADC));
	tree->SetBranchAddress("Energy", &(treeData.Energy));
	tree->SetBranchAddress("Amax", &(treeData.Amax));
	//	tree->SetBranchAddress("TDC", &(treeData.TDC));
	
	for(auto iEve = 0; iEve < nEvents; iEve++){
	  tree->GetEntry(iEve);

	  auto oneHit = new TreeData;
	  oneHit->domain = treeData.domain;
	  oneHit->FineTS = treeData.FineTS;
	  oneHit->ADC = treeData.ADC;
	  oneHit->Energy = treeData.Energy;
	  oneHit->Amax = treeData.Amax;
	  //	  oneHit->TDC = treeData.TDC;

	  fDataVec.push_back(oneHit);
	}
	
      } else {
	std::cout << "Skip this file" << std::endl;
      }

      input->Close();
    }
  }
}

void TSorter::SortData()
{
  std::cout << "Sorting....." << std::endl;
  std::sort(fDataVec.begin(), fDataVec.end(),                                                            
	    [](const TreeData *a, const TreeData *b) {                                                   
	      return a->FineTS < b->FineTS;                                                                
	    });
}

void TSorter::WriteData(TString fileName)
{
  std::cout << "Writing data" << std::endl;

  if(fileName == "") {
    //    fileName = std::string("rootfile/run") + std::to_string(fRunNumber) + std::string("_") + std::to_string(-1) + std::string(".root");
    fileName = std::string("rootfile/marged") + std::to_string(fRunNumber) + std::string("_") + std::to_string(-1) + std::string(".root");
  }
  
    auto output = new TFile(fileName, "RECREATE");

    auto tree = new TTree("tree", "tree");
    int domain;
    double FineTS;
    //    double TDC;
    int ADC;
    float Energy;
    float Amax;

    tree->Branch("domain", &domain, "domain/I");
    tree->Branch("ADC", &ADC, "ADC/I");
    tree->Branch("FineTS", &FineTS, "FineTS/D");
    tree->Branch("Energy", &Energy, "Energy/F");
    tree->Branch("Amax", &Amax, "Amax/F");
    //    tree->Branch("TDC", &TDC, "TDC/D");

            
    for(auto iEve = 0; iEve < fDataVec.size(); iEve++) {
      if(fDataVec.at(iEve)->FineTS < 100) continue;
      domain = fDataVec.at(iEve)->domain;
      FineTS = fDataVec.at(iEve)->FineTS;
      ADC = fDataVec.at(iEve)->ADC;
      Energy = fDataVec.at(iEve)->Energy;
      Amax = fDataVec.at(iEve)->Amax;
      //      TDC = fDataVec.at(iEve)->TDC;
	
      tree->Fill();
    }

    tree->Write();
    output->Close();
  
}
