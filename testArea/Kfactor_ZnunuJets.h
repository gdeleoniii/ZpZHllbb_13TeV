#include <iostream>
#include <string>
#include <TH1.h>
#include <TFile.h>

Double_t kfactorWeight(Float_t HT){

  const Double_t varBins[] = {100,200,400,600,10000000};

  TH1D* h = new TH1D("h","", 4, varBins);

  // for ZJetsToNuNu

  Double_t kfactor[4] = {1.626,
			 1.617,
			 1.459,
			 1.391};

  // HT: The scalar sum pt of the outgoing parton ( product of hard collisions, not including those from pileups)

  Double_t k1 = kfactor[h->FindBin(HT)-1];

  /*
    TFile* inf = new TFile("scalefactors_v4.root");
    TF1 fewk_z = (TF1*)inf->Get("z_ewkcorr/z_ewkcorr_func");
    TF1 fewk_w = (TF1*)inf->Get("w_ewkcorr/w_ewkcorr_func");
  */



  return k1*k2;

}
