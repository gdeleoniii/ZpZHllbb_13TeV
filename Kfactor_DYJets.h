#include <iostream>
#include <string>
#include <vector>
#include <TH1.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "untuplizer.h"

Double_t kfactorWeight(TreeReader &data, TF1* fewk_z){

  Int_t         nGenPar     = data.GetInt("nGenPar"); 
  Int_t*        genParId    = data.GetPtrInt("genParId");
  Int_t*        genMomParId = data.GetPtrInt("genMomParId");
  Int_t*        genParSt    = data.GetPtrInt("genParSt");
  Float_t       HT          = data.GetFloat("HT");
  TClonesArray* genParP4    = (TClonesArray*)data.GetPtrTObject("genParP4");

  // LO->NLO correction

  const Double_t varBins[] = {0,200,400,600,1e10};

  TH1D* h = new TH1D("h","", 4, varBins);

  // for DYJetsToLL

  Double_t kfactor[4] = {1.588,
			 1.438,
			 1.494,
			 1.139};

  // HT: The scalar sum pt of the outgoing parton ( product of hard collisions, not including those from pileups)

  Double_t k1 = kfactor[h->FindBin(HT)-1];

  h->Clear();  

  // NLO->NLO+EW correction

  /* 

     Remember put this two lines in your main function:

     TFile* f = TFile::Open("scalefactors_v4.root");
     TF1* fewk_z = (TF1*)(f->Get("z_ewkcorr/z_ewkcorr_func"));

  */

  vector<Int_t> goodLepID;

  for(Int_t ig = 0; ig < nGenPar; ig++){

    Int_t PID    = genParId[ig];
    Int_t momPID = genMomParId[ig];
    Int_t status = genParSt[ig];

    if( abs(PID) != 11 && 
	abs(PID) != 13 && 
	abs(PID) != 15 ) continue;

    if( (status != 1 && abs(PID) != 15) || 
	(status != 2 && abs(PID) == 15) ) continue;

    if( momPID != 23 && momPID != PID ) continue;

    goodLepID.push_back(ig);

  }

  TLorentzVector* l4_thisLep = (TLorentzVector*)genParP4->At(goodLepID[0]);
  TLorentzVector* l4_thatLep = (TLorentzVector*)genParP4->At(goodLepID[1]);

  TLorentzVector l4_z = (*l4_thisLep+*l4_thatLep);

  Double_t k2 = fewk_z->Eval(l4_z.Pt());

  return k1*k2;

}
