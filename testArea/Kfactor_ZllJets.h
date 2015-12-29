#include <iostream>
#include <string>
#include <TH1.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "../untuplizer.h"

Double_t kfactorWeight(TreeReader &data){

  Int_t    nGenPar       = data.GetInt("nGenPar"); 
  Int_t*   genParId      = data.GetPtrInt("genParId");
  Int_t*   genMomParId   = data.GetPtrInt("genMomParId");
  Int_t*   genParSt      = data.GetPtrInt("genParSt");
  Float_t  HT            = data.GetFloat("HT");
  TClonesArray* genParP4 = (TClonesArray*)data.GetPtrTObject("genParP4");

  const Double_t varBins[] = {100,200,400,600,10000000};

  TH1D* h = new TH1D("h","", 4, varBins);

  // for DYJetsToLL

  Double_t kfactor[4] = {1.588,
			 1.438,
			 1.494,
			 1.139};

  // HT: The scalar sum pt of the outgoing parton ( product of hard collisions, not including those from pileups)

  Double_t k1 = kfactor[h->FindBin(HT)-1];
  
  TFile* inf = new TFile("scalefactors_v4.root");
  TF1* fewk_z = (TF1*)inf->Get("z_ewkcorr/z_ewkcorr_func");


  for(Int_t ig = 0; ig < nGenPar; ig++){

    bool isEle    = ( abs(genParId[ig] == 11) && genParSt[ig] == 1 );
    bool isMu     = ( abs(genParId[ig] == 13) && genParSt[ig] == 1 );
    bool isMomEle = ( abs(genMomParId[ig] == 11) );
    bool isMomMu  = ( abs(genMomParId[ig] == 13) );
    bool isMomZ   = ( genMomParId[ig] == 23 );

    TLorentzVector* thisGen = (TLorentzVector*)genParP4->At(ig);

    Double_t momElept = -1;
    Double_t momMupt = -1;
    Double_t momZpt = -1;

    if( isMomEle ) momElept = thisGen->Pt();
    if( isMomMu )  momMupt = thisGen->Pt();
    if( isMomZ )   momZpt = thisGen->Pt();

    //if( ( isEle || isMu ) && ( isMomEle || isMomMu || isMomZ ) )


  }
  

  return k1;

}
