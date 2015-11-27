#include <vector>
#include <string>
#include <iostream>
#include <TMath.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "untuplizer.h"

bool isPassZmumu(TreeReader &data, vector<Int_t>& goodMuID){
  
  goodMuID.clear();

  Int_t    nMu         = data.GetInt("nMu");
  Int_t*   muCharge    = data.GetPtrInt("muCharge");
  Float_t* muMiniIsoEA = data.GetPtrFloat("muMiniIsoEA");
  TClonesArray* muP4   = (TClonesArray*) data.GetPtrTObject("muP4");
  vector<bool>& isHighPtMuon        = *((vector<bool>*) data.GetPtr("isHighPtMuon"));
  vector<bool>& isCustomTrackerMuon = *((vector<bool>*) data.GetPtr("isCustomTrackerMuon"));

  // select good muons
      
  std::vector<Int_t> goodMuons;
  bool hasTrigMuon = false;  

  for(Int_t im = 0; im < nMu; im++){
      
    TLorentzVector* myMu = (TLorentzVector*)muP4->At(im);

    if( muMiniIsoEA[im] > 0.2 ) continue;
    if( !isHighPtMuon[im] && !isCustomTrackerMuon[im] ) continue;
    if( fabs(myMu->Eta()) > 2.4 ) continue;
    if( myMu->Pt() < 20 ) continue;
    if( myMu->Pt() > 50 && fabs(myMu->Eta()) < 2.1 ) hasTrigMuon = true;

    goodMuons.push_back(im);

  }	

  if( !hasTrigMuon ) goodMuons.clear();

  // select good Z boson

  bool findMPair = false;
  TLorentzVector* thisMu = NULL;
  TLorentzVector* thatMu = NULL;

  for(unsigned int i = 0; i < goodMuons.size(); i++){

    Int_t im = goodMuons[i];
    thisMu = (TLorentzVector*)muP4->At(im);

    for(unsigned int j = 0; j < i; j++){

      Int_t jm = goodMuons[j];
      thatMu = (TLorentzVector*)muP4->At(jm);

      Float_t pt1  = thisMu->Pt();
      Float_t pt2  = thatMu->Pt();
      Float_t mll  = (*thisMu+*thatMu).M();
      Float_t ptll = (*thisMu+*thatMu).Pt();

      if( muCharge[im]*muCharge[jm] > 0 ) continue;
      if( mll < 70 || mll > 110 ) continue;
      if( ptll < 200 ) continue;
      if( !( (isHighPtMuon[im] && isCustomTrackerMuon[jm]) || (isHighPtMuon[jm] && isCustomTrackerMuon[im]) ) ) continue;

      if( !findMPair ){

	goodMuID.push_back( (pt1 > pt2) ? im : jm );
	goodMuID.push_back( (pt1 > pt2) ? jm : im );

      }

      findMPair = true;
      break;

    }
  }

  if( !findMPair ) return false;
  return true;

}
