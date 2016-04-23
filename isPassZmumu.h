#include <vector>
#include <iostream>
#include <TMath.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "untuplizer.h"

bool isPassZmumu(TreeReader &data, vector<Int_t>& goodMuID){
  
  goodMuID.clear();

  const Int_t    nMu          = data.GetInt("nMu");
  const Int_t*   muCharge     = data.GetPtrInt("muCharge");
  const Float_t* muTrkIso     = data.GetPtrFloat("muTrkIso");
  const Float_t* muInnerTrkPt = data.GetPtrFloat("muInnerTrkPt");
  const TClonesArray* muP4    = (TClonesArray*) data.GetPtrTObject("muP4");
  const vector<bool>& isHighPtMuon        = *((vector<bool>*) data.GetPtr("isHighPtMuon"));
  const vector<bool>& isCustomTrackerMuon = *((vector<bool>*) data.GetPtr("isCustomTrackerMuon"));

  // select good muons
      
  std::vector<Int_t> goodMuons;
  bool hasTrigMuon = false;  

  for( int im = nMu-1; im >= 0; --im ){
      
    TLorentzVector* myMu = (TLorentzVector*)muP4->At(im);

    if( !isHighPtMuon[im] && !isCustomTrackerMuon[im] ) continue;
    if( fabs(myMu->Eta()) > 2.4 ) continue;
    if( myMu->Pt() < 20 ) continue;
    if( myMu->Pt() > 50 && fabs(myMu->Eta()) < 2.1 ) hasTrigMuon = true;

    goodMuons.push_back(im);

  }	

  if( !hasTrigMuon ) goodMuons.clear();

  // select good Z boson

  bool findMPair = false;
  TLorentzVector *thisMu = NULL, *thatMu = NULL;

  for( unsigned int i = 0; i < goodMuons.size(); ++i ){

    int im = goodMuons[i];
    thisMu = (TLorentzVector*)muP4->At(im);

    for( unsigned int j = 0; j < i; ++j ){

      int jm = goodMuons[j];
      thatMu = (TLorentzVector*)muP4->At(jm);

      // if the two muons are far away, use regular isolation

      if( thisMu->DeltaR(*thatMu) > 0.3 && muTrkIso[im]/thisMu->Pt() > 0.1 ) continue;
      if( thisMu->DeltaR(*thatMu) > 0.3 && muTrkIso[jm]/thatMu->Pt() > 0.1 ) continue;

      // if the two muons are close, use corrected isolation

      if( thisMu->DeltaR(*thatMu) < 0.3 && (muTrkIso[im]-muInnerTrkPt[jm])/thisMu->Pt() > 0.1 ) continue;
      if( thisMu->DeltaR(*thatMu) < 0.3 && (muTrkIso[jm]-muInnerTrkPt[im])/thatMu->Pt() > 0.1 ) continue;

      if( muCharge[im]*muCharge[jm] > 0 ) continue;
      if( (*thisMu+*thatMu).M() < 70 || (*thisMu+*thatMu).M() > 110 ) continue;
      if( (*thisMu+*thatMu).Pt() < 200 ) continue;
      if( !( (isHighPtMuon[im] && isCustomTrackerMuon[jm]) || (isHighPtMuon[jm] && isCustomTrackerMuon[im]) ) ) continue;

      if( !findMPair ){

	goodMuID.push_back( (thisMu->Pt() > thatMu->Pt()) ? im : jm );
	goodMuID.push_back( (thisMu->Pt() > thatMu->Pt()) ? jm : im );

      }

      findMPair = true;
      break;

    }
  }

  return findMPair ? true : false;

}
