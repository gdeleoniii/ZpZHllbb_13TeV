#include <vector>
#include <iostream>
#include <TMath.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "untuplizer.h"

bool isPassZee(TreeReader &data, vector<Int_t>& goodEleID){

  goodEleID.clear();

  const Int_t    nEle       = data.GetInt("nEle");
  const Int_t*   eleCharge  = data.GetPtrInt("eleCharge");
  const Float_t* eleScEta   = data.GetPtrFloat("eleScEta");
  const TClonesArray* eleP4 = (TClonesArray*) data.GetPtrTObject("eleP4");
  const vector<bool>& eleIsPassLoose = *((vector<bool>*) data.GetPtr("eleIsPassLoose"));

  // select good electrons
        
  std::vector<Int_t> goodElectrons;

  for( int ie = nEle-1; ie >= 0; --ie ){

    TLorentzVector* myEle = (TLorentzVector*)eleP4->At(ie);

    if( fabs(eleScEta[ie]) > 1.4442 && fabs(eleScEta[ie]) < 1.566 ) continue;
    if( fabs(eleScEta[ie]) > 2.5 ) continue;
    if( fabs(myEle->Eta()) > 2.5 ) continue;
    if( myEle->Pt() < 20 ) continue;
    if( !eleIsPassLoose[ie] ) continue;
   
    goodElectrons.push_back(ie);

  } // end of ele loop

  // select good Z boson

  bool findEPair = false;
  TLorentzVector *thisEle = NULL, *thatEle = NULL;

  for( unsigned int i = 0; i < goodElectrons.size(); ++i ){

    int ie = goodElectrons[i];
    thisEle = (TLorentzVector*)eleP4->At(ie);

    for( unsigned int j = 0; j < i; ++j ){

      int je = goodElectrons[j];
      thatEle = (TLorentzVector*)eleP4->At(je);

      if( eleCharge[ie]*eleCharge[je] > 0 ) continue;   
      if( TMath::Max(thisEle->Pt(),thatEle->Pt()) < 135 ) continue;
      if( (*thisEle+*thatEle).M() < 70 || (*thisEle+*thatEle).M() > 110 ) continue;
      if( (*thisEle+*thatEle).Pt() < 200 ) continue;

      if( !findEPair ){

	goodEleID.push_back( (thisEle->Pt() > thatEle->Pt()) ? ie : je );
	goodEleID.push_back( (thisEle->Pt() > thatEle->Pt()) ? je : ie );

      }

      findEPair = true;
      break;

    }
  }

  return findEPair ? true : false;

}
