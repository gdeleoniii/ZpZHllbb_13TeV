#include <vector>
#include <string>
#include <iostream>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/untuplizer.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZmumu.h"
#include "BTagCalibrationStandalone.h"

float muBtagUnc(string inputFile, int cat, int jSF){

  // setup calibration and reader

  string region = (jSF==0) ? "central" : ((jSF==1) ? "up" : "down");

  BTagCalibration       calib ("csvv1", "CSVV1.csv");
  BTagCalibrationReader reader(BTagEntry::OP_LOOSE, region.data());
  
  reader.load(calib, BTagEntry::FLAV_B, "mujets");
  
  // read the ntuples (in pcncu)
  
  TreeReader data(inputFile.data());
  
  float btagScale = 1.;
  float passEvent = 0.;

  // begin of event loop

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    data.GetEntry(ev);

    TClonesArray*  muP4              = (TClonesArray*) data.GetPtrTObject("muP4");
    Int_t          FATnJet           = data.GetInt("FATnJet");    
    Int_t*         FATnSubSDJet      = data.GetPtrInt("FATnSubSDJet");
    Float_t*       FATjetPRmassCorr  = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    TClonesArray*  FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    vector<bool>&  FATjetPassIDLoose = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    vector<float>* FATsubjetSDCSV    = data.GetPtrVectorFloat("FATsubjetSDCSV", FATnJet);
    vector<int>*   FATsubjetSDHadronFlavor = data.GetPtrVectorInt("FATsubjetSDHadronFlavor", FATnJet);

    // select good reco level events     
    // select good leptons
      
    vector<int> goodLepID;

    if( !isPassZmumu(data,goodLepID) ) continue;

    TLorentzVector* thisLep = (TLorentzVector*)muP4->At(goodLepID[0]);
    TLorentzVector* thatLep = (TLorentzVector*)muP4->At(goodLepID[1]);

    // select good FATjet

    int goodFATJetID = -1;
    TLorentzVector thisJet(0,0,0,0);

    for( int ij = 0; ij < FATnJet; ++ij ){

      TLorentzVector* myJet = (TLorentzVector*)FATjetP4->At(ij);

      if( myJet->Pt() < 200 ) continue;
      if( fabs(myJet->Eta()) > 2.4 ) continue;
      if( !FATjetPassIDLoose[ij] ) continue;
      if( myJet->DeltaR(*thisLep) < 0.8 || myJet->DeltaR(*thatLep) < 0.8 ) continue;
      if( FATjetPRmassCorr[ij] < 105 || FATjetPRmassCorr[ij] > 135 ) continue;

      int nsubBjet = 0;

      for( int is = 0; is < FATnSubSDJet[ij]; ++is ){

	if( FATsubjetSDCSV[ij][is] > 0.605 ){

	  btagScale *= reader.eval_auto_bounds(region.data(), BTagEntry::FLAV_B, myJet->Eta(), myJet->Pt());
	  ++nsubBjet;

	}

	else{

	  if( FATsubjetSDHadronFlavor[ij][is] == 4 ) 

	    btagScale *= reader.eval_auto_bounds(region.data(), BTagEntry::FLAV_C, myJet->Eta(), myJet->Pt());

	  else if( FATsubjetSDHadronFlavor[ij][is] == 1 ||
		   FATsubjetSDHadronFlavor[ij][is] == 2 || 
		   FATsubjetSDHadronFlavor[ij][is] == 3 ||
		   FATsubjetSDHadronFlavor[ij][is] == 21 )	

	    btagScale *= reader.eval_auto_bounds(region.data(), BTagEntry::FLAV_UDSG, myJet->Eta(), myJet->Pt());

	}

      } // end of subjet for loop
 
      // b-tag cut
 
      if( cat == 1 && nsubBjet != 1 ) continue;
      if( cat == 2 && nsubBjet != 2 ) continue;
       
      goodFATJetID = ij;
      thisJet = *myJet;

      break;
 
    } // end of FatnJet loop
 
    if( goodFATJetID < 0 ) continue;

    if( (*thisLep+*thatLep+thisJet).M() < 750 ) continue;

    ++passEvent;

  } // end of event loop
  
  return passEvent*btagScale/(float)data.GetEntriesFast();

}
