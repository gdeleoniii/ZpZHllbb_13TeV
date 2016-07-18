#include <vector>
#include <string>
#include <iostream>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/untuplizer.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZee.h"
#include "BTagCalibrationStandalone.h"

float eleBtag(string inputFile, int cat){

  // setup calibration and reader

  BTagCalibration       calib   ("csvv1", "CSVV1.csv");
  BTagCalibrationReader reader  (BTagEntry::OP_LOOSE, "central");
  BTagCalibrationReader readerUp(BTagEntry::OP_LOOSE, "up");
  BTagCalibrationReader readerDw(BTagEntry::OP_LOOSE, "down");
  
  reader  .load(calib, BTagEntry::FLAV_B, "comb");
  readerUp.load(calib, BTagEntry::FLAV_B, "comb");
  readerDw.load(calib, BTagEntry::FLAV_B, "comb");
  
  // read the ntuples (in pcncu)
  
  TreeReader data(inputFile.data());
  
  int passEvent = 0;

  // begin of event loop

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    data.GetEntry(ev);

    TClonesArray*  eleP4             = (TClonesArray*) data.GetPtrTObject("eleP4");
    Int_t          FATnJet           = data.GetInt("FATnJet");    
    Int_t*         FATnSubSDJet      = data.GetPtrInt("FATnSubSDJet");
    Float_t*       FATjetPRmassCorr  = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    TClonesArray*  FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    vector<bool>&  FATjetPassIDLoose = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    vector<float>* FATsubjetSDCSV    = data.GetPtrVectorFloat("FATsubjetSDCSV", FATnJet);

    // select good reco level events     
    // select good leptons
      
    vector<int> goodLepID;

    if( !isPassZee(data,goodLepID) ) continue;

    TLorentzVector* thisLep = (TLorentzVector*)eleP4->At(goodLepID[0]);
    TLorentzVector* thatLep = (TLorentzVector*)eleP4->At(goodLepID[1]);

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

      float jetSF   = reader  .eval_auto_bounds("central", BTagEntry::FLAV_B, myJet->Eta(), myJet->Pt());
      float jetSFup = readerUp.eval_auto_bounds("up",      BTagEntry::FLAV_B, myJet->Eta(), myJet->Pt());
      float jetSFdw = readerDw.eval_auto_bounds("down",    BTagEntry::FLAV_B, myJet->Eta(), myJet->Pt());

      fprintf(stdout, "%g ** %g ** %g\n", jetSF, jetSFup, jetSFdw);

      int nsubBjet = 0;

      for( int is = 0; is < FATnSubSDJet[ij]; ++is ){

	if( FATsubjetSDCSV[ij][is] > 0.605 ) ++nsubBjet;

      }
 
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
  
  return (float)passEvent/(float)data.GetEntriesFast();

}