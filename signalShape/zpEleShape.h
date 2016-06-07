#include <vector>
#include <string>
#include <iostream>
#include <TH1F.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "../untuplizer.h"
#include "../isPassZee.h"

Float_t CrossSection(string token){

  ifstream textFile("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/xSec.txt");
  string thisSample;

  Float_t crosssection = 0., thisNum = 0.;

  while( textFile >> thisSample >> thisNum ){

    if( token.find(thisSample) != string::npos )
      crosssection = thisNum;

  }

  return crosssection;
  
}

TH1F* zpEleShape(string inputFile, int cat){

  // read the ntuples (in pcncu)

  TreeReader data(inputFile.data());
  
  TFile* f = new TFile(inputFile.data());
  TH1F* h_totalEvents = (TH1F*)f->Get("h_totalEv");

  // Declare the histogram
     
  TH1F* h_mZprime = new TH1F("h_mZprime", "mZprime", 100, 400, 5000);

  h_mZprime->Sumw2();
  h_mZprime->GetXaxis()->SetTitle("mZprime");

  // Calculate the scale correspond to inputFile

  Float_t scale = 2512./(h_totalEvents->Integral()/CrossSection(inputFile.data()));

  // begin of event loop

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    data.GetEntry(ev);

    Float_t        eventWeight       = data.GetFloat("ev_weight");
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

    h_mZprime->Fill((*thisLep+*thatLep+thisJet).M(),eventWeight*scale);

  } // end of event loop
  
  return h_mZprime;

  delete f;

}