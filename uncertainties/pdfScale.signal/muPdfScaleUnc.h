#include <vector>
#include <string>
#include <iostream>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/untuplizer.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZmumu.h"

float muPdfScaleUnc(string inputFile, int cat, int first, int last){

  // read the ntuples (in pcncu)

  TreeReader data(inputFile.data());

  int N = 1+last-first;
  float* efficiency = new float[N];
  float* passEvent  = new float[N];
  float  totalEvent = 1./data.GetEntriesFast();
  float  cpass = 0.;

  std::fill_n(efficiency,N,0.);
  std::fill_n(passEvent,N,0.);

  // begin of event loop

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    data.GetEntry(ev);

    Float_t*       pdfscaleSysWeight = data.GetPtrFloat("pdfscaleSysWeights");
    TClonesArray*  muP4              = (TClonesArray*) data.GetPtrTObject("muP4");
    Int_t          FATnJet           = data.GetInt("FATnJet");    
    Int_t*         FATnSubSDJet      = data.GetPtrInt("FATnSubSDJet");
    Float_t*       FATjetPRmassCorr  = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    TClonesArray*  FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    vector<bool>&  FATjetPassIDLoose = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    vector<float>* FATsubjetSDCSV    = data.GetPtrVectorFloat("FATsubjetSDCSV", FATnJet);

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

      goodFATJetID = ij;
      thisJet = *myJet;

      break;
 
    } // end of FatnJet loop
 
    if( goodFATJetID < 0 ) continue;

    if( (*thisLep+*thatLep+thisJet).M() < 750 ) continue;

    int nsubBjet = 0;

    for( int is = 0; is < FATnSubSDJet[goodFATJetID]; ++is ){
            
      if( FATsubjetSDCSV[goodFATJetID][is] > 0.605 ) ++nsubBjet;
      
    } // end of subjet loop
    
    // b-tag cut
    
    if( cat == 1 && nsubBjet != 1 ) continue;
    if( cat == 2 && nsubBjet != 2 ) continue;
        
    ++cpass;

    int i = 0;
    for( int n = first; n <= last; ++n ){
      passEvent[i] += pdfscaleSysWeight[n];
      ++i;
    }

  } // end of event loop

  for( int n = 0; n < N; ++n )
    efficiency[n] = passEvent[n]*totalEvent;

  return (first != 0) ?
    TMath::RMS(N, efficiency)/(cpass*totalEvent) : 
    TMath::Max(fabs(efficiency[2]-efficiency[0]), fabs(efficiency[1]-efficiency[0]));

  delete [] passEvent;
  delete [] efficiency;

}
