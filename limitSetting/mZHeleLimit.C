#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <TH1D.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "../untuplizer.h"
#include "../isPassZee.h"

void mZHeleLimit(string inputFile, string outputFile, Int_t btag=1){

  // read the ntuples (in pcncu)
  
  TreeReader data(inputFile.data());
  
  TFile* f = new TFile(inputFile.data());
  TH1D* h_totalEvents = (TH1D*)f->Get("h_totalEv");

  // Declare the histogram
     
  TH1D* h_mZprime = new TH1D("h_mZprime", "mZprime", 100, 400, 5000);

  h_mZprime->Sumw2();
  h_mZprime->GetXaxis()->SetTitle("mZprime");

  // begin of event loop

  fprintf(stdout, "Total events %lli\n", data.GetEntriesFast());

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    if( (unsigned)ev % 100000 == 0 )
      fprintf(stdout, "Still left events %lli\n", ev);

    data.GetEntry(ev);

    Float_t        eventWeight       = data.GetFloat("ev_weight");
    TClonesArray*  eleP4             = (TClonesArray*) data.GetPtrTObject("eleP4");
    Int_t          FATnJet           = data.GetInt("FATnJet");    
    Int_t*         FATnSubSDJet      = data.GetPtrInt("FATnSubSDJet");
    Float_t*       corrPRmass        = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    TClonesArray*  FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    vector<bool>&  FATjetPassIDLoose = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    vector<float>* FATsubjetSDCSV    = data.GetPtrVectorFloat("FATsubjetSDCSV", FATnJet);

    // select good leptons
      
    vector<Int_t> goodLepID;
    if( !isPassZee(data, goodLepID) ) continue;

    TLorentzVector* thisLep = (TLorentzVector*)eleP4->At(goodLepID[0]);
    TLorentzVector* thatLep = (TLorentzVector*)eleP4->At(goodLepID[1]);

    // select good FATjet

    Int_t goodFATJetID = -1;
    TLorentzVector* thisJet = NULL;

    for(Int_t ij = FATnJet-1; ij >= 0; --ij){

      thisJet = (TLorentzVector*)FATjetP4->At(ij);

      if( thisJet->Pt() < 200 ) continue;
      if( fabs(thisJet->Eta()) > 2.4 ) continue;
      if( corrPRmass[ij] < 105 || corrPRmass[ij] > 135 ) continue;  // choose events in signal region
      if( !FATjetPassIDLoose[ij] ) continue;
      if( thisJet->DeltaR(*thisLep) < 0.8 || thisJet->DeltaR(*thatLep) < 0.8 ) continue;

      Int_t nsubBjet = 0;

      for(Int_t is = FATnSubSDJet[ij]-1; is >= 0; --is){

	if( FATsubjetSDCSV[ij][is] > 0.605 ) ++nsubBjet;

      }

      // b-tag cut

      if( nsubBjet != btag ) continue;

      goodFATJetID = ij;
      break;

    }

    if( goodFATJetID < 0 ) continue;
    
    Float_t mllbb = (*thisLep+*thatLep+*thisJet).M();

    if( mllbb < 750 ) continue;

    h_mZprime->Fill(mllbb,eventWeight);

  } // end of event loop

  fprintf(stderr, "Processed all events\n");

  TFile* outFile = new TFile(Form("%s_mZHLimit.root",outputFile.c_str()), "recreate");

  h_mZprime    ->Write("mZprime");
  h_totalEvents->Write("totalEvents");

  outFile->Write();

  delete f;
  delete outFile;
  
}
