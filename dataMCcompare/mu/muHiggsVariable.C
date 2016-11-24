#include <vector>
#include <string>
#include <iostream>
#include <TH1D.h>
#include <TMath.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/untuplizer.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZmumu.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassJet.h"

void muHiggsVariable(string inputFile, string outputFile){

  // read the ntuples (in pcncu)

  TreeReader data(inputFile.data());
  
  // Declare the histogram
     
  TFile f(inputFile.data());
  TH1D* h_totalEvents = (TH1D*)f.Get("h_totalEv");

  TH1D* h_FATjetPRmassCorr = new TH1D("h_FATjetPRmassCorr", "HiggsMass", 20,    0,  200);

  h_FATjetPRmassCorr->Sumw2();
  h_FATjetPRmassCorr->GetXaxis()->SetTitle("HiggsMass");
    
  // begin of event loop

  fprintf(stdout, "Total events %lli\n", data.GetEntriesFast());

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    if( (unsigned)ev % 100000 == 0 )
      fprintf(stdout, "Still left events %lli\n", ev);

    data.GetEntry(ev);

    Float_t        eventWeight       = data.GetFloat("ev_weight");
    TClonesArray*  muP4             = (TClonesArray*) data.GetPtrTObject("muP4");
    Float_t*       FATjetPRmassCorr  = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    TClonesArray*  FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");

    // select good muons

    vector<Int_t> goodLepID;

    if( !isPassZmumu(data,goodLepID) ) continue;

    TLorentzVector* thisLep = (TLorentzVector*)muP4->At(goodLepID[0]);
    TLorentzVector* thatLep = (TLorentzVector*)muP4->At(goodLepID[1]);

    // select good FATjet (include signal region and side bands)

    Int_t goodFATJetID = -1;

    if( !isPassJet(data, &goodFATJetID, thisLep, thatLep, false, false, true) ) continue;

    TLorentzVector* thisJet = (TLorentzVector*)FATjetP4->At(goodFATJetID);

    Float_t mllbb;

    if( !noiseCleaning(thisLep, thatLep, thisJet, &mllbb) ) continue;

    h_FATjetPRmassCorr->Fill(FATjetPRmassCorr[goodFATJetID],eventWeight);

  } // end of event loop

  fprintf(stdout, "Processed all events\n");

  TFile* outFile = new TFile(Form("%s_muHiggsVariable.root",outputFile.data()), "recreate");
      
  h_FATjetPRmassCorr->Write("HiggsMass");
  h_totalEvents     ->Write("totalEvents");
  
  outFile->Write();

}
