#include <vector>
#include <string>
#include <iostream>
#include <TH1D.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "../untuplizer.h"
#include "../readSample.h"
#include "../dataFilter.h"
#include "../isPassZmumu.h"
#include "../isPassZee.h"
#include "../correctMCweight.h"

void pseudoMC(std::string inputFile, std::string outputFile){

  // read the ntuples (in pcncu)

  std::vector<string> infiles;
    
  readSample(inputFile, infiles);
  
  TreeReader data(infiles);

  // Declare the histogram

  const Double_t xmin  = 500;
  const Double_t xmax  = 5000;
  const Int_t    nBins = (xmax-xmin)/100;
     
  TH1D* h_ZprimeSign_pMC   = new TH1D("h_ZprimeSign_pMC",   "ZprimeSign_pMC",  nBins, xmin, xmax);
  TH1D* h_ZprimeSide_pMC   = new TH1D("h_ZprimeSide_pMC",   "ZprimeSide_pMC",  nBins, xmin, xmax);
  TH1D* h_ZprimeSign_pDA   = new TH1D("h_ZprimeSign_pDA",   "ZprimeSign_pDA",  nBins, xmin, xmax);
  TH1D* h_ZprimeSide_pDA   = new TH1D("h_ZprimeSide_pDA",   "ZprimeSide_pDA",  nBins, xmin, xmax);
  TH1D* h_PRmassNoSIG_pDA  = new TH1D("h_PRmassNoSIG_pDA",  "PRmassNoSIG_pDA",    40,   40,  240);
  TH1D* h_PRmassAllLow_pDA = new TH1D("h_PRmassAllLow_pDA", "PRmassAllLow_pDA",   48,    0,  240);
  TH1D* h_PRmassAll_pDA    = new TH1D("h_PRmassAll_pDA",    "PRmassAll_pDA",      40,   40,  240);
  TH1D* h_eventWeight_pMC  = new TH1D("h_eventWeight_pMC",  "eventWeight_pMC",     2,   -1,    1);
  TH1D* h_eventWeight_pDA  = new TH1D("h_eventWeight_pDA",  "eventWeight_pDA",     2,   -1,    1);

  h_ZprimeSign_pMC   ->Sumw2();
  h_ZprimeSide_pMC   ->Sumw2();
  h_ZprimeSign_pDA   ->Sumw2();
  h_ZprimeSide_pDA   ->Sumw2();
  h_PRmassNoSIG_pDA  ->Sumw2();
  h_PRmassAllLow_pDA ->Sumw2();
  h_PRmassAll_pDA    ->Sumw2();

  h_ZprimeSign_pMC   ->GetXaxis()->SetTitle("ZprimeSign");
  h_ZprimeSide_pMC   ->GetXaxis()->SetTitle("ZprimeSide");
  h_ZprimeSign_pDA   ->GetXaxis()->SetTitle("ZprimeSign");
  h_ZprimeSide_pDA   ->GetXaxis()->SetTitle("ZprimeSide");
  h_PRmassNoSIG_pDA  ->GetXaxis()->SetTitle("corrPRmass");
  h_PRmassAllLow_pDA ->GetXaxis()->SetTitle("corrPRmass");
  h_PRmassAll_pDA    ->GetXaxis()->SetTitle("corrPRmass");

  // begin of event loop

  for(Long64_t ev = 0; ev < data.GetEntriesFast(); ev++){

    if( ev % 1000000 == 0 )
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());

    data.GetEntry(ev);

    Int_t          nVtx              = data.GetInt("nVtx");
    Bool_t         isData            = data.GetBool("isData");
    Float_t        pu_nTrueInt       = data.GetFloat("pu_nTrueInt");
    TClonesArray*  muP4              = (TClonesArray*) data.GetPtrTObject("muP4");
    TClonesArray*  eleP4             = (TClonesArray*) data.GetPtrTObject("eleP4");
    Int_t          FATnJet           = data.GetInt("FATnJet");    
    Float_t*       corrPRmass        = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    TClonesArray*  FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    vector<bool>&  FATjetPassIDLoose = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));

    // remove event which is no hard interaction (noise)

    if( nVtx < 1 ) continue;

    // Correct the pile-up shape of MC

    Double_t eventWeight = correctMCWeight(isData, (Int_t)pu_nTrueInt);

    if( ev % 2 == 0 )
      h_eventWeight_pMC->Fill(0.,eventWeight);

    else if( ev % 2 == 1 )
      h_eventWeight_pDA->Fill(0.,eventWeight);


    // data trigger cut
      
    bool muTrigger  = TriggerStatus(data, "HLT_Mu45");
    bool eleTrigger = TriggerStatus(data, "HLT_Ele105");

    // select good leptons
      
    vector<Int_t> goodLepID;

    TLorentzVector* thisLep = NULL;
    TLorentzVector* thatLep = NULL;

    if( !isPassZee(data,goodLepID) && !isPassZmumu(data,goodLepID) ) continue;

    if( isPassZee(data,goodLepID) ){
	
      if( !eleTrigger ) continue;

      thisLep = (TLorentzVector*)eleP4->At(goodLepID[0]);   
      thatLep = (TLorentzVector*)eleP4->At(goodLepID[1]);   

    }
    
    else if( isPassZmumu(data,goodLepID) ){
      
      if( !muTrigger ) continue;

      thisLep =  (TLorentzVector*)muP4->At(goodLepID[0]);   
      thatLep =  (TLorentzVector*)muP4->At(goodLepID[1]);   
      
    }

    // select good FATjet

    Int_t goodFATJetID = -1;
    TLorentzVector* thisJet = NULL;

    for(Int_t ij = 0; ij < FATnJet; ij++){

      thisJet = (TLorentzVector*)FATjetP4->At(ij);

      if( thisJet->Pt() < 200 ) continue;
      if( fabs(thisJet->Eta()) > 2.4 ) continue;
      if( !FATjetPassIDLoose[ij] ) continue;
      if( thisJet->DeltaR(*thisLep) < 0.8 || thisJet->DeltaR(*thatLep) < 0.8 ) continue;

      goodFATJetID = ij;
      break;

    }

    if( goodFATJetID < 0 ) continue;

    Float_t mllbb = (*thisLep+*thatLep+*thisJet).M();  

    if( mllbb < 500 ) continue;

    if( ev % 2 == 0 ){

      if( corrPRmass[goodFATJetID] > 40 && !(corrPRmass[goodFATJetID] > 65 && corrPRmass[goodFATJetID] < 145) )
	h_ZprimeSide_pMC->Fill(mllbb,eventWeight);

      if( corrPRmass[goodFATJetID] > 105 && corrPRmass[goodFATJetID] < 135 )
	h_ZprimeSign_pMC->Fill(mllbb,eventWeight);
  
    }

    else if( ev % 2 == 1 ){

      h_PRmassAll_pDA   ->Fill(corrPRmass[goodFATJetID],eventWeight);
      h_PRmassAllLow_pDA->Fill(corrPRmass[goodFATJetID],eventWeight);

      if( corrPRmass[goodFATJetID] > 40 && !(corrPRmass[goodFATJetID] > 65 && corrPRmass[goodFATJetID] < 145) ){
        h_ZprimeSide_pDA ->Fill(mllbb,eventWeight);
        h_PRmassNoSIG_pDA->Fill(corrPRmass[goodFATJetID],eventWeight);
      }

      if( corrPRmass[goodFATJetID] > 105 && corrPRmass[goodFATJetID] < 135 )
	h_ZprimeSign_pDA->Fill(mllbb,eventWeight);

    }

  } // end of event loop

  fprintf(stderr, "Processed all events\n");

  TFile* outFile = new TFile(Form("%s_pseudoTest.root",outputFile.c_str()), "recreate");
  
  h_ZprimeSign_pMC   ->Write("ZprimeSign_pMC");
  h_ZprimeSide_pMC   ->Write("ZprimeSide_pMC");
  h_ZprimeSign_pDA   ->Write("ZprimeSign_pDA");
  h_ZprimeSide_pDA   ->Write("ZprimeSide_pDA");
  h_PRmassNoSIG_pDA  ->Write("corrPRmass_pDA");
  h_PRmassAllLow_pDA ->Write("corrPRmassAllLow_pDA");
  h_PRmassAll_pDA    ->Write("corrPRmassAll_pDA");
  h_eventWeight_pMC  ->Write("eventWeight_pMC");
  h_eventWeight_pDA  ->Write("eventWeight_pDA");

  outFile->Write();

}
