#include <vector>
#include <string>
#include <iostream>
#include <TH1D.h>
#include <TMath.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "../../untuplizer.h"
#include "../../readSample.h"
#include "../../dataFilter.h"
#include "../../isPassZmumu.h"
#include "../../correctMCweight.h"

void muAlphaRatio(std::string inputFile, std::string outputFile){

  // read the ntuples (in pcncu)

  std::vector<string> infiles;
    
  readSample(inputFile, infiles);
  
  TreeReader data(infiles);

  // Declare the histogram

  const Double_t varBins[] = {500,700,900,1100,1300,1500,1700,1900,2100,??????};
  Int_t nvarBins = sizeof(varBins)/sizeof(varBins[0])-1;
     
  TH1D* h_ZprimeSign    = new TH1D("h_ZprimeSign",    "ZprimeSign", nvarBins, varBins);
  TH1D* h_ZprimeSide    = new TH1D("h_ZprimeSide",    "ZprimeSide", nvarBins, varBins);
  TH1D* h_corrPRmass    = new TH1D("h_corrPRmass",    "corrPRmass",       40, 40, 240);
  TH1D* h_corrPRmassAll = new TH1D("h_corrPRmassAll", "corrPRmassAll",    48,  0, 240);
  TH1D* h_eventWeight   = new TH1D("h_eventWeight",   "eventWeight",       2, -1,   1);

  h_ZprimeSign    ->Sumw2();
  h_ZprimeSide    ->Sumw2();
  h_corrPRmass    ->Sumw2();
  h_corrPRmassAll ->Sumw2();

  h_ZprimeSign    ->GetXaxis()->SetTitle("ZprimeSign");
  h_ZprimeSide    ->GetXaxis()->SetTitle("ZprimeSide");
  h_corrPRmass    ->GetXaxis()->SetTitle("corrPRmass");
  h_corrPRmassAll ->GetXaxis()->SetTitle("corrPRmass");

  // begin of event loop

  for( Long64_t ev = 0; ev < data.GetEntriesFast(); ev++ ){

    if( ev % 1000000 == 0 )
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());

    data.GetEntry(ev);

    Int_t          nVtx              = data.GetInt("nVtx");
    Bool_t         isData            = data.GetBool("isData");
    TClonesArray*  muP4              = (TClonesArray*) data.GetPtrTObject("muP4");
    Int_t          FATnJet           = data.GetInt("FATnJet");    
    Int_t*         FATnSubSDJet      = data.GetPtrInt("FATnSubSDJet");
    Float_t*       corrPRmass        = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    TClonesArray*  FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    vector<bool>&  FATjetPassIDLoose = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
  
    // remove event which is no hard interaction (noise)

    if( nVtx < 1 ) continue;

    // Correct the pile-up shape of MC

    Double_t eventWeight = correctMCWeight(isData, nVtx);
    
    h_eventWeight->Fill(0.,eventWeight);

    // data filter and trigger cut
      
    bool muTrigger = TriggerStatus(data, "HLT_Mu45");
    bool CSCT      = FilterStatus(data, "Flag_CSCTightHaloFilter");
    bool eeBadSc   = FilterStatus(data, "Flag_eeBadScFilter");
    bool Noise     = FilterStatus(data, "Flag_HBHENoiseFilter");
    bool NoiseIso  = FilterStatus(data, "Flag_HBHENoiseIsoFilter");

    if( !muTrigger ) continue;
    if( isData && !CSCT ) continue;
    if( isData && !eeBadSc ) continue;
    if( isData && !Noise ) continue;
    if( isData && !NoiseIso ) continue;

    // select good muons
      
    vector<Int_t> goodMuID;
    if( !isPassZmumu(data, goodMuID) ) continue;

    TLorentzVector* thisMu = (TLorentzVector*)muP4->At(goodMuID[0]);
    TLorentzVector* thatMu = (TLorentzVector*)muP4->At(goodMuID[1]);

    // select good FATjet

    Int_t goodFATJetID = -1;
    TLorentzVector* thisJet = NULL;

    for(Int_t ij = 0; ij < FATnJet; ij++){

      thisJet = (TLorentzVector*)FATjetP4->At(ij);

      if( thisJet->Pt() < 200 ) continue;
      if( fabs(thisJet->Eta()) > 2.4 ) continue;
      if( !FATjetPassIDLoose[ij] ) continue;
      if( thisJet->DeltaR(*thisMu) < 0.8 || thisJet->DeltaR(*thatMu) < 0.8 ) continue;

      goodFATJetID = ij;
      break;

    }

    if( goodFATJetID < 0 ) continue;

    Float_t mllbb = (*thisMu+*thatMu+*thisJet).M();  

    if( mllbb < 500 ) continue;

    h_corrPRmassAll->Fill(corrPRmass[goodFATJetID],eventWeight);

    if( corrPRmass[goodFATJetID] > 40 && !(corrPRmass[goodFATJetID] > 65 && corrPRmass[goodFATJetID] < 145) ){
      h_ZprimeSide->Fill(mllbb,eventWeight);
      h_corrPRmass->Fill(corrPRmass[goodFATJetID],eventWeight);
    }

    if( corrPRmass[goodFATJetID] > 105 && corrPRmass[goodFATJetID] < 135 )
      h_ZprimeSign->Fill(mllbb,eventWeight);

  } // end of event loop

  fprintf(stderr, "Processed all events\n");

  TFile* outFile = new TFile(Form("%s_muAlphaR.root",outputFile.c_str()), "recreate");
  
  h_ZprimeSign    ->Write("ZprimeSign");
  h_ZprimeSide    ->Write("ZprimeSide");
  h_corrPRmass    ->Write("corrPRmass");
  h_corrPRmassAll ->Write("corrPRmassAll");
  h_eventWeight   ->Write("eventWeight");

  outFile->Write();

}
