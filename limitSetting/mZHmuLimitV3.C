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
#include "../correctMCweight.h"

void mZHmuLimitV3(std::string inputFile, std::string outputFile){

  // read the ntuples (in pcncu)

  std::vector<string> infiles;

  readSample(inputFile, infiles);
  
  TreeReader data(infiles);
  
  // Declare the histogram
     
  TH1D* h_mZprime     = new TH1D("h_mZprime",     "mZprime",     100, 400, 5000);
  TH1D* h_eventWeight = new TH1D("h_eventWeight", "eventWeight",   2,  -1,    1);

  h_mZprime->Sumw2();
  h_mZprime->GetXaxis()->SetTitle("mZprime");
      
  // begin of event loop

  for( Long64_t ev = 0; ev < data.GetEntriesFast(); ev++){

    if( ev % 1000000 == 0 )
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());

    data.GetEntry(ev);

    Int_t          nVtx              = data.GetInt("nVtx");
    Bool_t         isData            = data.GetBool("isData");
    Float_t        pu_nTrueInt       = data.GetFloat("pu_nTrueInt");
    TClonesArray*  muP4              = (TClonesArray*) data.GetPtrTObject("muP4");
    Int_t          FATnJet           = data.GetInt("FATnJet");    
    Int_t*         FATnSubSDJet      = data.GetPtrInt("FATnSubSDJet");
    Float_t*       corrPRmass        = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    TClonesArray*  FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    vector<bool>&  FATjetPassIDLoose = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    vector<float>* FATsubjetSDPx     = data.GetPtrVectorFloat("FATsubjetSDPx", FATnJet);
    vector<float>* FATsubjetSDPy     = data.GetPtrVectorFloat("FATsubjetSDPy", FATnJet);
    vector<float>* FATsubjetSDPz     = data.GetPtrVectorFloat("FATsubjetSDPz", FATnJet);
    vector<float>* FATsubjetSDE      = data.GetPtrVectorFloat("FATsubjetSDE", FATnJet);
    vector<float>* FATsubjetSDCSV    = data.GetPtrVectorFloat("FATsubjetSDCSV", FATnJet);

    // remove event which is no hard interaction (noise)

    if( nVtx < 1 ) continue;

    // Correct the pile-up shape of MC

    Double_t eventWeight = correctMCWeight(isData, (Int_t)pu_nTrueInt);
    
    h_eventWeight->Fill(0.,eventWeight);

    // data filter (to filter non-collision bkg (ECAL/HCAL noise)) and trigger cut
      
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

    Int_t goodJetID = -1;
    TLorentzVector* thisJet = NULL;

    for(Int_t ij = 0; ij < FATnJet; ij++){

      thisJet = (TLorentzVector*)FATjetP4->At(ij);

      if( thisJet->Pt() < 200 ) continue;
      if( fabs(thisJet->Eta()) > 2.4 ) continue;
      if( corrPRmass[ij] < 105 || corrPRmass[ij] > 135 ) continue;
      if( !FATjetPassIDLoose[ij] ) continue;
      if( thisJet->DeltaR(*thisMu) < 0.8 || thisJet->DeltaR(*thatMu) < 0.8 ) continue;

      Int_t nsubBjet = 0;

      for(Int_t is = 0; is < FATnSubSDJet[ij]; is++){

	if( FATsubjetSDCSV[ij][is] > 0.605 ) nsubBjet++;

      }

      if( nsubBjet < 1 ) continue;

      goodJetID = ij;
      break;

    }

    if( goodJetID < 0 ) continue;
    
    Float_t mllbb = (*thisMu+*thatMu+*thisJet).M();

    if( mllbb < 700 ) continue;

    h_mZprime->Fill(mllbb,eventWeight);

  } // end of event loop

  fprintf(stderr, "Processed all events\n");

  TFile* outFile = new TFile(Form("%s_mZHmuSignal_loose.root",outputFile.c_str()), "recreate");

  h_mZprime    ->Write("mZprime");
  h_eventWeight->Write("eventWeight");

  outFile->Write();
  
}
