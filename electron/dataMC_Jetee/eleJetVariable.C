#include <vector>
#include <string>
#include <iostream>
#include <TH1D.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "../../untuplizer.h"
#include "../../readSample.h"
#include "../../dataFilter.h"
#include "../../isPassZee.h"
#include "../../correctMCweight.h"

void eleJetVariable(std::string inputFile, std::string outputFile){

  // read the ntuples (in pcncu)

  std::vector<string> infiles;
    
  readSample(inputFile, infiles);
  
  TreeReader data(infiles);
  
  // Declare the histogram

  Int_t nBin = 20;

  TH1D* h_nVtx             = new TH1D("h_nVtx",             "nVtx",              100,  -1,   99);     
  TH1D* h_FATjetPt         = new TH1D("h_FATjetPt",         "FATjetPt",         nBin, 100, 1000);
  TH1D* h_FATjetEta        = new TH1D("h_FATjetEta",        "FATjetEta",        nBin,  -4,    4);
  TH1D* h_FATjetCISVV2     = new TH1D("h_FATjetCISVV2",     "FATjetCISVV2",     nBin,   0,  1.2);
  TH1D* h_FATjetSDmass     = new TH1D("h_FATjetSDmass",     "FATjetSDmass",     nBin,   0,  200);
  TH1D* h_FATjetPRmass     = new TH1D("h_FATjetPRmass",     "FATjetPRmass",     nBin,   0,  200);
  TH1D* h_FATjetPRmassCorr = new TH1D("h_FATjetPRmassCorr", "FATjetPRmassCorr", nBin,   0,  200);
  TH1D* h_FATjetTau1       = new TH1D("h_FATjetTau1",       "FATjetTau1",       nBin,   0,    1);
  TH1D* h_FATjetTau2       = new TH1D("h_FATjetTau2",       "FATjetTau2",       nBin,   0,    1);
  TH1D* h_FATjetTau2dvTau1 = new TH1D("h_FATjetTau2dvTau1", "FATjetTau2dvTau1", nBin,   0,    1);
  TH1D* h_FATsubjetPt      = new TH1D("h_FATsubjetPt",      "FATsubjetPt",      nBin,   0,  800);
  TH1D* h_FATsubjetEta     = new TH1D("h_FATsubjetEta",     "FATsubjetEta",     nBin,  -4,    4);  
  TH1D* h_FATsubjetSDCSV   = new TH1D("h_FATsubjetSDCSV",   "FATsubjetSDCSV",   nBin,   0,  1.2);
  TH1D* h_eventWeight      = new TH1D("h_eventWeight",      "eventWeight",         2,  -1,    1);

  h_nVtx            ->Sumw2();
  h_FATjetPt        ->Sumw2();   
  h_FATjetEta       ->Sumw2();
  h_FATjetCISVV2    ->Sumw2();
  h_FATjetSDmass    ->Sumw2();
  h_FATjetPRmass    ->Sumw2();
  h_FATjetPRmassCorr->Sumw2();
  h_FATjetTau1      ->Sumw2();
  h_FATjetTau2      ->Sumw2();
  h_FATjetTau2dvTau1->Sumw2();
  h_FATsubjetPt     ->Sumw2();
  h_FATsubjetEta    ->Sumw2();
  h_FATsubjetSDCSV  ->Sumw2();

  h_nVtx            ->GetXaxis()->SetTitle("nVtx");
  h_FATjetPt        ->GetXaxis()->SetTitle("FATjetPt");
  h_FATjetEta       ->GetXaxis()->SetTitle("FATjetEta");
  h_FATjetCISVV2    ->GetXaxis()->SetTitle("FATjetCISVV2");
  h_FATjetSDmass    ->GetXaxis()->SetTitle("FATjetSDmass");
  h_FATjetPRmass    ->GetXaxis()->SetTitle("FATjetPRmass");
  h_FATjetPRmassCorr->GetXaxis()->SetTitle("FATjetPRmassL2L3Corr");
  h_FATjetTau1      ->GetXaxis()->SetTitle("FATjetTau1");
  h_FATjetTau2      ->GetXaxis()->SetTitle("FATjetTau2");
  h_FATjetTau2dvTau1->GetXaxis()->SetTitle("FATjetTau2dvTau1");
  h_FATsubjetPt     ->GetXaxis()->SetTitle("FATsubjetPt");
  h_FATsubjetEta    ->GetXaxis()->SetTitle("FATsubjetEta");
  h_FATsubjetSDCSV  ->GetXaxis()->SetTitle("FATsubjetSDCSV");
    
  // begin of event loop

  for( Long64_t ev = 0; ev < data.GetEntriesFast(); ev++ ){

    if( ev % 1000000 == 0 )
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());

    data.GetEntry(ev);

    Int_t          nVtx              = data.GetInt("nVtx");
    Bool_t         isData            = data.GetBool("isData");
    TClonesArray*  eleP4             = (TClonesArray*) data.GetPtrTObject("eleP4");
    Int_t          FATnJet           = data.GetInt("FATnJet");    
    Int_t*         FATnSubSDJet      = data.GetPtrInt("FATnSubSDJet");
    Float_t        pu_nTrueInt       = data.GetFloat("pu_nTrueInt");
    Float_t*       FATjetCISVV2      = data.GetPtrFloat("FATjetCISVV2");
    Float_t*       FATjetSDmass      = data.GetPtrFloat("FATjetSDmass");
    Float_t*       FATjetPRmass      = data.GetPtrFloat("FATjetPRmass");
    Float_t*       FATjetPRmassCorr  = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    Float_t*       FATjetTau1        = data.GetPtrFloat("FATjetTau1");
    Float_t*       FATjetTau2        = data.GetPtrFloat("FATjetTau2");
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

    // data filter and trigger cut
      
    bool eleTrigger = TriggerStatus(data, "HLT_Ele105");
    bool CSCT       = FilterStatus(data, "Flag_CSCTightHaloFilter");
    bool eeBadSc    = FilterStatus(data, "Flag_eeBadScFilter");
    bool Noise      = FilterStatus(data, "Flag_HBHENoiseFilter");
    bool NoiseIso   = FilterStatus(data, "Flag_HBHENoiseIsoFilter");

    if( !eleTrigger ) continue;
    if( isData && !CSCT ) continue;
    if( isData && !eeBadSc ) continue;
    if( isData && !Noise ) continue;
    if( isData && !NoiseIso ) continue;

    // select good Electrons

    vector<Int_t> goodEleID;
    if( !isPassZee(data, goodEleID) ) continue;

    TLorentzVector* thisEle = (TLorentzVector*)eleP4->At(goodEleID[0]);
    TLorentzVector* thatEle = (TLorentzVector*)eleP4->At(goodEleID[1]);

    // select good FATjet

    Int_t goodFATJetID = -1;
    TLorentzVector* thisJet = NULL;

    for(Int_t ij = 0; ij < FATnJet; ij++){

      thisJet = (TLorentzVector*)FATjetP4->At(ij);

      if( thisJet->Pt() < 200 ) continue;
      if( fabs(thisJet->Eta()) > 2.4 ) continue;
      if( !FATjetPassIDLoose[ij] ) continue;
      if( thisJet->DeltaR(*thisEle) < 0.8 || thisJet->DeltaR(*thatEle) < 0.8 ) continue;
      if( FATjetPRmassCorr[ij] > 65 && FATjetPRmassCorr[ij] < 145 ) continue;

      goodFATJetID = ij;
      break;

    }

    if( goodFATJetID < 0 ) continue;

    h_nVtx            ->Fill(nVtx,eventWeight);
    h_FATjetPt        ->Fill(thisJet->Pt(),eventWeight);
    h_FATjetEta       ->Fill(thisJet->Eta(),eventWeight);
    h_FATjetCISVV2    ->Fill(FATjetCISVV2[goodFATJetID],eventWeight);
    h_FATjetSDmass    ->Fill(FATjetSDmass[goodFATJetID],eventWeight);
    h_FATjetPRmass    ->Fill(FATjetPRmass[goodFATJetID],eventWeight);
    h_FATjetPRmassCorr->Fill(FATjetPRmassCorr[goodFATJetID],eventWeight);
    h_FATjetTau1      ->Fill(FATjetTau1[goodFATJetID],eventWeight);
    h_FATjetTau2      ->Fill(FATjetTau2[goodFATJetID],eventWeight);
    h_FATjetTau2dvTau1->Fill(FATjetTau2[goodFATJetID]/FATjetTau1[goodFATJetID],eventWeight);

    for(Int_t is = 0; is < FATnSubSDJet[goodFATJetID]; is++){

      TLorentzVector l4_subjet(0,0,0,0);

      l4_subjet.SetPxPyPzE(FATsubjetSDPx[goodFATJetID][is],
			   FATsubjetSDPy[goodFATJetID][is],
			   FATsubjetSDPz[goodFATJetID][is],
			   FATsubjetSDE [goodFATJetID][is]);

      h_FATsubjetPt   ->Fill(l4_subjet.Pt(),eventWeight);
      h_FATsubjetEta  ->Fill(l4_subjet.Eta(),eventWeight);
      h_FATsubjetSDCSV->Fill(FATsubjetSDCSV[goodFATJetID][is],eventWeight);

    }

  } // end of event loop

  fprintf(stderr, "Processed all events\n");

  TFile* outFile = new TFile(Form("%s_jeteeVariable.root",outputFile.c_str()), "recreate");

  h_nVtx            ->Write("nVtx");
  h_FATjetPt        ->Write("FATjetPt");
  h_FATjetEta       ->Write("FATjetEta");
  h_FATjetCISVV2    ->Write("FATjetCISVV2");
  h_FATjetSDmass    ->Write("FATjetSDmass");
  h_FATjetPRmass    ->Write("FATjetPRmass");
  h_FATjetPRmassCorr->Write("FATjetPRmassCorr");
  h_FATjetTau1      ->Write("FATjetTau1");
  h_FATjetTau2      ->Write("FATjetTau2");
  h_FATjetTau2dvTau1->Write("FATjetTau2dvTau1");
  h_FATsubjetPt     ->Write("FATsubjetPt");
  h_FATsubjetEta    ->Write("FATsubjetEta");
  h_FATsubjetSDCSV  ->Write("FATsubjetSDCSV");
  h_eventWeight     ->Write("eventWeight");
  
  outFile->Write();
  
}
