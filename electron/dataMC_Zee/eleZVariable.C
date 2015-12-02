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
#include "../../isPassZee.h"
#include "../../correctMCweight.h"

void eleZVariable(std::string inputFile, std::string outputFile){

  // read the ntuples (in pcncu)

  std::vector<string> infiles;
 
  readSample(inputFile, infiles);
  
  TreeReader data(infiles);
  
  // Declare the histogram
     
  TH1D* h_Zmass         = new TH1D("h_Zmass",         "Zmass",         30, 60,  120);
  TH1D* h_Zpt           = new TH1D("h_Zpt",           "Zpt",           50,  0, 1000);
  TH1D* h_Zeta          = new TH1D("h_Zeta",          "Zeta",          40, -4,    4);
  TH1D* h_ZRapidity     = new TH1D("h_ZRapidity",     "ZRapidity",     40, -4,    4);
  TH1D* h_leadElePt     = new TH1D("h_leadElePt",     "leadElePt",     16,  0,  800);
  TH1D* h_leadEleEta    = new TH1D("h_leadEleEta",    "leadEleEta",    40, -4,    4);
  TH1D* h_subleadElePt  = new TH1D("h_subleadElePt",  "subleadElePt",  25,  0,  500);
  TH1D* h_subleadEleEta = new TH1D("h_subleadEleEta", "subleadEleEta", 40, -4,    4);
  TH1D* h_eventWeight   = new TH1D("h_eventWeight",   "eventWeight",    2, -1,    1);

  h_Zmass        ->Sumw2();
  h_Zpt          ->Sumw2();
  h_Zeta         ->Sumw2();
  h_ZRapidity    ->Sumw2();
  h_leadElePt    ->Sumw2();
  h_leadEleEta   ->Sumw2();
  h_subleadElePt ->Sumw2();
  h_subleadEleEta->Sumw2(); 
  
  h_Zmass        ->GetXaxis()->SetTitle("Zmass"); 
  h_Zpt          ->GetXaxis()->SetTitle("Zpt");   
  h_Zeta         ->GetXaxis()->SetTitle("Zeta");    
  h_ZRapidity    ->GetXaxis()->SetTitle("ZRapidity");
  h_leadElePt    ->GetXaxis()->SetTitle("leadElePt");  
  h_leadEleEta   ->GetXaxis()->SetTitle("leadEleEta");
  h_subleadElePt ->GetXaxis()->SetTitle("subleadElePt");   
  h_subleadEleEta->GetXaxis()->SetTitle("subleadEleEta"); 
    
  // begin of event loop

  for( Long64_t ev = 0; ev < data.GetEntriesFast(); ev++ ){

    if( ev % 1000000 == 0 )
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());

    data.GetEntry(ev);

    Int_t   nVtx        = data.GetInt("nVtx");
    Bool_t  isData      = data.GetBool("isData");
    Float_t pu_nTrueInt = data.GetFloat("pu_nTrueInt");
    TClonesArray* eleP4 = (TClonesArray*) data.GetPtrTObject("eleP4");

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

    // select good electrons

    vector<Int_t> goodEleID;
    if( !isPassZee(data, goodEleID) ) continue;

    TLorentzVector* thisEle = (TLorentzVector*)eleP4->At(goodEleID[0]);
    TLorentzVector* thatEle = (TLorentzVector*)eleP4->At(goodEleID[1]);

    TLorentzVector l4_Z = (*thisEle+*thatEle);

    h_Zmass    ->Fill(l4_Z.M(),eventWeight);
    h_Zpt      ->Fill(l4_Z.Pt(),eventWeight);
    h_Zeta     ->Fill(l4_Z.Eta(),eventWeight);
    h_ZRapidity->Fill(l4_Z.Rapidity(),eventWeight);

    if( thisEle->Pt() > thatEle->Pt() ){

      h_leadElePt    ->Fill(thisEle->Pt(),eventWeight);
      h_leadEleEta   ->Fill(thisEle->Eta(),eventWeight);
      h_subleadElePt ->Fill(thatEle->Pt(),eventWeight);
      h_subleadEleEta->Fill(thatEle->Eta(),eventWeight);

    }else{

      h_leadElePt    ->Fill(thatEle->Pt(),eventWeight);
      h_leadEleEta   ->Fill(thatEle->Eta(),eventWeight);
      h_subleadElePt ->Fill(thisEle->Pt(),eventWeight);
      h_subleadEleEta->Fill(thisEle->Eta(),eventWeight);

    }

  } // end of event loop

  fprintf(stderr, "Processed all events\n");

  TFile* outFile = new TFile(Form("%s_ZeeVariable.root",outputFile.c_str()), "recreate");
      
  h_Zmass        ->Write("Zmass");
  h_Zpt          ->Write("Zpt");
  h_Zeta         ->Write("Zeta");
  h_ZRapidity    ->Write("ZRapidity");
  h_leadElePt    ->Write("leadElePt");
  h_leadEleEta   ->Write("leadEleEta");
  h_subleadElePt ->Write("subleadElePt");
  h_subleadEleEta->Write("subleadEleEta");
  h_eventWeight  ->Write("eventWeight");
  
  outFile->Write();
  
}
