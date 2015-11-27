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
#include "../../isPassZmumu.h"

void muZVariable(std::string inputFile, std::string outputFile){

  // read the ntuples (in pcncu)

  std::vector<string> infiles;

  readSample(inputFile, infiles);
  
  TreeReader data(infiles);
  
  // Declare the histogram
     
  Int_t nBin = 20;

  TH1D* h_Zmass        = new TH1D("h_Zmass", "Zmass", nBin, 60, 120);
  TH1D* h_Zpt          = new TH1D("h_Zpt", "Zpt", nBin, 0, 1000);
  TH1D* h_Zeta         = new TH1D("h_Zeta", "Zeta", nBin, -4, 4);
  TH1D* h_ZRapidity    = new TH1D("h_ZRapidity", "ZRapidity", nBin, -4, 4);
  TH1D* h_leadMuPt     = new TH1D("h_leadMuPt", "leadMuPt", nBin, 0, 1000);
  TH1D* h_leadMuEta    = new TH1D("h_leadMuEta", "leadMuEta", nBin, -4, 4);
  TH1D* h_subleadMuPt  = new TH1D("h_subleadMuPt", "subleadMuPt", nBin, 0, 500);
  TH1D* h_subleadMuEta = new TH1D("h_subleadMuEta", "subleadMuEta", nBin, -4, 4);
  TH1D* h_eventWeight  = new TH1D("h_eventWeight", "eventWeight", 100, -1, 1);

  h_Zmass       ->Sumw2();
  h_Zpt         ->Sumw2();
  h_Zeta        ->Sumw2();
  h_ZRapidity   ->Sumw2();
  h_leadMuPt    ->Sumw2();
  h_leadMuEta   ->Sumw2();
  h_subleadMuPt ->Sumw2();
  h_subleadMuEta->Sumw2(); 
  
  h_Zmass       ->GetXaxis()->SetTitle("Zmass"); 
  h_Zpt         ->GetXaxis()->SetTitle("Zpt");   
  h_Zeta        ->GetXaxis()->SetTitle("Zeta");    
  h_ZRapidity   ->GetXaxis()->SetTitle("ZRapidity");
  h_leadMuPt    ->GetXaxis()->SetTitle("leadMuPt");  
  h_leadMuEta   ->GetXaxis()->SetTitle("leadMuEta");
  h_subleadMuPt ->GetXaxis()->SetTitle("subleadMuPt");   
  h_subleadMuEta->GetXaxis()->SetTitle("subleadMuEta"); 
  h_eventWeight ->GetXaxis()->SetTitle("eventWeight");
    
  // begin of event loop

  for( Long64_t ev = 0; ev < data.GetEntriesFast(); ev++ ){

    if( ev % 1000000 == 0 )
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());

    data.GetEntry(ev);

    Int_t    nVtx      = data.GetInt("nVtx");
    Float_t  mcWeight  = data.GetFloat("mcWeight");    
    TClonesArray* muP4 = (TClonesArray*) data.GetPtrTObject("muP4");

    Double_t eventWeight = mcWeight;
    if( inputFile.find("DYJets") != std::string::npos ){
      if( eventWeight > 0 ) eventWeight = 1;
      else if( eventWeight < 0 ) eventWeight = -1;
    }
    else
      eventWeight = 1;
    
    h_eventWeight->Fill(0.,eventWeight);

    if( nVtx < 1 ) continue;
    
    // data trigger cut (muon channel)

    std::string* trigName = data.GetPtrString("hlt_trigName");
    vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
    const Int_t nsize = data.GetPtrStringSize();    
    bool passTrigger = false;
    
    for(Int_t it = 0; it < nsize; it++){
    
      std::string thisTrig = trigName[it];
      bool results = trigResult[it];
      
      if( thisTrig.find("HLT_Mu45") != std::string::npos && results==1 ){
	passTrigger = true;
	break;
      }
      
    }

    if( !passTrigger ) continue;

    // select good muons
      
    vector<Int_t> goodMuID;
    if( !isPassZmumu(data, goodMuID) ) continue;

    TLorentzVector* thisMu = (TLorentzVector*)muP4->At(goodMuID[0]);
    TLorentzVector* thatMu = (TLorentzVector*)muP4->At(goodMuID[1]);

    TLorentzVector l4_Z = (*thisMu+*thatMu);

    h_Zmass    ->Fill(l4_Z.M(),eventWeight);
    h_Zpt      ->Fill(l4_Z.Pt(),eventWeight);
    h_Zeta     ->Fill(l4_Z.Eta(),eventWeight);
    h_ZRapidity->Fill(l4_Z.Rapidity(),eventWeight);

    if( thisMu->Pt() > thatMu->Pt() ){

      h_leadMuPt    ->Fill(thisMu->Pt(),eventWeight);
      h_leadMuEta   ->Fill(thisMu->Eta(),eventWeight);
      h_subleadMuPt ->Fill(thatMu->Pt(),eventWeight);
      h_subleadMuEta->Fill(thatMu->Eta(),eventWeight);

    }else{

      h_leadMuPt    ->Fill(thatMu->Pt(),eventWeight);
      h_leadMuEta   ->Fill(thatMu->Eta(),eventWeight);
      h_subleadMuPt ->Fill(thisMu->Pt(),eventWeight);
      h_subleadMuEta->Fill(thisMu->Eta(),eventWeight);

    }

  } // end of event loop

  fprintf(stderr, "Processed all events\n");

  std::string h_name[9] = {"Zmass","Zpt","Zeta","ZRapidity","leadMuPt","leadMuEta",
			   "subleadMuPt","subleadMuEta","eventWeight"};

  TFile* outFile = new TFile(Form("%s_ZmumuVariable.root",outputFile.c_str()), "recreate");
      
  h_Zmass       ->Write(h_name[0].data());  
  h_Zpt         ->Write(h_name[1].data());  
  h_Zeta        ->Write(h_name[2].data());    
  h_ZRapidity   ->Write(h_name[3].data());
  h_leadMuPt    ->Write(h_name[4].data()); 
  h_leadMuEta   ->Write(h_name[5].data());   
  h_subleadMuPt ->Write(h_name[6].data());
  h_subleadMuEta->Write(h_name[7].data());    
  h_eventWeight ->Write(h_name[8].data());
  
  outFile->Write();
  
}
