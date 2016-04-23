#include <vector>
#include <string>
#include <iostream>
#include <TH1D.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "../../untuplizer.h"
#include "../../isPassZmumu.h"

void muZVariable(std::string inputFile, std::string outputFile){

  // read the ntuples (in pcncu)

  TreeReader data(inputFile.data());
  
  // Declare the histogram

  TFile* f = new TFile(inputFile.data());
  TH1D* h_totalEvents = (TH1D*)f->Get("h_totalEv");
     
  TH1D* h_Zmass        = new TH1D("h_Zmass",        "Zmass",        30, 60,  120);
  TH1D* h_Zpt          = new TH1D("h_Zpt",          "Zpt",          50,  0, 1000);
  TH1D* h_Zeta         = new TH1D("h_Zeta",         "Zeta",         40, -4,    4);
  TH1D* h_ZRapidity    = new TH1D("h_ZRapidity",    "ZRapidity",    40, -4,    4);
  TH1D* h_leadMuPt     = new TH1D("h_leadMuPt",     "leadMuPt",     50,  0, 1000);
  TH1D* h_leadMuEta    = new TH1D("h_leadMuEta",    "leadMuEta",    40, -4,    4);
  TH1D* h_subleadMuPt  = new TH1D("h_subleadMuPt",  "subleadMuPt",  25,  0,  500);
  TH1D* h_subleadMuEta = new TH1D("h_subleadMuEta", "subleadMuEta", 40, -4,    4);

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
    
  // begin of event loop

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    if( (unsigned)ev % 500000 == 0 )
      fprintf(stdout, "Still left events %lli of %lli\n", ev, data.GetEntriesFast());

    data.GetEntry(ev);

    Float_t eventWeight = data.GetFloat("ev_weight");
    TClonesArray* muP4  = (TClonesArray*) data.GetPtrTObject("muP4");

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

  TFile* outFile = new TFile(Form("%s_ZmumuVariable.root",outputFile.c_str()), "recreate");
      
  h_Zmass       ->Write("Zmass");
  h_Zpt         ->Write("Zpt");
  h_Zeta        ->Write("Zeta");
  h_ZRapidity   ->Write("ZRapidity");
  h_leadMuPt    ->Write("leadMuPt");
  h_leadMuEta   ->Write("leadMuEta");
  h_subleadMuPt ->Write("subleadMuPt");
  h_subleadMuEta->Write("subleadMuEta");
  h_totalEvents ->Write("totalEvents");
  
  outFile->Write();

  delete f;
  delete outFile;
  
}
