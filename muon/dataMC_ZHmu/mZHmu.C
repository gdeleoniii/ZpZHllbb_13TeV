#include <vector>
#include <string>
#include <iostream>
#include <TH1D.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "../../untuplizer.h"
#include "../../isPassZmumu.h"

void mZHmu(std::string inputFile, std::string outputFile){

  // read the ntuples (in pcncu)

  TreeReader data(inputFile.data());
  
  // Declare the histogram

  TFile* f = new TFile(inputFile.data());
  TH1D* h_totalEvents = (TH1D*)f->Get("h_totalEv");

  TH1D* h_mZprime          = new TH1D("h_mZprime",          "mZprime",          20,   0, 5000);
  TH1D* h_mZ               = new TH1D("h_mZ",               "mZ",               30,  60,  120);
  TH1D* h_ptZ              = new TH1D("h_ptZ",              "ptZ",              50,   0, 1000);
  TH1D* h_FATjetPt         = new TH1D("h_FATjetPt",         "FATjetPt",         20, 100, 1000);
  TH1D* h_FATjetSDmass     = new TH1D("h_FATjetSDmass",     "FATjetSDmass",     15,  50,  200);
  TH1D* h_FATjetPRmass     = new TH1D("h_FATjetPRmass",     "FATjetPRmass",     15,  50,  200);
  TH1D* h_FATjetTau2dvTau1 = new TH1D("h_FATjetTau2dvTau1", "FATjetTau2dvTau1", 20,   0,    1);

  h_mZprime         ->Sumw2();
  h_mZ              ->Sumw2();
  h_ptZ             ->Sumw2();
  h_FATjetPt        ->Sumw2();   
  h_FATjetSDmass    ->Sumw2();
  h_FATjetPRmass    ->Sumw2();
  h_FATjetTau2dvTau1->Sumw2();

  h_mZprime         ->GetXaxis()->SetTitle("mZprime");
  h_mZ              ->GetXaxis()->SetTitle("mZ");
  h_ptZ             ->GetXaxis()->SetTitle("ptZ");
  h_FATjetPt        ->GetXaxis()->SetTitle("FATjetPt");
  h_FATjetSDmass    ->GetXaxis()->SetTitle("FATjetSDmass");
  h_FATjetPRmass    ->GetXaxis()->SetTitle("FATjetPRmass");
  h_FATjetTau2dvTau1->GetXaxis()->SetTitle("FATjetTau2dvTau1");
    
  // begin of event loop

  for( Long64_t ev = 0; ev < data.GetEntriesFast(); ev++ ){

    if( ev % 1000000 == 0 )
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());

    data.GetEntry(ev);

    Float_t        eventWeight       = data.GetFloat("ev_weight"); 
    TClonesArray*  muP4              = (TClonesArray*) data.GetPtrTObject("muP4");
    Int_t          FATnJet           = data.GetInt("FATnJet");    
    Int_t*         FATnSubSDJet      = data.GetPtrInt("FATnSubSDJet");
    Float_t*       FATjetSDmass      = data.GetPtrFloat("FATjetSDmass");
    Float_t*       corrPRmass        = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    Float_t*       FATjetTau1        = data.GetPtrFloat("FATjetTau1");
    Float_t*       FATjetTau2        = data.GetPtrFloat("FATjetTau2");
    TClonesArray*  FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    vector<bool>&  FATjetPassIDLoose = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    vector<float>* FATsubjetSDCSV    = data.GetPtrVectorFloat("FATsubjetSDCSV", FATnJet);

    // select good muons
      
    vector<Int_t> goodMuID;

    if( !isPassZmumu(data, goodMuID) ) continue;

    TLorentzVector* thisMu = (TLorentzVector*)muP4->At(goodMuID[0]);
    TLorentzVector* thatMu = (TLorentzVector*)muP4->At(goodMuID[1]);
 
    Float_t mll  = (*thisMu+*thatMu).M();
    Float_t ptll = (*thisMu+*thatMu).Pt();

    h_mZ ->Fill(mll,eventWeight);
    h_ptZ->Fill(ptll,eventWeight);

    // select good FATjet
    
    Int_t goodFATJetID = -1;
    TLorentzVector* thisJet = NULL;

    for(Int_t ij = 0; ij < FATnJet; ij++){

      thisJet = (TLorentzVector*)FATjetP4->At(ij);
      
      if( thisJet->Pt() < 200 ) continue;
      if( fabs(thisJet->Eta()) > 2.4 ) continue;
      if( corrPRmass[ij] < 105 || corrPRmass[ij] > 135 ) continue;
      if( !FATjetPassIDLoose[ij] ) continue;
      if( FATnSubSDJet[ij] != 2 ) continue;
      if( FATsubjetSDCSV[ij][0] < 0.605 || FATsubjetSDCSV[ij][1] < 0.605 ) continue;
      if( thisJet->DeltaR(*thisMu) < 0.8 || thisJet->DeltaR(*thatMu) < 0.8 ) continue;
      
      goodFATJetID = ij;
      break;

    }

    if( goodFATJetID < 0 ) continue; 
    
    h_FATjetPt        ->Fill(thisJet->Pt(),eventWeight);
    h_FATjetSDmass    ->Fill(FATjetSDmass[goodFATJetID],eventWeight);
    h_FATjetPRmass    ->Fill(corrPRmass[goodFATJetID],eventWeight);
    h_FATjetTau2dvTau1->Fill(FATjetTau2[goodFATJetID]/FATjetTau1[goodFATJetID],eventWeight);

    Float_t mllbb = (*thisMu+*thatMu+*thisJet).M();

    if( mllbb < 700 ) continue;

    h_mZprime->Fill(mllbb,eventWeight);

  } // end of event loop

  fprintf(stderr, "Processed all events\n");

  TFile* outFile = new TFile(Form("%s_mZHmu.root",outputFile.c_str()), "recreate");

  h_mZprime         ->Write("mZprime");
  h_mZ              ->Write("mZ");
  h_ptZ             ->Write("ptZ");
  h_FATjetPt        ->Write("FATjetPt");
  h_FATjetSDmass    ->Write("FATjetSDmass");
  h_FATjetPRmass    ->Write("FATjetPRmass");
  h_FATjetTau2dvTau1->Write("FATjetTau2dvTau1");
  h_totalEvents     ->Write("totalEvents");

  outFile->Write();

  delete f;
  delete outFile;
  
}
