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

void muZjetVariable(string inputFile, string outputFile){

  // read the ntuples (in pcncu)

  TreeReader data(inputFile.data());
  
  // Declare the histogram
     
  TFile f(inputFile.data());
  TH1D* h_totalEvents = (TH1D*)f.Get("h_totalEv");

  TH1D* h_Zmass            = new TH1D("h_Zmass",            "Zmass",            30,   60,  120);
  TH1D* h_Zpt              = new TH1D("h_Zpt",              "Zpt",              50,    0, 1000);
  TH1D* h_Zeta             = new TH1D("h_Zeta",             "Zeta",             40,   -4,    4);
  TH1D* h_ZRapidity        = new TH1D("h_ZRapidity",        "ZRapidity",        40,   -4,    4);
  TH1D* h_leadLepPt        = new TH1D("h_leadLepPt",        "leadLepPt",        16,    0,  800);
  TH1D* h_leadLepEta       = new TH1D("h_leadLepEta",       "leadLepEta",       40,   -4,    4);
  TH1D* h_subleadLepPt     = new TH1D("h_subleadLepPt",     "subleadLepPt",     25,    0,  500);
  TH1D* h_subleadLepEta    = new TH1D("h_subleadLepEta",    "subleadLepEta",    40,   -4,    4);
  TH1D* h_nVtx             = new TH1D("h_nVtx",             "nVtx",             30, -0.5, 29.5);     
  TH1D* h_FATjetPt         = new TH1D("h_FATjetPt",         "FATjetPt",         20,  100, 1000);
  TH1D* h_FATjetEta        = new TH1D("h_FATjetEta",        "FATjetEta",        20,   -4,    4);
  TH1D* h_FATjetCISVV2     = new TH1D("h_FATjetCISVV2",     "FATjetCISVV2",     20,    0,  1.2);
  TH1D* h_FATjetSDmass     = new TH1D("h_FATjetSDmass",     "FATjetSDmass",     20,    0,  200);
  TH1D* h_FATjetPRmass     = new TH1D("h_FATjetPRmass",     "FATjetPRmass",     20,    0,  200);
  TH1D* h_FATjetPRmassCorr = new TH1D("h_FATjetPRmassCorr", "FATjetPRmassCorr", 20,    0,  200);
  TH1D* h_FATsubjetPt      = new TH1D("h_FATsubjetPt",      "FATsubjetPt",      20,    0,  800);
  TH1D* h_FATsubjetEta     = new TH1D("h_FATsubjetEta",     "FATsubjetEta",     20,   -4,    4);  
  TH1D* h_FATsubjetSDCSV   = new TH1D("h_FATsubjetSDCSV",   "FATsubjetSDCSV",   20,    0,  1.2);

  h_Zmass           ->Sumw2();
  h_Zpt             ->Sumw2();
  h_Zeta            ->Sumw2();
  h_ZRapidity       ->Sumw2();
  h_leadLepPt       ->Sumw2();
  h_leadLepEta      ->Sumw2();
  h_subleadLepPt    ->Sumw2();
  h_subleadLepEta   ->Sumw2(); 
  h_nVtx            ->Sumw2();
  h_FATjetPt        ->Sumw2();   
  h_FATjetEta       ->Sumw2();
  h_FATjetCISVV2    ->Sumw2();
  h_FATjetSDmass    ->Sumw2();
  h_FATjetPRmass    ->Sumw2();
  h_FATjetPRmassCorr->Sumw2();
  h_FATsubjetPt     ->Sumw2();
  h_FATsubjetEta    ->Sumw2();
  h_FATsubjetSDCSV  ->Sumw2();
  
  h_Zmass           ->GetXaxis()->SetTitle("Zmass (GeV)"); 
  h_Zpt             ->GetXaxis()->SetTitle("Zpt (GeV)");   
  h_Zeta            ->GetXaxis()->SetTitle("Zeta");    
  h_ZRapidity       ->GetXaxis()->SetTitle("ZRapidity");
  h_leadLepPt       ->GetXaxis()->SetTitle("leadLepPt (GeV)");  
  h_leadLepEta      ->GetXaxis()->SetTitle("leadLepEta");
  h_subleadLepPt    ->GetXaxis()->SetTitle("subleadLepPt (GeV)");   
  h_subleadLepEta   ->GetXaxis()->SetTitle("subleadLepEta"); 
  h_nVtx            ->GetXaxis()->SetTitle("nVtx");
  h_FATjetPt        ->GetXaxis()->SetTitle("FATjetPt (GeV)");
  h_FATjetEta       ->GetXaxis()->SetTitle("FATjetEta");
  h_FATjetCISVV2    ->GetXaxis()->SetTitle("FATjetCISVV2");
  h_FATjetSDmass    ->GetXaxis()->SetTitle("FATjetSDmass (GeV)");
  h_FATjetPRmass    ->GetXaxis()->SetTitle("FATjetPRmass (GeV)");
  h_FATjetPRmassCorr->GetXaxis()->SetTitle("FATjetPRmassL2L3Corr (GeV)");
  h_FATsubjetPt     ->GetXaxis()->SetTitle("FATsubjetPt (GeV)");
  h_FATsubjetEta    ->GetXaxis()->SetTitle("FATsubjetEta");
  h_FATsubjetSDCSV  ->GetXaxis()->SetTitle("FATsubjetSDCSV");
    
  // begin of event loop

  fprintf(stdout, "Total events %lli\n", data.GetEntriesFast());

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    if( (unsigned)ev % 100000 == 0 )
      fprintf(stdout, "Still left events %lli\n", ev);

    data.GetEntry(ev);

    Int_t          nVtx              = data.GetInt("nVtx");
    Float_t        eventWeight       = data.GetFloat("ev_weight");
    TClonesArray*  muP4              = (TClonesArray*) data.GetPtrTObject("muP4");
    Int_t          FATnJet           = data.GetInt("FATnJet");    
    Int_t*         FATnSubSDJet      = data.GetPtrInt("FATnSubSDJet");
    Float_t*       FATjetCISVV2      = data.GetPtrFloat("FATjetCISVV2");
    Float_t*       FATjetSDmass      = data.GetPtrFloat("FATjetSDmass");
    Float_t*       FATjetPRmass      = data.GetPtrFloat("FATjetPRmass");
    Float_t*       FATjetPRmassCorr  = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    TClonesArray*  FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    vector<float>* FATsubjetSDPx     = data.GetPtrVectorFloat("FATsubjetSDPx", FATnJet);
    vector<float>* FATsubjetSDPy     = data.GetPtrVectorFloat("FATsubjetSDPy", FATnJet);
    vector<float>* FATsubjetSDPz     = data.GetPtrVectorFloat("FATsubjetSDPz", FATnJet);
    vector<float>* FATsubjetSDE      = data.GetPtrVectorFloat("FATsubjetSDE", FATnJet);
    vector<float>* FATsubjetSDCSV    = data.GetPtrVectorFloat("FATsubjetSDCSV", FATnJet);

    // select good electrons

    vector<Int_t> goodLepID;

    if( !isPassZmumu(data,goodLepID) ) continue;

    TLorentzVector* thisLep = (TLorentzVector*)muP4->At(goodLepID[0]);
    TLorentzVector* thatLep = (TLorentzVector*)muP4->At(goodLepID[1]);

    h_Zmass    ->Fill((*thisLep+*thatLep).M(),eventWeight);
    h_Zpt      ->Fill((*thisLep+*thatLep).Pt(),eventWeight);
    h_Zeta     ->Fill((*thisLep+*thatLep).Eta(),eventWeight);
    h_ZRapidity->Fill((*thisLep+*thatLep).Rapidity(),eventWeight);

    if( thisLep->Pt() > thatLep->Pt() ){

      h_leadLepPt    ->Fill(thisLep->Pt(),eventWeight);
      h_leadLepEta   ->Fill(thisLep->Eta(),eventWeight);
      h_subleadLepPt ->Fill(thatLep->Pt(),eventWeight);
      h_subleadLepEta->Fill(thatLep->Eta(),eventWeight);

    }else{

      h_leadLepPt    ->Fill(thatLep->Pt(),eventWeight);
      h_leadLepEta   ->Fill(thatLep->Eta(),eventWeight);
      h_subleadLepPt ->Fill(thisLep->Pt(),eventWeight);
      h_subleadLepEta->Fill(thisLep->Eta(),eventWeight);

    }

    // select good FATjet (inclusive, only side bands)

    Int_t goodFATJetID = -1;

    if( !isPassJet(data, &goodFATJetID, thisLep, thatLep, false, true) ) continue;

    TLorentzVector* thisJet = (TLorentzVector*)FATjetP4->At(goodFATJetID);

    Float_t mllbb;

    if( !noiseCleaning(thisLep, thatLep, thisJet, &mllbb) ) continue;

    h_nVtx            ->Fill(nVtx,eventWeight);
    h_FATjetPt        ->Fill(thisJet->Pt(),eventWeight);
    h_FATjetEta       ->Fill(thisJet->Eta(),eventWeight);
    h_FATjetCISVV2    ->Fill(FATjetCISVV2[goodFATJetID],eventWeight);
    h_FATjetSDmass    ->Fill(FATjetSDmass[goodFATJetID],eventWeight);
    h_FATjetPRmass    ->Fill(FATjetPRmass[goodFATJetID],eventWeight);
    h_FATjetPRmassCorr->Fill(FATjetPRmassCorr[goodFATJetID],eventWeight);

    for(Int_t is = 0; is < FATnSubSDJet[goodFATJetID]; ++is){

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

  fprintf(stdout, "Processed all events\n");

  TFile* outFile = new TFile(Form("%s_muZjetVariable.root",outputFile.data()), "recreate");
      
  h_Zmass           ->Write("Zmass");
  h_Zpt             ->Write("Zpt");
  h_Zeta            ->Write("Zeta");
  h_ZRapidity       ->Write("ZRapidity");
  h_leadLepPt       ->Write("leadLepPt");
  h_leadLepEta      ->Write("leadLepEta");
  h_subleadLepPt    ->Write("subleadLepPt");
  h_subleadLepEta   ->Write("subleadLepEta");
  h_nVtx            ->Write("nVtx");
  h_FATjetPt        ->Write("FATjetPt");
  h_FATjetEta       ->Write("FATjetEta");
  h_FATjetCISVV2    ->Write("FATjetCISVV2");
  h_FATjetSDmass    ->Write("FATjetSDmass");
  h_FATjetPRmass    ->Write("FATjetPRmass");
  h_FATjetPRmassCorr->Write("FATjetPRmassCorr");
  h_FATsubjetPt     ->Write("FATsubjetPt");
  h_FATsubjetEta    ->Write("FATsubjetEta");
  h_FATsubjetSDCSV  ->Write("FATsubjetSDCSV");
  h_totalEvents     ->Write("totalEvents");
  
  outFile->Write();

}
