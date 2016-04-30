#include <vector>
#include <string>
#include <iostream>
#include <TH1D.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "../../untuplizer.h"

void muTrackerVariable(std::string inputFile, std::string outputFile){

  // read the ntuples (in pcncu)

  TreeReader data(inputFile.data());

  // Declare the histogram 

  TFile* f = new TFile(inputFile.data());

  TH1D* h_totalEvents       = (TH1D*)f->Get("h_totalEv");

  TH1D* h_muHits            = new TH1D("h_muHits",            "muHits",             60,  -0.5, 59.5);
  TH1D* h_muMatches         = new TH1D("h_muMatches",         "muMatches",           7,  -0.5,  6.5);
  TH1D* h_muTrkLayers       = new TH1D("h_muTrkLayers",       "muTrkLayers",        18,  -0.5, 17.5);
  TH1D* h_muPixelHits       = new TH1D("h_muPixelHits",       "muPixelHits",         8,  -0.5,  7.5);
  TH1D* h_muTrkPtErrdvTrkPt = new TH1D("h_muTrkPtErrdvTrkPt", "muTrkPtErrdvTrkPt",  20,     0, 0.15);
  TH1D* h_mudxy             = new TH1D("h_mudxy",             "mudxy",              20, -0.01, 0.01);
  TH1D* h_mudz              = new TH1D("h_mudz",              "mudz",               20, -0.05, 0.05);  
  TH1D* h_muMiniIsoEA       = new TH1D("h_muMiniIsoEA",       "muMiniIsoEA",        20,     0, 0.15);

  h_muHits           ->Sumw2();
  h_muMatches        ->Sumw2();
  h_muTrkLayers      ->Sumw2();
  h_muPixelHits      ->Sumw2();
  h_muTrkPtErrdvTrkPt->Sumw2();
  h_mudxy            ->Sumw2();
  h_mudz             ->Sumw2();
  h_muMiniIsoEA      ->Sumw2();

  h_muHits           ->GetXaxis()->SetTitle("muHits");
  h_muMatches        ->GetXaxis()->SetTitle("muMatches");
  h_muTrkLayers      ->GetXaxis()->SetTitle("muTrkLayers");
  h_muPixelHits      ->GetXaxis()->SetTitle("muPixelHits");
  h_muTrkPtErrdvTrkPt->GetXaxis()->SetTitle("muTrkPtErrdvTrkPt");
  h_mudxy            ->GetXaxis()->SetTitle("mudxy");
  h_mudz             ->GetXaxis()->SetTitle("mudz");
  h_muMiniIsoEA      ->GetXaxis()->SetTitle("muMiniIsoEA");

  // begin of event loop

  fprintf(stdout, "Total events %lli\n", data.GetEntriesFast());

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    if( (unsigned)ev % 10000 == 0 )
      fprintf(stdout, "Still left events %lli\n", ev);

    data.GetEntry(ev);

    Float_t  eventWeight = data.GetFloat("ev_weight");
    Int_t    nMu         = data.GetInt("nMu");
    Int_t*   muHits      = data.GetPtrInt("muHits");
    Int_t*   muMatches   = data.GetPtrInt("muMatches");
    Int_t*   muTrkLayers = data.GetPtrInt("muTrkLayers");
    Int_t*   muPixelHits = data.GetPtrInt("muPixelHits");
    Float_t* muTrkPtErr  = data.GetPtrFloat("muTrkPtErr");	
    Float_t* muTrkPt     = data.GetPtrFloat("muTrkPt");
    Float_t* mudxy       = data.GetPtrFloat("mudxy");
    Float_t* mudz        = data.GetPtrFloat("mudz");
    Float_t* muMiniIsoEA = data.GetPtrFloat("muMiniIsoEA");
    TClonesArray* muP4   = (TClonesArray*) data.GetPtrTObject("muP4");
    vector<bool>& isGlobalMuon  = *((vector<bool>*) data.GetPtr("isGlobalMuon"));
    vector<bool>& isTrackerMuon = *((vector<bool>*) data.GetPtr("isTrackerMuon"));

    // choosing muon pair

    vector<Int_t> muId;

    bool findMPair = false;

    for(Int_t im = 0; im < nMu; ++im){
      for(Int_t jm = 0; jm < im; ++jm){

	TLorentzVector* thisMu = (TLorentzVector*)muP4->At(im);
	TLorentzVector* thatMu = (TLorentzVector*)muP4->At(jm);
	
	if( !isGlobalMuon[im] || !isTrackerMuon[im] ) continue;
	if( !isGlobalMuon[jm] || !isTrackerMuon[jm] ) continue;
	if( thisMu->Pt() < 50 || thatMu->Pt() < 50  ) continue;
	if( fabs(thisMu->Eta()) > 2.1  || fabs(thatMu->Eta()) > 2.1   ) continue;
	if( (*thisMu+*thatMu).M() < 60 || (*thisMu+*thatMu).M() > 120 ) continue;
	if( (*thisMu+*thatMu).Pt() < 150 ) continue;

	if( !findMPair ){
	  muId.push_back(im);
	  muId.push_back(jm);
	} 

	break;
	
      }
    }

    if( !findMPair ) continue;

    // muon selections and cuts

    for(Int_t im = 0; im < 2; ++im){

      if( isTrackerMuon[muId[im]] ){ // customizeTrackerMuon selections and cuts

	for(Int_t flag = 0; flag <= 7; ++flag){

	  if( muMatches[muId[im]]   <  2   && flag != 0 ) continue;
	  if( muTrkLayers[muId[im]] <  6   && flag != 1 ) continue;
	  if( muPixelHits[muId[im]] <  1   && flag != 2 ) continue;
	  if( muTrkPtErr[muId[im]]/muTrkPt[muId[im]] > 0.3 && flag != 3 ) continue;
	  if( fabs(mudxy[muId[im]]) >  0.2 && flag != 4 ) continue;
	  if( fabs(mudz[muId[im]])  >  0.5 && flag != 5 ) continue;
	  if( muMiniIsoEA[muId[im]] >  0.2 && flag != 6 ) continue;

	  switch(flag){

	  case 0: h_muMatches   ->Fill(muMatches[muId[im]],eventWeight);   break;
	  case 1: h_muTrkLayers ->Fill(muTrkLayers[muId[im]],eventWeight); break;
	  case 2: h_muPixelHits ->Fill(muPixelHits[muId[im]],eventWeight); break;
	  case 3: h_muTrkPtErrdvTrkPt ->Fill(muTrkPtErr[muId[im]]/muTrkPt[muId[im]],eventWeight); break;
	  case 4: h_mudxy       ->Fill(mudxy[muId[im]],eventWeight);       break;
	  case 5: h_mudz        ->Fill(mudz[muId[im]],eventWeight);        break;
	  case 6: h_muMiniIsoEA ->Fill(muMiniIsoEA[muId[im]],eventWeight); break;
	  case 7: h_muHits      ->Fill(muHits[muId[im]],eventWeight);      break;
    
	  } // end of switch
	
	} // end of flag loop

      } // end of customizeTrackerMuon
	
    } // end of two muons loop
      
  } // end of event loop

  fprintf(stdout, "Processed all events\n");
    
  TFile* outFile = new TFile(Form("%s_muTrackerVariable.root",outputFile.c_str()), "recreate");
      
  h_muMatches        ->Write("muMatches");
  h_muTrkLayers      ->Write("muTrkLayers");
  h_muPixelHits      ->Write("muPixelHits");
  h_muTrkPtErrdvTrkPt->Write("muTrkPtErrdvTrkPt");
  h_mudxy            ->Write("mudxy");
  h_mudz             ->Write("mudz");
  h_muMiniIsoEA      ->Write("muMiniIsoEA");
  h_muHits           ->Write("muHits");
  h_totalEvents      ->Write("totalEvents");

  outFile->Write();

  delete outFile;
  delete f;
  
}
