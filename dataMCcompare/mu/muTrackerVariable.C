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

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    if( (unsigned)ev % 500000 == 0 )
      fprintf(stdout, "Still left events %lli of %lli\n", ev, data.GetEntriesFast());

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

    Int_t muId[2] = {-1,-1};

    for(Int_t ie = 0; ie < nMu; ie++){

      TLorentzVector* thisMu = (TLorentzVector*)muP4->At(ie);

      if( !isGlobalMuon[ie] || !isTrackerMuon[ie] ) continue;
      if( thisMu->Pt() <= 50 ) continue;
      if( fabs(thisMu->Eta()) >= 2.1 ) continue;    

      for(Int_t je = 0; je < ie; je++){

	TLorentzVector* thatMu = (TLorentzVector*)muP4->At(je);
	
	if( !isGlobalMuon[je] || !isTrackerMuon[je] ) continue;
	if( thatMu->Pt() <= 50 ) continue;
	if( fabs(thatMu->Eta()) >= 2.1 ) continue;
		
	Float_t mll  = (*thisMu+*thatMu).M();
	Float_t ptll = (*thisMu+*thatMu).Pt();

	if( mll < 60 || mll > 120 ) continue;		
	if( ptll < 150 ) continue;

	muId[0] = ie;
	muId[1] = je;

	break;
	
      }
    }

    if(	muId[0] < 0 || muId[1] < 0 ) continue;

    // muon selections and cuts

    for(Int_t ie = 0; ie < 2; ie++){

      if( isTrackerMuon[muId[ie]] ){ // customizeTrackerMuon selections and cuts

	for(Int_t flag = 0; flag <= 7; flag++){

	  if( muMatches[muId[ie]]   <  2   && flag != 0 ) continue;
	  if( muTrkLayers[muId[ie]] <  6   && flag != 1 ) continue;
	  if( muPixelHits[muId[ie]] <  1   && flag != 2 ) continue;
	  if( muTrkPtErr[muId[ie]]/muTrkPt[muId[ie]] > 0.3 && flag != 3 ) continue;
	  if( fabs(mudxy[muId[ie]]) >  0.2 && flag != 4 ) continue;
	  if( fabs(mudz[muId[ie]])  >  0.5 && flag != 5 ) continue;
	  if( muMiniIsoEA[muId[ie]] >  0.2 && flag != 6 ) continue;

	  switch(flag){

	  case 0: h_muMatches   ->Fill(muMatches[muId[ie]],eventWeight);   break;
	  case 1: h_muTrkLayers ->Fill(muTrkLayers[muId[ie]],eventWeight); break;
	  case 2: h_muPixelHits ->Fill(muPixelHits[muId[ie]],eventWeight); break;
	  case 3: h_muTrkPtErrdvTrkPt ->Fill(muTrkPtErr[muId[ie]]/muTrkPt[muId[ie]],eventWeight); break;
	  case 4: h_mudxy       ->Fill(mudxy[muId[ie]],eventWeight);       break;
	  case 5: h_mudz        ->Fill(mudz[muId[ie]],eventWeight);        break;
	  case 6: h_muMiniIsoEA ->Fill(muMiniIsoEA[muId[ie]],eventWeight); break;
	  case 7: h_muHits      ->Fill(muHits[muId[ie]],eventWeight);      break;
    
	  } // end of switch
	
	} // end of flag loop

      } // end of customizeTrackerMuon
	
    } // end of two muons loop
      
  } // end of event loop

  fprintf(stderr, "Processed all events\n");
    
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
