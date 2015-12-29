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
#include "../../pileupMCweight.h"

void muVariable(std::string inputFile, std::string outputFile){

  // read the ntuples (in pcncu)

  std::vector<string> infiles;

  readSample(inputFile, infiles);
    
  TreeReader data(infiles);
  
  // Declare the histogram (hightptMuon, customizeTrackerMuon)

  TH1D* h_muHits[2]; 
  TH1D* h_muMatches[2];
  TH1D* h_muTrkLayers[2];  
  TH1D* h_muPixelHits[2];
  TH1D* h_muTrkPtErrdvTrkPt[2];
  TH1D* h_mudxy[2]; 
  TH1D* h_mudz[2];   
  TH1D* h_muMiniIsoEA[2]; 
  TH1D* h_eventWeight[2];

  for(Int_t i = 0; i < 2; i++){

    h_muHits[i]            = new TH1D(Form("h_muHits%d",i),            "muHits",             60,  -0.5, 59.5);
    h_muMatches[i]         = new TH1D(Form("h_muMatches%d",i),         "muMatches",           7,  -0.5,  6.5);
    h_muTrkLayers[i]       = new TH1D(Form("h_muTrkLayers%d",i),       "muTrkLayers",        18,  -0.5, 17.5);
    h_muPixelHits[i]       = new TH1D(Form("h_muPixelHits%d",i),       "muPixelHits",         8,  -0.5,  7.5);
    h_muTrkPtErrdvTrkPt[i] = new TH1D(Form("h_muTrkPtErrdvTrkPt%d",i), "muTrkPtErrdvTrkPt",  20,     0, 0.15);
    h_mudxy[i]             = new TH1D(Form("h_mudxy%d",i),             "mudxy",              20, -0.01, 0.01);
    h_mudz[i]              = new TH1D(Form("h_mudz%d",i),              "mudz",               20, -0.05, 0.05);  
    h_muMiniIsoEA[i]       = new TH1D(Form("h_muMiniIsoEA%d",i),       "muMiniIsoEA",        20,     0, 0.15);
    h_eventWeight[i]       = new TH1D(Form("h_eventWeight%d",i),       "eventWeight",       100,    -1,    1);

    h_muHits[i]           ->Sumw2();
    h_muMatches[i]        ->Sumw2();
    h_muTrkLayers[i]      ->Sumw2();
    h_muPixelHits[i]      ->Sumw2();
    h_muTrkPtErrdvTrkPt[i]->Sumw2();
    h_mudxy[i]            ->Sumw2();
    h_mudz[i]             ->Sumw2();
    h_muMiniIsoEA[i]      ->Sumw2();

    h_muHits[i]           ->GetXaxis()->SetTitle("muHits");
    h_muMatches[i]        ->GetXaxis()->SetTitle("muMatches");
    h_muTrkLayers[i]      ->GetXaxis()->SetTitle("muTrkLayers");
    h_muPixelHits[i]      ->GetXaxis()->SetTitle("muPixelHits");
    h_muTrkPtErrdvTrkPt[i]->GetXaxis()->SetTitle("muTrkPtErrdvTrkPt");
    h_mudxy[i]            ->GetXaxis()->SetTitle("mudxy");
    h_mudz[i]             ->GetXaxis()->SetTitle("mudz");
    h_muMiniIsoEA[i]      ->GetXaxis()->SetTitle("muMiniIsoEA");

  }

  // begin of event loop

  for( Long64_t ev = 0; ev < data.GetEntriesFast(); ev++ ){

    if( ev % 1000000 == 0 )
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());

    data.GetEntry(ev);

    Int_t    nVtx        = data.GetInt("nVtx");
    Int_t    nMu         = data.GetInt("nMu");
    Int_t*   muHits      = data.GetPtrInt("muHits");
    Int_t*   muMatches   = data.GetPtrInt("muMatches");
    Int_t*   muTrkLayers = data.GetPtrInt("muTrkLayers");
    Int_t*   muPixelHits = data.GetPtrInt("muPixelHits");
    Bool_t   isData      = data.GetBool("isData");
    Float_t  pu_nTrueInt = data.GetFloat("pu_nTrueInt");
    Float_t* muTrkPtErr  = data.GetPtrFloat("muTrkPtErr");	
    Float_t* muTrkPt     = data.GetPtrFloat("muTrkPt");
    Float_t* mudxy       = data.GetPtrFloat("mudxy");
    Float_t* mudz        = data.GetPtrFloat("mudz");
    Float_t* muMiniIsoEA = data.GetPtrFloat("muMiniIsoEA");
    TClonesArray* muP4   = (TClonesArray*) data.GetPtrTObject("muP4");
    vector<bool>& isGlobalMuon  = *((vector<bool>*) data.GetPtr("isGlobalMuon"));
    vector<bool>& isTrackerMuon = *((vector<bool>*) data.GetPtr("isTrackerMuon"));

    // remove event which is no hard interaction (noise)

    if( nVtx < 1 ) continue;

    // Correct the pile-up shape of MC

    Double_t eventWeight = pileupWeight(isData, (Int_t)pu_nTrueInt);
    
    h_eventWeight[0]->Fill(0.,eventWeight);
    h_eventWeight[1]->Fill(0.,eventWeight);

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

      if( isGlobalMuon[muId[ie]] ){ // highptMuon selections and cuts

	for(Int_t flag = 0; flag <= 7; flag++){

	  if( muHits[muId[ie]]      <  1   && flag != 0 ) continue;
	  if( muMatches[muId[ie]]   <  2   && flag != 1 ) continue;
	  if( muTrkLayers[muId[ie]] <  6   && flag != 2 ) continue;
	  if( muPixelHits[muId[ie]] <  1   && flag != 3 ) continue;
	  if( muTrkPtErr[muId[ie]]/muTrkPt[muId[ie]] > 0.3 && flag != 4 ) continue;
	  if( fabs(mudxy[muId[ie]]) >  0.2 && flag != 5 ) continue;
	  if( fabs(mudz[muId[ie]])  >  0.5 && flag != 6 ) continue;
	  if( muMiniIsoEA[muId[ie]] >  0.2 && flag != 7 ) continue;
	  	  	    
	  switch(flag){

	  case 0: h_muHits[0]      ->Fill(muHits[muId[ie]],eventWeight);      break;
	  case 1: h_muMatches[0]   ->Fill(muMatches[muId[ie]],eventWeight);   break;
	  case 2: h_muTrkLayers[0] ->Fill(muTrkLayers[muId[ie]],eventWeight); break;
	  case 3: h_muPixelHits[0] ->Fill(muPixelHits[muId[ie]],eventWeight); break;
	  case 4: h_muTrkPtErrdvTrkPt[0] ->Fill(muTrkPtErr[muId[ie]]/muTrkPt[muId[ie]],eventWeight); break;
	  case 5: h_mudxy[0]       ->Fill(mudxy[muId[ie]],eventWeight);       break;	 
	  case 6: h_mudz[0]        ->Fill(mudz[muId[ie]],eventWeight);        break;
	  case 7: h_muMiniIsoEA[0] ->Fill(muMiniIsoEA[muId[ie]],eventWeight); break;
	    
	  } // end of switch
    
	} // end of flag loop

      } // end of barrel


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

	  case 0: h_muMatches[1]   ->Fill(muMatches[muId[ie]],eventWeight);   break;
          case 1: h_muTrkLayers[1] ->Fill(muTrkLayers[muId[ie]],eventWeight); break;
          case 2: h_muPixelHits[1] ->Fill(muPixelHits[muId[ie]],eventWeight); break;
          case 3: h_muTrkPtErrdvTrkPt[1] ->Fill(muTrkPtErr[muId[ie]]/muTrkPt[muId[ie]],eventWeight); break;
          case 4: h_mudxy[1]       ->Fill(mudxy[muId[ie]],eventWeight);       break;
          case 5: h_mudz[1]        ->Fill(mudz[muId[ie]],eventWeight);        break;
          case 6: h_muMiniIsoEA[1] ->Fill(muMiniIsoEA[muId[ie]],eventWeight); break;
	  case 7: h_muHits[1]      ->Fill(muHits[muId[ie]],eventWeight);      break;
    
	  } // end of switch
	
	} // end of flag loop

      } // end of endcap
	
    } // end of two muctrons loop
      
  } // end of event loop

  fprintf(stderr, "Processed all events\n");
    
  TFile* outFile[2];

  std::string region[2] = {"highptMuon","customizeTrackerMuon"};

  for(Int_t i = 0; i < 2; i++){

    outFile[i] = new TFile(Form("%s_%s.root",outputFile.c_str(),region[i].c_str()), "recreate");
      
    h_muMatches[i]        ->Write("muMatches");
    h_muTrkLayers[i]      ->Write("muTrkLayers");
    h_muPixelHits[i]      ->Write("muPixelHits");
    h_muTrkPtErrdvTrkPt[i]->Write("muTrkPtErrdvTrkPt");
    h_mudxy[i]            ->Write("mudxy");
    h_mudz[i]             ->Write("mudz");
    h_muMiniIsoEA[i]      ->Write("muMiniIsoEA");
    h_muHits[i]           ->Write("muHits");
    h_eventWeight[i]      ->Write("eventWeight");

    outFile[i]->Write();

  }
  
}
