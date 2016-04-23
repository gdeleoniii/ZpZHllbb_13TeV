#include <vector>
#include <string>
#include <iostream>
#include <TH1D.h>
#include <TFile.h>
#include <TProfile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "../../untuplizer.h"

void eleEndcapVariable(std::string inputFile, std::string outputFile){

  // read the ntuples (in pcncu)

  TreeReader data(inputFile.data());
  
  // Declare the histogram

  TFile* f = new TFile(inputFile.data());
  TH1D* h_totalEvents             = (TH1D*)f->Get("h_totalEv");

  TH1D* h_eleEtaseedAtVtx         = new TH1D("h_eleEtaseedAtVtx",         "eleEtaseedAtVtx",         20,  -0.01,  0.01);
  TH1D* h_eledPhiAtVtx            = new TH1D("h_eledPhiAtVtx",            "eledPhiAtVtx",            20,  -0.03,  0.03);
  TH1D* h_eleHoverE               = new TH1D("h_eleHoverE",               "eleHoverE",               20,      0,  0.05);
  TH1D* h_eleSigmaIEtaIEtaFull5x5 = new TH1D("h_eleSigmaIEtaIEtaFull5x5", "eleSigmaIEtaIEtaFull5x5", 40,      0,  0.06);
  TH1D* h_eleFull5x5E2x5dvE5x5    = new TH1D("h_eleFull5x5E2x5dvE5x5",    "eleFull5x5E2x5dvE5x5",    20,      0,     1);
  TH1D* h_eleFull5x5E1x5dvE5x5    = new TH1D("h_eleFull5x5E1x5dvE5x5",    "eleFull5x5E1x5dvE5x5",    20,      0,     1);
  TH1D* h_eleMissHits             = new TH1D("h_eleMissHits",             "eleMissHits",              6,   -0.5,   5.5);
  TH1D* h_eleD0                   = new TH1D("h_eleD0",                   "eleD0",                   20, -0.015, 0.015);  
  TH1D* h_eleMiniIsoEA            = new TH1D("h_eleMiniIsoEA",            "eleMiniIsoEA",            20,      0,  0.12);

  h_eleEtaseedAtVtx        ->Sumw2();
  h_eledPhiAtVtx           ->Sumw2();
  h_eleHoverE              ->Sumw2();
  h_eleSigmaIEtaIEtaFull5x5->Sumw2();
  h_eleFull5x5E2x5dvE5x5   ->Sumw2();
  h_eleFull5x5E1x5dvE5x5   ->Sumw2();
  h_eleMissHits            ->Sumw2(); 
  h_eleD0                  ->Sumw2();
  h_eleMiniIsoEA           ->Sumw2();

  h_eleEtaseedAtVtx        ->GetXaxis()->SetTitle("eleEtaseedAtVtx"); 
  h_eledPhiAtVtx           ->GetXaxis()->SetTitle("eledPhiAtVtx");   
  h_eleHoverE              ->GetXaxis()->SetTitle("eleHoverE");    
  h_eleSigmaIEtaIEtaFull5x5->GetXaxis()->SetTitle("eleSigmaIEtaIEtaFull5x5");  
  h_eleFull5x5E2x5dvE5x5   ->GetXaxis()->SetTitle("eleFull5x5E2x5dvE5x5");
  h_eleFull5x5E1x5dvE5x5   ->GetXaxis()->SetTitle("eleFull5x5E1x5dvE5x5");   
  h_eleMissHits            ->GetXaxis()->SetTitle("eleMissHits"); 
  h_eleD0                  ->GetXaxis()->SetTitle("eleD0");      
  h_eleMiniIsoEA           ->GetXaxis()->SetTitle("eleMiniIsoEA"); 

  // begin of event loop

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    if( (unsigned)ev % 500000 == 0 )
      fprintf(stdout, "Still left events %lli of %lli\n", ev, data.GetEntriesFast());

    data.GetEntry(ev);

    Int_t    nEle                    = data.GetInt("nEle");
    Int_t*   eleMissHits             = data.GetPtrInt("eleMissHits");
    Float_t  eventWeight             = data.GetFloat("ev_weight");
    Float_t* eleScEn                 = data.GetPtrFloat("eleScEn");
    Float_t* eleScEt                 = data.GetPtrFloat("eleScEt");
    Float_t* eleScEta                = data.GetPtrFloat("eleScEta");
    Float_t* eleEtaseedAtVtx         = data.GetPtrFloat("eleEtaseedAtVtx");	
    Float_t* eledPhiAtVtx            = data.GetPtrFloat("eledPhiAtVtx");
    Float_t* eleHoverE               = data.GetPtrFloat("eleHoverE");
    Float_t* eleSigmaIEtaIEtaFull5x5 = data.GetPtrFloat("eleSigmaIEtaIEtaFull5x5");
    Float_t* eleE1x5Full5x5          = data.GetPtrFloat("eleE1x5Full5x5");
    Float_t* eleE2x5Full5x5          = data.GetPtrFloat("eleE2x5Full5x5");
    Float_t* eleE5x5Full5x5          = data.GetPtrFloat("eleE5x5Full5x5");   
    Float_t* eleD0                   = data.GetPtrFloat("eleD0");
    Float_t* eleMiniIsoEA            = data.GetPtrFloat("eleMiniIsoEA");
    TClonesArray* eleP4              = (TClonesArray*) data.GetPtrTObject("eleP4");
    vector<bool>& eleEcalDrivenSeed  = *((vector<bool>*) data.GetPtr("eleEcalDrivenSeed"));

    // choosing electron

    vector<Int_t> eleId;

    bool findEPair = false;
    
    for(Int_t ie = 0; ie < nEle; ie++){

      if( !(fabs(eleScEta[ie]) < 1.442 || fabs(eleScEta[ie]) > 1.566) ) continue;
      if( fabs(eleScEta[ie]) > 2.5 ) continue;

      TLorentzVector* thisEle = (TLorentzVector*)eleP4->At(ie);
      
      for(Int_t je = 0; je < ie; je++){

	if( eleScEt[ie] <= 35 || eleScEt[je] <= 35 ) continue;
	if( !eleEcalDrivenSeed[ie] || !eleEcalDrivenSeed[je] ) continue;
	if( !(fabs(eleScEta[je]) < 1.442 || fabs(eleScEta[je]) > 1.566) ) continue;
	if( fabs(eleScEta[je]) > 2.5 ) continue;

	TLorentzVector* thatEle = (TLorentzVector*)eleP4->At(je);
	
	Float_t mll  = (*thisEle+*thatEle).M();
	Float_t ptll = (*thisEle+*thatEle).Pt();

	if( mll < 60 || mll > 120 ) continue;		
	if( ptll < 150 ) continue;

	if( !findEPair ){

	  eleId.push_back(ie);
	  eleId.push_back(je);

	} 

	findEPair = true;
	break;
	
      }
    }

    if(	!findEPair ) continue;
    
    // electron selections and cuts

    for(Int_t ie = 0; ie < 2; ie++){

      Float_t E = eleScEn[eleId[ie]];

      if( fabs(eleScEta[eleId[ie]]) > 1.566 && fabs(eleScEta[eleId[ie]]) < 2.5 ){ // endcap selections and cuts

	for(Int_t flag = 0; flag <= 7; flag++){

	  if( fabs(eleEtaseedAtVtx[eleId[ie]])   >= 0.006    && flag != 0 ) continue;
	  if( fabs(eledPhiAtVtx[eleId[ie]])      >= 0.06     && flag != 1 ) continue;
	  if( eleHoverE[eleId[ie]]               >= 5/E+0.05 && flag != 2 ) continue;
	  if( eleSigmaIEtaIEtaFull5x5[eleId[ie]] >= 0.03     && flag != 3 ) continue;
	  if( eleMissHits[eleId[ie]]             >  1        && flag != 4 ) continue;
	  if( fabs(eleD0[eleId[ie]])             >= 0.05     && flag != 5 ) continue;
	  if( eleMiniIsoEA[eleId[ie]]            >= 0.1      && flag != 6 ) continue;
	      
	  switch(flag){

	  case 0: h_eleEtaseedAtVtx ->Fill(eleEtaseedAtVtx[eleId[ie]],eventWeight); break;
	  case 1: h_eledPhiAtVtx    ->Fill(eledPhiAtVtx[eleId[ie]],eventWeight);    break;
	  case 2: h_eleHoverE       ->Fill(eleHoverE[eleId[ie]],eventWeight);       break;
	  case 3: h_eleSigmaIEtaIEtaFull5x5 ->Fill(eleSigmaIEtaIEtaFull5x5[eleId[ie]],eventWeight); break;
	  case 4: h_eleMissHits     ->Fill(eleMissHits[eleId[ie]],eventWeight);     break;
	  case 5: h_eleD0           ->Fill(eleD0[eleId[ie]],eventWeight);           break;
	  case 6: h_eleMiniIsoEA    ->Fill(eleMiniIsoEA[eleId[ie]],eventWeight);    break;  
	  case 7: 
	    h_eleFull5x5E2x5dvE5x5  ->Fill(eleE2x5Full5x5[eleId[ie]]/eleE5x5Full5x5[eleId[ie]],eventWeight);
	    h_eleFull5x5E1x5dvE5x5  ->Fill(eleE1x5Full5x5[eleId[ie]]/eleE5x5Full5x5[eleId[ie]],eventWeight);
	    break;
	        
	  } // end of switch
	  
	} // end of flag loop

      } // end of endcap
	
    } // end of two electrons loop
      
  } // end of event loop

  fprintf(stderr, "Processed all events\n");
    
  TFile* outFile = new TFile(Form("%s_eleEndcapVariable.root",outputFile.c_str()), "recreate");
      
  h_eleEtaseedAtVtx        ->Write("eleEtaseedAtVtx");
  h_eledPhiAtVtx           ->Write("eledPhiAtVtx");
  h_eleHoverE              ->Write("eleHoverE");
  h_eleSigmaIEtaIEtaFull5x5->Write("eleSigmaIEtaIEtaFull5x5");
  h_eleFull5x5E2x5dvE5x5   ->Write("eleFull5x5E2x5dvE5x5");
  h_eleFull5x5E1x5dvE5x5   ->Write("eleFull5x5E1x5dvE5x5");
  h_eleMissHits            ->Write("eleMissHits");
  h_eleD0                  ->Write("eleD0");
  h_eleMiniIsoEA           ->Write("eleMiniIsoEA");
  h_totalEvents            ->Write("totalEvents");

  outFile->Write();

  delete outFile;  
  delete f;

}
