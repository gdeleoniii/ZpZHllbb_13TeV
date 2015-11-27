#include <vector>
#include <string>
#include <iostream>
#include <TH1D.h>
#include <TMath.h>
#include <TFile.h>
#include <TProfile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "../../untuplizer.h"
#include "../../readSample.h"

void eleVariable(std::string inputFile, std::string outputFile){

  // read the ntuples (in pcncu)

  std::vector<string> infiles;

  readSample(inputFile, infiles);
  
  TreeReader data(infiles);
  
  // Declare the histogram (barrel,endcap)

  TH1D* h_eleEtaseedAtVtx[2]; 
  TH1D* h_eledPhiAtVtx[2];
  TH1D* h_eleHoverE[2];  
  TH1D* h_eleSigmaIEtaIEtaFull5x5[2];
  TH1D* h_eleFull5x5E2x5dvE5x5[2];
  TH1D* h_eleFull5x5E1x5dvE5x5[2]; 
  TH1D* h_eleMissHits[2]; 
  TH1D* h_eleD0[2];   
  TH1D* h_eleMiniIsoEA[2]; 
  TH1D* h_eventWeight[2];

  for(Int_t i = 0; i < 2; i++){

    Int_t nBin = 20;    

    h_eleEtaseedAtVtx[i]         = new TH1D(Form("h_eleEtaseedAtVtx%d",i), "eleEtaseedAtVtx", nBin, -0.01, 0.01);
    h_eledPhiAtVtx[i]            = new TH1D(Form("h_eledPhiAtVtx%d",i), "eledPhiAtVtx", nBin, -0.03, 0.03);
    h_eleHoverE[i]               = new TH1D(Form("h_eleHoverE%d",i), "eleHoverE", nBin, 0, 0.05);
    h_eleSigmaIEtaIEtaFull5x5[i] = new TH1D(Form("h_eleSigmaIEtaIEtaFull5x5%d",i), "eleSigmaIEtaIEtaFull5x5", nBin, 0, 0.015);
    h_eleFull5x5E2x5dvE5x5[i]    = new TH1D(Form("h_eleFull5x5E2x5dvE5x5%d",i), "eleFull5x5E2x5dvE5x5", nBin, 0.8, 1.2);
    h_eleFull5x5E1x5dvE5x5[i]    = new TH1D(Form("h_eleFull5x5E1x5dvE5x5%d",i), "eleFull5x5E1x5dvE5x5", nBin, 0, 1);
    h_eleMissHits[i]             = new TH1D(Form("h_eleMissHits%d",i), "eleMissHits", 6, -0.5, 5.5);
    h_eleD0[i]                   = new TH1D(Form("h_eleD0%d",i), "eleD0", nBin, -0.015, 0.015);  
    h_eleMiniIsoEA[i]            = new TH1D(Form("h_eleMiniIsoEA%d",i), "eleMiniIsoEA", nBin, 0, 0.12);
    h_eventWeight[i]             = new TH1D(Form("h_eventWeight%d",i), "eventWeight", 100, -1, 1);

    h_eleEtaseedAtVtx[i]        ->Sumw2();
    h_eledPhiAtVtx[i]           ->Sumw2();
    h_eleHoverE[i]              ->Sumw2();
    h_eleSigmaIEtaIEtaFull5x5[i]->Sumw2();
    h_eleFull5x5E2x5dvE5x5[i]   ->Sumw2();
    h_eleFull5x5E1x5dvE5x5[i]   ->Sumw2();
    h_eleMissHits[i]            ->Sumw2(); 
    h_eleD0[i]                  ->Sumw2();
    h_eleMiniIsoEA[i]           ->Sumw2();

    h_eleEtaseedAtVtx[i]        ->GetXaxis()->SetTitle("eleEtaseedAtVtx"); 
    h_eledPhiAtVtx[i]           ->GetXaxis()->SetTitle("eledPhiAtVtx");   
    h_eleHoverE[i]              ->GetXaxis()->SetTitle("eleHoverE");    
    h_eleSigmaIEtaIEtaFull5x5[i]->GetXaxis()->SetTitle("eleSigmaIEtaIEtaFull5x5");  
    h_eleFull5x5E2x5dvE5x5[i]   ->GetXaxis()->SetTitle("eleFull5x5E2x5dvE5x5");
    h_eleFull5x5E1x5dvE5x5[i]   ->GetXaxis()->SetTitle("eleFull5x5E1x5dvE5x5");   
    h_eleMissHits[i]            ->GetXaxis()->SetTitle("eleMissHits"); 
    h_eleD0[i]                  ->GetXaxis()->SetTitle("eleD0");      
    h_eleMiniIsoEA[i]           ->GetXaxis()->SetTitle("eleMiniIsoEA"); 
    h_eventWeight[i]            ->GetXaxis()->SetTitle("eventWeight");

  }

  // begin of event loop

  for( Long64_t ev = 0; ev < data.GetEntriesFast(); ev++ ){

    if( ev % 1000000 == 0 )
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());

    data.GetEntry(ev);

    Int_t    nVtx                    = data.GetInt("nVtx");
    Int_t    nEle                    = data.GetInt("nEle");
    Int_t*   eleMissHits             = data.GetPtrInt("eleMissHits");
    Float_t  mcWeight                = data.GetFloat("mcWeight");
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

    Double_t eventWeight = mcWeight;
    if( inputFile.find("DYJets") != std::string::npos ){
      if( eventWeight > 0 ) eventWeight = 1;
      else if( eventWeight < 0 ) eventWeight = -1;
    }
    else
      eventWeight = 1;
    
    h_eventWeight[0]->Fill(0.,eventWeight);
    h_eventWeight[1]->Fill(0.,eventWeight);

    if( nVtx < 1 ) continue;

    // data trigger cut (electron channel)

    std::string* trigName = data.GetPtrString("hlt_trigName");
    vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
    const Int_t nsize = data.GetPtrStringSize();    
    bool passTrigger = false;
    
    for(Int_t it = 0; it < nsize; it++){
    
      std::string thisTrig = trigName[it];
      bool results = trigResult[it];
      
      if( thisTrig.find("HLT_Ele105") != std::string::npos && results==1 ){
	passTrigger = true;
	break;
      }
      
    }

    if( !passTrigger ) continue;

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

      if( fabs(eleScEta[eleId[ie]]) < 1.4442 ){ // barrel selections and cuts

	for(Int_t flag = 0; flag <= 7; flag++){

	  if( fabs(eleEtaseedAtVtx[eleId[ie]])   >= 0.004     && flag != 0 ) continue;
	  if( fabs(eledPhiAtVtx[eleId[ie]])      >= 0.06      && flag != 1 ) continue;
	  if( eleHoverE[eleId[ie]]               >= 1/E+0.05  && flag != 2 ) continue;
	  if( eleMissHits[eleId[ie]]             >  1         && flag != 3 ) continue;
	  if( fabs(eleD0[eleId[ie]])             >= 0.02      && flag != 4 ) continue;
	  if( eleMiniIsoEA[eleId[ie]]            >= 0.1       && flag != 5 ) continue;
	  if( (eleE2x5Full5x5[eleId[ie]]/eleE5x5Full5x5[eleId[ie]] <= 0.94  &&
	       eleE1x5Full5x5[eleId[ie]]/eleE5x5Full5x5[eleId[ie]] <= 0.83) && flag != 6 )
	    continue;
	  	    
	  switch(flag){

	  case 0: h_eleEtaseedAtVtx[0] ->Fill(eleEtaseedAtVtx[eleId[ie]],eventWeight); break;
	  case 1: h_eledPhiAtVtx[0]    ->Fill(eledPhiAtVtx[eleId[ie]],eventWeight);    break;
	  case 2: h_eleHoverE[0]       ->Fill(eleHoverE[eleId[ie]],eventWeight);       break;
	  case 3: h_eleMissHits[0]     ->Fill(eleMissHits[eleId[ie]],eventWeight);     break;
	  case 4: h_eleD0[0]           ->Fill(eleD0[eleId[ie]],eventWeight);           break;
	  case 5: h_eleMiniIsoEA[0]    ->Fill(eleMiniIsoEA[eleId[ie]],eventWeight);    break;	 
	  case 6:
	    h_eleFull5x5E2x5dvE5x5[0]  ->Fill(eleE2x5Full5x5[eleId[ie]]/eleE5x5Full5x5[eleId[ie]],eventWeight);
	    h_eleFull5x5E1x5dvE5x5[0]  ->Fill(eleE1x5Full5x5[eleId[ie]]/eleE5x5Full5x5[eleId[ie]],eventWeight);
	    break;
	  case 7: h_eleSigmaIEtaIEtaFull5x5[0] ->Fill(eleSigmaIEtaIEtaFull5x5[eleId[ie]],eventWeight); break;
	    
	  } // end of switch
    
	} // end of flag loop

      } // end of barrel


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

	  case 0: h_eleEtaseedAtVtx[1] ->Fill(eleEtaseedAtVtx[eleId[ie]],eventWeight); break;
	  case 1: h_eledPhiAtVtx[1]    ->Fill(eledPhiAtVtx[eleId[ie]],eventWeight);    break;
	  case 2: h_eleHoverE[1]       ->Fill(eleHoverE[eleId[ie]],eventWeight);       break;
	  case 3: h_eleSigmaIEtaIEtaFull5x5[1] ->Fill(eleSigmaIEtaIEtaFull5x5[eleId[ie]],eventWeight); break;
	  case 4: h_eleMissHits[1]     ->Fill(eleMissHits[eleId[ie]],eventWeight);     break;
	  case 5: h_eleD0[1]           ->Fill(eleD0[eleId[ie]],eventWeight);           break;
	  case 6: h_eleMiniIsoEA[1]    ->Fill(eleMiniIsoEA[eleId[ie]],eventWeight);    break;	  
	  case 7: 
	    h_eleFull5x5E2x5dvE5x5[1]  ->Fill(eleE2x5Full5x5[eleId[ie]]/eleE5x5Full5x5[eleId[ie]],eventWeight);
	    h_eleFull5x5E1x5dvE5x5[1]  ->Fill(eleE1x5Full5x5[eleId[ie]]/eleE5x5Full5x5[eleId[ie]],eventWeight);
	    break;
		    
	  } // end of switch
	
	} // end of flag loop

      } // end of endcap
	
    } // end of two electrons loop
      
  } // end of event loop

  fprintf(stderr, "Processed all events\n");
    
  TFile* outFile[2];
  
  std::string h_name[10] = {"eleEtaseedAtVtx","eledPhiAtVtx","eleHoverE","eleSigmaIEtaIEtaFull5x5",
			    "eleFull5x5E2x5dvE5x5","eleFull5x5E1x5dvE5x5","eleMissHits","eleD0",
			    "eleMiniIsoEA","eventWeight"};

  std::string region[2] = {"barrel","endcap"};

  for(Int_t i = 0; i < 2; i++){

    outFile[i] = new TFile(Form("%s_%s.root",outputFile.c_str(),region[i].c_str()), "recreate");
      
    h_eleEtaseedAtVtx[i]        ->Write(h_name[0].data());  
    h_eledPhiAtVtx[i]           ->Write(h_name[1].data());  
    h_eleHoverE[i]              ->Write(h_name[2].data());    
    h_eleSigmaIEtaIEtaFull5x5[i]->Write(h_name[3].data()); 
    h_eleFull5x5E2x5dvE5x5[i]   ->Write(h_name[4].data());   
    h_eleFull5x5E1x5dvE5x5[i]   ->Write(h_name[5].data());    
    h_eleMissHits[i]            ->Write(h_name[6].data());   
    h_eleD0[i]                  ->Write(h_name[7].data());      
    h_eleMiniIsoEA[i]           ->Write(h_name[8].data());
    h_eventWeight[i]            ->Write(h_name[9].data());

    outFile[i]->Write();

  }
  
}
