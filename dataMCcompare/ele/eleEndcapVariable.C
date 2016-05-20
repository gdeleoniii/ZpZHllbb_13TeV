#include <vector>
#include <string>
#include <iostream>
#include <TH1D.h>
#include <TFile.h>
#include <TProfile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "../../untuplizer.h"

float EA(float eta){

  float ea = 0.0;

  if     ( eta > 0.000 && eta < 1.000 ) ea = 0.1752;
  else if( eta > 1.000 && eta < 1.479 ) ea = 0.1862;
  else if( eta > 1.479 && eta < 2.000 ) ea = 0.1411;
  else if( eta > 2.000 && eta < 2.200 ) ea = 0.1534;
  else if( eta > 2.200 && eta < 2.300 ) ea = 0.1903;
  else if( eta > 2.300 && eta < 2.400 ) ea = 0.2243;
  else if( eta > 2.400 && eta < 5.000 ) ea = 0.2687;

  return ea;

}

void eleEndcapVariable(std::string inputFile, std::string outputFile){

  // read the ntuples (in pcncu)

  TreeReader data(inputFile.data());
  
  // Declare the histogram

  TFile* f = new TFile(inputFile.data());
  TH1D* h_totalEvents             = (TH1D*)f->Get("h_totalEv");

  TH1D* h_eleSigmaIEtaIEtaFull5x5 = new TH1D("h_eleSigmaIEtaIEtaFull5x5", "eleSigmaIEtaIEtaFull5x5", 40,      0,  0.06);
  TH1D* h_eleEtaseedAtVtx         = new TH1D("h_eleEtaseedAtVtx",         "eleEtaseedAtVtx",         20,  -0.01,  0.01);
  TH1D* h_eledPhiAtVtx            = new TH1D("h_eledPhiAtVtx",            "eledPhiAtVtx",            20,  -0.03,  0.03);
  TH1D* h_eleHoverE               = new TH1D("h_eleHoverE",               "eleHoverE",               20,      0,  0.05);
  TH1D* h_eleRelIsoWithEA         = new TH1D("h_eleRelIsoWithEA",         "eleRelIsoWithEA",         20,      0,  0.12);
  TH1D* h_eleEoverPInv            = new TH1D("h_eleEoverPInv",            "eleEoverPInv",            20,      0,  0.12);
  TH1D* h_eleD0                   = new TH1D("h_eleD0",                   "eleD0",                   20, -0.015, 0.015);
  TH1D* h_eleDz                   = new TH1D("h_eleDz",                   "eleDz",                   20, -0.015, 0.015);
  TH1D* h_eleMissHits             = new TH1D("h_eleMissHits",             "eleMissHits",              6,   -0.5,   5.5);

  h_eleSigmaIEtaIEtaFull5x5->Sumw2();
  h_eleEtaseedAtVtx        ->Sumw2();
  h_eledPhiAtVtx           ->Sumw2();
  h_eleHoverE              ->Sumw2();
  h_eleRelIsoWithEA        ->Sumw2();
  h_eleEoverPInv           ->Sumw2();
  h_eleD0                  ->Sumw2();
  h_eleDz                  ->Sumw2();
  h_eleMissHits            ->Sumw2(); 

  h_eleSigmaIEtaIEtaFull5x5->GetXaxis()->SetTitle("eleSigmaIEtaIEtaFull5x5");
  h_eleEtaseedAtVtx        ->GetXaxis()->SetTitle("eleEtaseedAtVtx"); 
  h_eledPhiAtVtx           ->GetXaxis()->SetTitle("eledPhiAtVtx");   
  h_eleHoverE              ->GetXaxis()->SetTitle("eleHoverE");    
  h_eleRelIsoWithEA        ->GetXaxis()->SetTitle("eleRelIsoWithEA");
  h_eleEoverPInv           ->GetXaxis()->SetTitle("eleEoverPInv");
  h_eleD0                  ->GetXaxis()->SetTitle("eleD0");
  h_eleDz                  ->GetXaxis()->SetTitle("eleDz");
  h_eleMissHits            ->GetXaxis()->SetTitle("eleMissHits"); 

  // begin of event loop

  fprintf(stdout, "Total events %lli\n", data.GetEntriesFast());

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    if( (unsigned)ev % 100000 == 0 )
      fprintf(stdout, "Still left events %lli\n", ev);

    data.GetEntry(ev);

    Int_t    nEle                    = data.GetInt("nEle");
    Int_t*   eleMissHits             = data.GetPtrInt("eleMissHits");
    Float_t  eventWeight             = data.GetFloat("ev_weight");
    Float_t  eleRho                  = data.GetFloat("eleRho");
    Float_t* eleChHadIso             = data.GetPtrFloat("eleChHadIso");
    Float_t* eleNeHadIso             = data.GetPtrFloat("eleNeHadIso");
    Float_t* eleGamIso               = data.GetPtrFloat("eleGamIso");
    Float_t* eleScEt                 = data.GetPtrFloat("eleScEt");
    Float_t* eleScEta                = data.GetPtrFloat("eleScEta");
    Float_t* eleEtaseedAtVtx         = data.GetPtrFloat("eleEtaseedAtVtx");	
    Float_t* eledPhiAtVtx            = data.GetPtrFloat("eledPhiAtVtx");
    Float_t* eleHoverE               = data.GetPtrFloat("eleHoverE");
    Float_t* eleSigmaIEtaIEtaFull5x5 = data.GetPtrFloat("eleSigmaIEtaIEtaFull5x5");
    Float_t* eleEoverPInv            = data.GetPtrFloat("eleEoverPInv");
    Float_t* eleD0                   = data.GetPtrFloat("eleD0");
    Float_t* eleDz                   = data.GetPtrFloat("eleDz");
    TClonesArray* eleP4              = (TClonesArray*) data.GetPtrTObject("eleP4");
    vector<bool>& eleConvVeto        = *((vector<bool>*) data.GetPtr("eleConvVeto"));
    vector<bool>& eleInBarrel        = *((vector<bool>*) data.GetPtr("eleInBarrel"));
    vector<bool>& eleInEndcap        = *((vector<bool>*) data.GetPtr("eleInEndcap"));

    // choosing electron

    vector<Int_t> eleId;

    bool findEPair = false;
    
    for(Int_t ie = 0; ie < nEle; ++ie){
      for(Int_t je = 0; je < ie; ++je){

	TLorentzVector* thisEle = (TLorentzVector*)eleP4->At(ie);
	TLorentzVector* thatEle = (TLorentzVector*)eleP4->At(je);

	if( !eleInBarrel[ie] && !eleInEndcap[ie] ) continue;
	if( !eleInBarrel[je] && !eleInEndcap[je] ) continue;
	if( !eleConvVeto[ie] || !eleConvVeto[je] ) continue;
	if( eleScEt[ie] < 35 || eleScEt[je] < 35 ) continue;
	if( (*thisEle+*thatEle).M()  < 60  ) continue;
	if( (*thisEle+*thatEle).M()  > 120 ) continue;		
	if( (*thisEle+*thatEle).Pt() < 150 ) continue;

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

    for(Int_t ie = 0; ie < 2; ++ie){

      Float_t eleRelIsoWithEA = eleChHadIso[eleId[ie]] + TMath::Max(0.0f, eleNeHadIso[eleId[ie]]+eleGamIso[eleId[ie]]-eleRho*EA(eleScEta[eleId[ie]]));

      if( eleInEndcap[eleId[ie]] ){ // endcap selections and cuts

	for(Int_t flag = 0; flag <= 8; ++flag){

	  if( eleSigmaIEtaIEtaFull5x5[eleId[ie]] >= 0.0301  && flag != 0 ) continue;
	  if( fabs(eleEtaseedAtVtx[eleId[ie]])   >= 0.00814 && flag != 1 ) continue;
	  if( fabs(eledPhiAtVtx[eleId[ie]])      >= 0.182   && flag != 2 ) continue;
	  if( eleHoverE[eleId[ie]]               >= 0.0897  && flag != 3 ) continue;
	  if( eleRelIsoWithEA                    >= 0.121   && flag != 4 ) continue;
	  if( eleEoverPInv[eleId[ie]]            >= 0.126   && flag != 5 ) continue;
	  if( fabs(eleD0[eleId[ie]])             >= 0.118   && flag != 6 ) continue;
	  if( fabs(eleDz[eleId[ie]])             >= 0.822   && flag != 7 ) continue;
	  if( eleMissHits[eleId[ie]]             >  1       && flag != 8 ) continue;
	  	    
	  switch(flag){

	  case 0: h_eleSigmaIEtaIEtaFull5x5 ->Fill(eleSigmaIEtaIEtaFull5x5[eleId[ie]], eventWeight); break;
	  case 1: h_eleEtaseedAtVtx         ->Fill(eleEtaseedAtVtx[eleId[ie]],         eventWeight); break;
	  case 2: h_eledPhiAtVtx            ->Fill(eledPhiAtVtx[eleId[ie]],            eventWeight); break;
	  case 3: h_eleHoverE               ->Fill(eleHoverE[eleId[ie]],               eventWeight); break;
	  case 4: h_eleRelIsoWithEA         ->Fill(eleRelIsoWithEA,                    eventWeight); break;
	  case 5: h_eleEoverPInv            ->Fill(eleEoverPInv[eleId[ie]],            eventWeight); break;
	  case 6: h_eleD0                   ->Fill(eleD0[eleId[ie]],                   eventWeight); break;
	  case 7: h_eleDz                   ->Fill(eleDz[eleId[ie]],                   eventWeight); break;
	  case 8: h_eleMissHits             ->Fill(eleMissHits[eleId[ie]],             eventWeight); break;
	    
	  } // end of switch
    
	} // end of flag loop

      } // end of endcap
	
    } // end of two electrons loop
      
  } // end of event loop

  fprintf(stdout, "Processed all events\n");
    
  TFile* outFile = new TFile(Form("%s_eleEndcapVariable.root",outputFile.c_str()), "recreate");
      
  h_eleSigmaIEtaIEtaFull5x5->Write("eleSigmaIEtaIEtaFull5x5");
  h_eleEtaseedAtVtx        ->Write("eleEtaseedAtVtx");
  h_eledPhiAtVtx           ->Write("eledPhiAtVtx");
  h_eleHoverE              ->Write("eleHoverE");
  h_eleRelIsoWithEA        ->Write("eleRelIsoWithEA");
  h_eleEoverPInv           ->Write("eleEoverPInv");
  h_eleD0                  ->Write("eleD0");
  h_eleDz                  ->Write("eleDz");
  h_eleMissHits            ->Write("eleMissHits");
  h_totalEvents            ->Write("totalEvents");

  outFile->Write();

  delete outFile;  
  delete f;

}
