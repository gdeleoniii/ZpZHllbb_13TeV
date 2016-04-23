#include <vector>
#include <string>
#include <iostream>
#include <TH1D.h>
#include <TGraph.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "../untuplizer.h"
#include "../isPassZee.h"

float getEfficiency(string inputFile){

  // read the ntuples (in pcncu)

  TreeReader data(inputFile.data());
  
  int passEvent = 0;

  // begin of event loop

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    data.GetEntry(ev);

    TClonesArray*  eleP4             = (TClonesArray*) data.GetPtrTObject("eleP4");
    Int_t          FATnJet           = data.GetInt("FATnJet");    
    Int_t*         FATnSubSDJet      = data.GetPtrInt("FATnSubSDJet");
    Float_t*       FATjetPRmassCorr  = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    TClonesArray*  FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    vector<bool>&  FATjetPassIDLoose = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    vector<float>* FATsubjetSDCSV    = data.GetPtrVectorFloat("FATsubjetSDCSV", FATnJet);
    vector<float>* FATsubjetSDPx     = data.GetPtrVectorFloat("FATsubjetSDPx", FATnJet);
    vector<float>* FATsubjetSDPy     = data.GetPtrVectorFloat("FATsubjetSDPy", FATnJet);
    vector<float>* FATsubjetSDPz     = data.GetPtrVectorFloat("FATsubjetSDPz", FATnJet);
    vector<float>* FATsubjetSDE      = data.GetPtrVectorFloat("FATsubjetSDE", FATnJet);

    // select good reco level events     
    // select good leptons
      
    vector<int> goodLepID;

    if( !isPassZee(data,goodLepID) ) continue;

    TLorentzVector* thisLep = (TLorentzVector*)eleP4->At(goodLepID[0]);
    TLorentzVector* thatLep = (TLorentzVector*)eleP4->At(goodLepID[1]);

    // select good FATjet

    int goodFATJetID = -1;
    TLorentzVector* thisJet = NULL;

    for( int ij = 0; ij < FATnJet; ++ij ){

      thisJet = (TLorentzVector*)FATjetP4->At(ij);

      if( thisJet->Pt() < 200 ) continue;
      if( fabs(thisJet->Eta()) > 2.4 ) continue;
      if( !FATjetPassIDLoose[ij] ) continue;
      if( thisJet->DeltaR(*thisLep) < 0.8 || thisJet->DeltaR(*thatLep) < 0.8 ) continue;
      if( FATjetPRmassCorr[ij] < 105 || FATjetPRmassCorr[ij] > 135 ) continue;

      int nsubBjet = 0;

      for( int is = 0; is < FATnSubSDJet[ij]; ++is ){

	if( FATsubjetSDCSV[ij][is] > 0.605 ) ++nsubBjet;

      }
 
      Float_t subjetDeltaR = -1;
 
      if( nsubBjet == 2 ){
 
	TLorentzVector l4_subjet0(0,0,0,0);
	TLorentzVector l4_subjet1(0,0,0,0);
 
	l4_subjet0.SetPxPyPzE(FATsubjetSDPx[ij][0],
			      FATsubjetSDPy[ij][0],
			      FATsubjetSDPz[ij][0],
			      FATsubjetSDE [ij][0]);
 
	l4_subjet1.SetPxPyPzE(FATsubjetSDPx[ij][1],
			      FATsubjetSDPy[ij][1],
			      FATsubjetSDPz[ij][1],
			      FATsubjetSDE [ij][1]);
 
	subjetDeltaR = l4_subjet0.DeltaR(l4_subjet1);
 
      }
 
      // deltaR depends b-tag cut
 
      if( !(subjetDeltaR < 0.3 && nsubBjet > 0) ) continue;
      //if( !(subjetDeltaR > 0.3 && nsubBjet > 1) ) continue;
       
      goodFATJetID = ij;

      break;
 
    } // end of FatnJet loop
 
    if( goodFATJetID < 0 ) continue;

    if( (*thisLep+*thatLep+*thisJet).M() < 750 ) continue;

    ++passEvent;

  } // end of event loop
  
  return (float)passEvent/(float)data.GetEntriesFast();

}

void signalEfficiency(){

  Float_t x_mzh[13] = {600,800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500};
  Float_t y_eff[13];

  for( int i = 0; i < 13; ++i )
    y_eff[i] = getEfficiency(Form("/data7/htong/skim_samples/ele/skim_ele_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",(Int_t)x_mzh[i]));

  TGraph g_eff(13, x_mzh, y_eff);

  TCanvas c("c", "", 0, 0, 800, 600);

  c.cd();
  g_eff.Draw();
  c.Print("eleIDeff.pdf");

}
