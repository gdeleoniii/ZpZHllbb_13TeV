#include <vector>
#include <string>
#include <iostream>
#include <TH1.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/untuplizer.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZee.h"

TGraphAsymmErrors* eleBtagEff(string inputFile, int cat){
  
  // read the ntuples (in pcncu)
  
  TreeReader data(inputFile.data());

  TH1F* h_jetPtnoCSV = new TH1F("h_jetPtnoCSV", "", 100, 0, 3000);
  TH1F* h_jetPtwtCSV = new TH1F("h_jetPtwtCSV", "", 100, 0, 3000);

  // begin of event loop

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    data.GetEntry(ev);

    TClonesArray*  eleP4              = (TClonesArray*) data.GetPtrTObject("eleP4");
    Int_t          FATnJet           = data.GetInt("FATnJet");    
    Int_t*         FATnSubSDJet      = data.GetPtrInt("FATnSubSDJet");
    Float_t*       FATjetPRmassCorr  = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    TClonesArray*  FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    vector<bool>&  FATjetPassIDLoose = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    vector<float>* FATsubjetSDCSV    = data.GetPtrVectorFloat("FATsubjetSDCSV", FATnJet);

    // select good reco level events     
    // select good leptons
      
    vector<int> goodLepID;

    if( !isPassZee(data,goodLepID) ) continue;

    TLorentzVector* thisLep = (TLorentzVector*)eleP4->At(goodLepID[0]);
    TLorentzVector* thatLep = (TLorentzVector*)eleP4->At(goodLepID[1]);

    // select good FATjet

    bool bTag = false;
    int goodFATJetID = -1;
    TLorentzVector thisJet(0,0,0,0);

    for( int ij = 0; ij < FATnJet; ++ij ){

      TLorentzVector* myJet = (TLorentzVector*)FATjetP4->At(ij);

      if( myJet->Pt() < 200 ) continue;
      if( fabs(myJet->Eta()) > 2.4 ) continue;
      if( !FATjetPassIDLoose[ij] ) continue;
      if( myJet->DeltaR(*thisLep) < 0.8 || myJet->DeltaR(*thatLep) < 0.8 ) continue;
      if( FATjetPRmassCorr[ij] < 105 || FATjetPRmassCorr[ij] > 135 ) continue;

      int nsubBjet = 0;

      for( int is = 0; is < FATnSubSDJet[ij]; ++is ){

	if( FATsubjetSDCSV[ij][is] > 0.605 ) ++nsubBjet;

      } // end of subjet for loop
 
      // b-tag cut
 
      if( cat == 1 && nsubBjet == 1 ) bTag = true;
      if( cat == 2 && nsubBjet == 2 ) bTag = true;
       
      goodFATJetID = ij;
      thisJet = *myJet;

      break;
 
    } // end of FatnJet loop
 
    if( goodFATJetID < 0 ) continue;

    if( (*thisLep+*thatLep+thisJet).M() < 750 ) continue;

    if( bTag ) h_jetPtwtCSV->Fill(thisJet.Pt());

    h_jetPtnoCSV->Fill(thisJet.Pt());

  } // end of event loop
  
  // Divide two histograms to get the efficiency

  TGraphAsymmErrors* g_bTagEff = new TGraphAsymmErrors();

  g_bTagEff->Divide(h_jetPtwtCSV, h_jetPtnoCSV, "B");
  g_bTagEff->SetMarkerStyle(8);
  g_bTagEff->SetMaximum(1.3);
  g_bTagEff->GetYaxis()->SetTitle("Efficiency");  
  g_bTagEff->GetXaxis()->SetTitle("p_{T jet} [GeV]");

  delete h_jetPtwtCSV;
  delete h_jetPtnoCSV;

  return g_bTagEff;

}
