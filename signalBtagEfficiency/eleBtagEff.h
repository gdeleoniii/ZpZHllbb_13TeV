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

TGraphAsymmErrors* eleBtagEff(string inputFile, string flavor){
  
  // read the ntuples (in pcncu)
  
  TreeReader data(inputFile.data());

  float varBins[] = {30,50,70,100,140,200,300,670,2000};
  int   nvarBins  = sizeof(varBins)/sizeof(varBins[0])-1;

  TH1F* h_jetPtnoCSV = new TH1F("h_jetPtnoCSV", "", nvarBins, varBins);
  TH1F* h_jetPtwtCSV = new TH1F("h_jetPtwtCSV", "", nvarBins, varBins);

  // begin of event loop

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    data.GetEntry(ev);

    Float_t        eventWeight       = data.GetFloat("ev_weight");
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
    vector<int>*   FATsubjetFlavor   = data.GetPtrVectorInt("FATsubjetSDHadronFlavor", FATnJet);

    // select good reco level events     
    // select good leptons
      
    vector<int> goodLepID;

    if( !isPassZee(data,goodLepID) ) continue;

    TLorentzVector* thisLep = (TLorentzVector*)eleP4->At(goodLepID[0]);
    TLorentzVector* thatLep = (TLorentzVector*)eleP4->At(goodLepID[1]);

    // select good FATjet

    int goodFATJetID = -1;
    TLorentzVector thisJet(0,0,0,0);

    for( int ij = 0; ij < FATnJet; ++ij ){

      TLorentzVector* myJet = (TLorentzVector*)FATjetP4->At(ij);

      if( myJet->Pt() < 200 ) continue;
      if( fabs(myJet->Eta()) > 2.4 ) continue;
      if( !FATjetPassIDLoose[ij] ) continue;
      if( myJet->DeltaR(*thisLep) < 0.8 || myJet->DeltaR(*thatLep) < 0.8 ) continue;
      if( FATjetPRmassCorr[ij] < 105 || FATjetPRmassCorr[ij] > 135 ) continue;
        
      goodFATJetID = ij;
      thisJet = *myJet;

      break;
 
    } // end of FatnJet loop
 
    if( goodFATJetID < 0 ) continue;

    if( (*thisLep+*thatLep+thisJet).M() < 750 ) continue;

    for( int is = 0; is < FATnSubSDJet[goodFATJetID]; ++is ){
      
      if( flavor == "udsg" && (FATsubjetFlavor[goodFATJetID][is] == 4 || FATsubjetFlavor[goodFATJetID][is] == 5) ) continue;
      if( flavor == "c"    && FATsubjetFlavor[goodFATJetID][is] != 4 ) continue;
      if( flavor == "b"    && FATsubjetFlavor[goodFATJetID][is] != 5 ) continue;
 
      TLorentzVector thisSubJet;

      thisSubJet.SetPxPyPzE(FATsubjetSDPx[goodFATJetID][is],
			    FATsubjetSDPy[goodFATJetID][is],
			    FATsubjetSDPz[goodFATJetID][is],
			    FATsubjetSDE[goodFATJetID][is]);
 
      h_jetPtnoCSV->Fill(thisSubJet.Pt(),eventWeight);     
    
      if( FATsubjetSDCSV[goodFATJetID][is] > 0.605 )
	h_jetPtwtCSV->Fill(thisSubJet.Pt(),eventWeight);
 
    } // end of subjet loop                           

  } // end of event loop
  
  // Divide two histograms to get the efficiency

  TGraphAsymmErrors* g_bTagEff = new TGraphAsymmErrors();

  g_bTagEff->BayesDivide(h_jetPtwtCSV, h_jetPtnoCSV, "B");
  g_bTagEff->SetMarkerStyle(8);
  g_bTagEff->SetMinimum(0);
  g_bTagEff->SetMaximum(1.3);
  g_bTagEff->GetYaxis()->SetTitle("Efficiency");  
  g_bTagEff->GetXaxis()->SetTitle("p_{T SubJet} [GeV]");

  delete h_jetPtwtCSV;
  delete h_jetPtnoCSV;

  return g_bTagEff;

}
