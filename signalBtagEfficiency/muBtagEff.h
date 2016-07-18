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
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZmumu.h"

TGraphAsymmErrors* muBtagEff(string flavor){
  
  // read the ntuples (in pcncu)
  
  TreeReader data("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_13TeV-madgraph.root");

  TH1F* h_jetPtnoCSV = new TH1F("h_jetPtnoCSV", "", 100, 0, 2000);
  TH1F* h_jetPtwtCSV = new TH1F("h_jetPtwtCSV", "", 100, 0, 2000);

  // begin of event loop

  fprintf(stdout, "Total events %lli\n", data.GetEntriesFast());

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    if( (unsigned)ev % 50000 == 0 )
      fprintf(stdout, "Still left events %lli\n", ev);

    data.GetEntry(ev);

    TClonesArray*  muP4              = (TClonesArray*) data.GetPtrTObject("muP4");
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

    if( !isPassZmumu(data,goodLepID) ) continue;

    TLorentzVector* thisLep = (TLorentzVector*)muP4->At(goodLepID[0]);
    TLorentzVector* thatLep = (TLorentzVector*)muP4->At(goodLepID[1]);

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

      if( flavor == "udsg" && 
	  ( FATsubjetFlavor[goodFATJetID][is] != 1  || 
	    FATsubjetFlavor[goodFATJetID][is] != 2  || 
	    FATsubjetFlavor[goodFATJetID][is] != 3  || 
	    FATsubjetFlavor[goodFATJetID][is] != 21 )) continue;
      
      if( flavor == "c" && FATsubjetFlavor[goodFATJetID][is] != 4 ) continue;
      if( flavor == "b" && FATsubjetFlavor[goodFATJetID][is] != 5 ) continue;
      
      TLorentzVector thisSubJet;
      thisSubJet.SetPxPyPzE(FATsubjetSDPx[goodFATJetID][is],
			    FATsubjetSDPy[goodFATJetID][is],
			    FATsubjetSDPz[goodFATJetID][is],
			    FATsubjetSDE[goodFATJetID][is]);
 
      h_jetPtnoCSV->Fill(thisSubJet.Pt());     
    
      if( FATsubjetSDCSV[goodFATJetID][is] > 0.605 )
	h_jetPtwtCSV->Fill(thisSubJet.Pt());
 
    } // end of subjet loop                           

  } // end of event loop
  
  fprintf(stdout, "Processed all events\n");

  // Divide two histograms to get the efficiency

  TGraphAsymmErrors* g_bTagEff = new TGraphAsymmErrors();

  g_bTagEff->Divide(h_jetPtwtCSV, h_jetPtnoCSV, "B");
  g_bTagEff->SetMarkerStyle(8);
  g_bTagEff->SetMinimum(0);
  g_bTagEff->SetMaximum(1.3);
  g_bTagEff->GetYaxis()->SetTitle("Efficiency");  
  g_bTagEff->GetXaxis()->SetTitle("p_{T SubJet} [GeV]");

  delete h_jetPtwtCSV;
  delete h_jetPtnoCSV;

  return g_bTagEff;

}
