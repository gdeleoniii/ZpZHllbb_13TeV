#include <vector>
#include <string>
#include <iostream>
#include <TH1F.h>
#include <TH2F.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/untuplizer.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZee.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZmumu.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassJet.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/leptonWeight.h"

TGraphAsymmErrors* btaggingEff(string inputFile, string channel){
  
  float varBins[] = {30,50,70,100,140,200,300,670,2000};
  int   nvarBins  = sizeof(varBins)/sizeof(varBins[0])-1;

  TH1F* h_jetPtnoCSV = new TH1F("h_jetPtnoCSV", "", nvarBins, varBins);
  TH1F* h_jetPtwtCSV = new TH1F("h_jetPtwtCSV", "", nvarBins, varBins);

  // to read lepton scale factor / trigger

  TFile* f_ele    = TFile::Open("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/leptonSFroot/CutBasedID_LooseWP_fromTemplates_withSyst_Final.txt_SF2D.root");
  TFile* f_muScal = TFile::Open("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/leptonSFroot/MuonHighPt_Z_RunCD_Reco74X_Dec17.root");
  TFile* f_muTrig = TFile::Open("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/leptonSFroot/SingleMuonTrigger_Z_RunCD_Reco74X_Dec1.root");

  TH2F* h2_ele    = (TH2F*)(f_ele->Get("EGamma_SF2D"));
  TH2F* h2_muPt20 = (TH2F*)(f_muScal->Get("HighPtID_PtEtaBins_Pt20/abseta_pTtuneP_ratio"));
  TH2F* h2_muPt53 = (TH2F*)(f_muScal->Get("HighPtID_PtEtaBins_Pt53/abseta_pTtuneP_ratio"));
  TH2F* h2_muRunD = (TH2F*)(f_muTrig->Get("runD_Mu45_eta2p1_PtEtaBins/abseta_pt_ratio"));

  // read the ntuples (in pcncu)

  TreeReader data(inputFile.data());

  // begin of event loop

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    data.GetEntry(ev);

    Float_t        eventWeight     = data.GetFloat("ev_weight");
    TClonesArray*  muP4            = (TClonesArray*) data.GetPtrTObject("muP4");
    TClonesArray*  eleP4           = (TClonesArray*) data.GetPtrTObject("eleP4");
    TClonesArray*  FATjetP4        = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    Int_t          FATnJet         = data.GetInt("FATnJet");    
    Int_t*         FATnSubSDJet    = data.GetPtrInt("FATnSubSDJet");
    vector<float>* FATsubjetSDCSV  = data.GetPtrVectorFloat("FATsubjetSDCSV", FATnJet);
    vector<float>* FATsubjetSDPx   = data.GetPtrVectorFloat("FATsubjetSDPx", FATnJet);
    vector<float>* FATsubjetSDPy   = data.GetPtrVectorFloat("FATsubjetSDPy", FATnJet);
    vector<float>* FATsubjetSDPz   = data.GetPtrVectorFloat("FATsubjetSDPz", FATnJet);
    vector<float>* FATsubjetSDE    = data.GetPtrVectorFloat("FATsubjetSDE", FATnJet);
    vector<int>*   FATsubjetFlavor = data.GetPtrVectorInt("FATsubjetSDHadronFlavor", FATnJet);
    vector<bool>&  isHighPtMuon    = *((vector<bool>*) data.GetPtr("isHighPtMuon"));

    // select good reco level events     
    // select good leptons
      
    vector<int> goodLepID;

    if( channel == "ele" && !isPassZee(data,goodLepID)   ) continue;
    if( channel == "mu"  && !isPassZmumu(data,goodLepID) ) continue;

    TLorentzVector* thisLep = (channel=="ele") ? (TLorentzVector*)eleP4->At(goodLepID[0]) : (TLorentzVector*)muP4->At(goodLepID[0]);
    TLorentzVector* thatLep = (channel=="ele") ? (TLorentzVector*)eleP4->At(goodLepID[1]) : (TLorentzVector*)muP4->At(goodLepID[1]);

    // calculate lepton weight
    
    float thisLepWeight=1, thatLepWeight=1;

    if( channel == "ele" ){
    
      thisLepWeight = leptonWeight(h2_ele, thisLep, false);
      thatLepWeight = leptonWeight(h2_ele, thatLep, false);

    }

    else if( channel == "mu" ){

      if( isHighPtMuon[goodLepID[0]] ) 
	thisLepWeight = (thisLep->Pt() < 53) ? leptonWeight(h2_muPt20, thisLep) : leptonWeight(h2_muPt53, thisLep);
      
      else thisLepWeight = 1;
      
      if( isHighPtMuon[goodLepID[1]] )
	thatLepWeight = (thatLep->Pt() < 53) ? leptonWeight(h2_muPt20, thatLep) : leptonWeight(h2_muPt53, thatLep);
      
      else thatLepWeight = 1;
      
    }

    // calculate trigger weight for muon

    float muTrigWeight=1;

    if( channel=="mu" ){

      if( thisLep->Pt() > thatLep->Pt() )
	muTrigWeight = (fabs(thisLep->Eta()) < 2.1) ? leptonWeight(h2_muRunD, thisLep) : 1;
      else 
	muTrigWeight = (fabs(thatLep->Eta()) < 2.1) ? leptonWeight(h2_muRunD, thatLep) : 1;

    }

    else muTrigWeight = 1;

    // select good FATjet

    int goodFATJetID = -1;

    if( !isPassJet(data, &goodFATJetID, thisLep, thatLep) ) continue;

    TLorentzVector* thisJet = (TLorentzVector*)FATjetP4->At(goodFATJetID);

    float mllbb;

    noiseCleaning(&mllbb, thisLep, thatLep, thisJet);

    // b-tag efficiency part (b flavor only)

    for( int is = 0; is < FATnSubSDJet[goodFATJetID]; ++is ){
      
      if( FATsubjetFlavor[goodFATJetID][is] != 5 ) continue;
 
      TLorentzVector thisSubJet;

      thisSubJet.SetPxPyPzE(FATsubjetSDPx[goodFATJetID][is],
			    FATsubjetSDPy[goodFATJetID][is],
			    FATsubjetSDPz[goodFATJetID][is],
			    FATsubjetSDE[goodFATJetID][is]);
 
      h_jetPtnoCSV->Fill(thisSubJet.Pt(), eventWeight*thisLepWeight*thatLepWeight*muTrigWeight);     
    
      if( FATsubjetSDCSV[goodFATJetID][is] > 0.605 )
	h_jetPtwtCSV->Fill(thisSubJet.Pt(), eventWeight*thisLepWeight*thatLepWeight*muTrigWeight);
 
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
