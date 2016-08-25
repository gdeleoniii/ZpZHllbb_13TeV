#include <vector>
#include <string>
#include <iostream>
#include <TH1F.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/untuplizer.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZee.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZmumu.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassJet.h"

TGraphAsymmErrors* bTagEff(string inputFile, string channel, string flavor){
  
  // read the ntuples (in pcncu)
  
  TreeReader data(inputFile.data());

  float varBins[] = {30,50,70,100,140,200,300,670,2000};
  int   nvarBins  = sizeof(varBins)/sizeof(varBins[0])-1;

  TH1F* h_jetPtnoCSV = new TH1F("h_jetPtnoCSV", "", nvarBins, varBins);
  TH1F* h_jetPtwtCSV = new TH1F("h_jetPtwtCSV", "", nvarBins, varBins);

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

    // select good reco level events     
    // select good leptons
      
    vector<int> goodLepID;

    if( channel == "ele" && !isPassZee(data,goodLepID)   ) continue;
    if( channel == "mu"  && !isPassZmumu(data,goodLepID) ) continue;

    TLorentzVector* thisLep = (channel=="ele") ? (TLorentzVector*)eleP4->At(goodLepID[0]) : (TLorentzVector*)muP4->At(goodLepID[0]);
    TLorentzVector* thatLep = (channel=="ele") ? (TLorentzVector*)eleP4->At(goodLepID[1]) : (TLorentzVector*)muP4->At(goodLepID[1]);

    // select good FATjet

    int goodFATJetID = -1;

    if( !isPassJet(data, &goodFATJetID, thisLep, thatLep) ) continue;

    TLorentzVector* thisJet = (TLorentzVector*)FATjetP4->At(goodFATJetID);

    float mllbb;

    noiseCleaning(&mllbb, thisLep, thatLep, thisJet);

    // b-tag efficiency part

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
