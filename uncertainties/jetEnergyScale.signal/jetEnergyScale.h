R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagCalhead/BTagCalibrationStandalone_cpp.so)
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/untuplizer.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZee.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZmumu.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassJet.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagCalhead/BTagCalibrationStandalone.h"

float jetEnergyScale(string inputFile, string js, string channel, int cat, int mzh){

  // setup calibration and reader

  BTagCalibration calib("csvv1", "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/CSVV1.csv");

  BTagCalibrationReader reader_l(BTagEntry::OP_LOOSE, "central");
  BTagCalibrationReader reader_c(BTagEntry::OP_LOOSE, "central");
  BTagCalibrationReader reader_b(BTagEntry::OP_LOOSE, "central");

  reader_l.load(calib, BTagEntry::FLAV_UDSG, "comb");
  reader_c.load(calib, BTagEntry::FLAV_C, "mujets");
  reader_b.load(calib, BTagEntry::FLAV_B, "mujets");

  // to read b-tag effinciency 

  TFile* f_l = TFile::Open(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagEffroot/%s_udsgflavor_zjetsBtagEff.root", channel.data()));
  TFile* f_c = TFile::Open(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagEffroot/%s_cflavor_zjetsBtagEff.root",    channel.data()));
  TFile* f_b = TFile::Open(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagEffroot/%s_bflavor_signalBtagEff.root",   channel.data()));
  
  TH1F* h_l = (TH1F*)(f_l->Get(Form("%s_udsgflavor",  channel.data())));
  TH1F* h_c = (TH1F*)(f_c->Get(Form("%s_cflavor",     channel.data())));
  TH1F* h_b = (TH1F*)(f_b->Get(Form("%s_bflavor_m%i", channel.data(), mzh)));

  // read the ntuples (in pcncu)

  TreeReader data(inputFile.data());
  
  TFile f(inputFile.data());

  float totalEvent = ((TH1D*)f.Get("h_totalEv"))->Integral();
  float passEvent  = 0.;

  // begin of event loop

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    data.GetEntry(ev);

    Float_t       eventWeight       = data.GetFloat("ev_weight");
    TClonesArray* muP4              = (TClonesArray*) data.GetPtrTObject("muP4");
    TClonesArray* eleP4             = (TClonesArray*) data.GetPtrTObject("eleP4");
    TClonesArray* FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    Int_t         FATnJet           = data.GetInt("FATnJet");    
    Float_t*      FATjetPRmassCorr  = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    Float_t*      FATjetCorrUncUp   = data.GetPtrFloat("FATjetCorrUncUp");
    Float_t*      FATjetCorrUncDown = data.GetPtrFloat("FATjetCorrUncDown");
    vector<bool>& FATjetPassIDLoose = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));

    // select good reco level events     
    // select good leptons
      
    vector<int> goodLepID;

    if( channel == "ele" && !isPassZee(data,goodLepID)   ) continue;
    if( channel == "mu"  && !isPassZmumu(data,goodLepID) ) continue;

    TLorentzVector* thisLep = (channel=="ele") ? (TLorentzVector*)eleP4->At(goodLepID[0]) : (TLorentzVector*)muP4->At(goodLepID[0]);
    TLorentzVector* thatLep = (channel=="ele") ? (TLorentzVector*)eleP4->At(goodLepID[1]) : (TLorentzVector*)muP4->At(goodLepID[1]);

    // select good FATjet

    int goodFATJetID = -1;
    TLorentzVector thisJet(0,0,0,0);
    
    for( int ij = 0; ij < FATnJet; ++ij ){

      TLorentzVector* myJet = (TLorentzVector*)FATjetP4->At(ij);

      *myJet *= (js=="central") ? 1 : ( (js=="up") ? (1+FATjetCorrUncUp[ij]) : (1-FATjetCorrUncDown[ij]) );

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

    // b-tag cut

    int nsubBjet = 0;

    float btagWeight = bTagWeight(data, goodFATJetID, &nsubBjet, h_l, h_c, h_b, reader_l, reader_c, reader_b);
    
    if( cat == 1 && nsubBjet != 1 ) continue;
    if( cat == 2 && nsubBjet != 2 ) continue;

    passEvent += eventWeight * btagWeight;

  } // end of event loop
  
  return passEvent/totalEvent;

}
