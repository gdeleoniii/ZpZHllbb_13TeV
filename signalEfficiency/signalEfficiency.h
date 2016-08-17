R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagCalhead/BTagCalibrationStandalone_cpp.so)
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/untuplizer.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZee.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZmumu.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassJet.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagCalhead/BTagCalibrationStandalone.h"

float signalEfficiency(string inputFile, string channel, int cat, int mzh){

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
  TFile* f_c = TFile::Open(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagEffroot/%s_cflavor_zjetsBtagEff.root", channel.data()));
  TFile* f_b = TFile::Open(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagEffroot/%s_bflavor_signalBtagEff.root", channel.data()));
  
  TGraphAsymmErrors* g_l = (TGraphAsymmErrors*)(f_l->Get(Form("%s_udsgflavor", channel.data())));
  TGraphAsymmErrors* g_c = (TGraphAsymmErrors*)(f_c->Get(Form("%s_cflavor", channel.data())));
  TGraphAsymmErrors* g_b = (TGraphAsymmErrors*)(f_b->Get(Form("%s_bflavor_m%i", channel.data(), mzh)));

  // read the ntuples (in pcncu)

  TreeReader data(inputFile.data());
  
  TFile f(inputFile.data());

  float totalEvent = ((TH1D*)f.Get("h_totalEv"))->Integral();
  float passEvent  = 0.;

  // begin of event loop

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    data.GetEntry(ev);

    Float_t       eventWeight = data.GetFloat("ev_weight");
    TClonesArray* muP4        = (TClonesArray*) data.GetPtrTObject("muP4");
    TClonesArray* eleP4       = (TClonesArray*) data.GetPtrTObject("eleP4");
    TClonesArray* FATjetP4    = (TClonesArray*) data.GetPtrTObject("FATjetP4");

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

    if( (*thisLep+*thatLep+*thisJet).M() < 750 ) continue;

    // b-tag cut

    int nsubBjet = 0;

    float btagWeight = 1; bTagWeight(data, goodFATJetID, &nsubBjet, reader_l, reader_c, reader_b, g_l, g_c, g_b);

    if( cat == 1 && nsubBjet != 1 ) continue;
    if( cat == 2 && nsubBjet != 2 ) continue;

    passEvent += eventWeight * btagWeight;

  } // end of event loop
  
  return passEvent/totalEvent;

}
