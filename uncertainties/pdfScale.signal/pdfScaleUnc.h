R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagCalhead/BTagCalibrationStandalone_cpp.so)
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/untuplizer.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZee.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZmumu.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassJet.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/leptonWeight.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagCalhead/BTagCalibrationStandalone.h"

float pdfScaleUnc(string inputFile, string channel, int cat, int mzh, int first, int last, bool onlyCentral=false){

  // setup calibration and reader

  BTagCalibration calib("csvv1", "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/CSVV1.csv");

  BTagCalibrationReader reader_l(BTagEntry::OP_LOOSE, "central");
  BTagCalibrationReader reader_c(BTagEntry::OP_LOOSE, "central");
  BTagCalibrationReader reader_b(BTagEntry::OP_LOOSE, "central");

  reader_l.load(calib, BTagEntry::FLAV_UDSG, "comb");
  reader_c.load(calib, BTagEntry::FLAV_C,    "mujets");
  reader_b.load(calib, BTagEntry::FLAV_B,    "mujets");

  // to read b-tag effinciency 

  TFile* f_l = TFile::Open(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagEffroot/%s_udsgflavor_zjetsBtagEff.root", channel.data()));
  TFile* f_c = TFile::Open(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagEffroot/%s_cflavor_zjetsBtagEff.root",    channel.data()));
  TFile* f_b = TFile::Open(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagEffroot/%s_bflavor_signalBtagEff.root",   channel.data()));
  
  TH1F* h_l = (TH1F*)(f_l->Get(Form("%s_udsgflavor",  channel.data())));
  TH1F* h_c = (TH1F*)(f_c->Get(Form("%s_cflavor",     channel.data())));
  TH1F* h_b = (TH1F*)(f_b->Get(Form("%s_bflavor_m%i", channel.data(), mzh)));

  // to read lepton scale factor

  TFile* f_ele = TFile::Open("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/leptonSFroot/CutBasedID_LooseWP_fromTemplates_withSyst_Final.txt_SF2D.root");
  TFile* f_mu  = TFile::Open("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/leptonSFroot/MuonHighPt_Z_RunCD_Reco74X_Dec17.root");

  TH2F* h2_ele    = (TH2F*)(f_ele->Get("EGamma_SF2D"));
  TH2F* h2_muPt20 = (TH2F*)(f_mu->Get("HighPtID_PtEtaBins_Pt20/abseta_pTtuneP_ratio"));
  TH2F* h2_muPt53 = (TH2F*)(f_mu->Get("HighPtID_PtEtaBins_Pt53/abseta_pTtuneP_ratio"));

  // read the ntuples (in pcncu)

  TreeReader data(inputFile.data());

  TFile f(inputFile.data());

  int N = 1+last-first;
  float* efficiency = new float[N-1];
  float* passEvent  = new float[N];
  float  totalEvent = ((TH1D*)f.Get("h_totalEv"))->Integral();
  float  efficiencyCentral = -1;

  std::fill_n(efficiency,N,0.);
  std::fill_n(passEvent,N,0.);

  // begin of event loop

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    data.GetEntry(ev);

    Float_t*      pdfscaleSysWeight = data.GetPtrFloat("pdfscaleSysWeights");
    Float_t       eventWeight       = data.GetFloat("ev_weight");
    TClonesArray* muP4              = (TClonesArray*) data.GetPtrTObject("muP4");
    TClonesArray* eleP4             = (TClonesArray*) data.GetPtrTObject("eleP4");
    TClonesArray* FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    vector<bool>& isHighPtMuon = *((vector<bool>*) data.GetPtr("isHighPtMuon"));

    // select good reco level events     
    // select good leptons
      
    vector<int> goodLepID;

    if( channel == "ele" && !isPassZee(data,goodLepID)   ) continue;
    if( channel == "mu"  && !isPassZmumu(data,goodLepID) ) continue;

    TLorentzVector* thisLep = (channel=="ele") ? (TLorentzVector*)eleP4->At(goodLepID[0]) : (TLorentzVector*)muP4->At(goodLepID[0]);
    TLorentzVector* thatLep = (channel=="ele") ? (TLorentzVector*)eleP4->At(goodLepID[1]) : (TLorentzVector*)muP4->At(goodLepID[1]);

    // calculate lepton weight
    
    float thisLepWeight, thatLepWeight;

    if( channel == "ele" ){
    
      thisLepWeight = leptonWeight(h2_ele, thisLep);
      thatLepWeight = leptonWeight(h2_ele, thatLep);

    }

    else if( channel == "mu" ){

      if( isHighPtMuon[goodLepID[0]] ) 
	thisLepWeight = (thisLep->Pt() < 53) ? leptonWeight(h2_muPt20, thisLep) : leptonWeight(h2_muPt53, thisLep);
      
      else thisLepWeight = 1;
      
      if( isHighPtMuon[goodLepID[1]] )
	thatLepWeight = (thatLep->Pt() < 53) ? leptonWeight(h2_muPt20, thatLep) : leptonWeight(h2_muPt53, thatLep);
      
      else thatLepWeight = 1;
      
    }

    // select good FATjet

    int goodFATJetID = -1;

    if( !isPassJet(data, &goodFATJetID, thisLep, thatLep) ) continue;

    TLorentzVector* thisJet = (TLorentzVector*)FATjetP4->At(goodFATJetID);

    float mllbb;

    noiseCleaning(&mllbb, thisLep, thatLep, thisJet);

    // b-tag cut

    int nsubBjet = 0;

    float btagWeight = bTagWeight(data, goodFATJetID, &nsubBjet, h_l, h_c, h_b, reader_l, reader_c, reader_b);
    
    if( cat == 1 && nsubBjet != 1 ) continue;
    if( cat == 2 && nsubBjet != 2 ) continue;
        
    int iw = N-1;
    for( int nw = last; nw >= first; --nw ){
      passEvent[iw] += (pdfscaleSysWeight[nw]/pdfscaleSysWeight[first]) * eventWeight * btagWeight * thisLepWeight * thatLepWeight;
      --iw;
    }

  } // end of event loop

  // Calculate signal efficiency

  for( int nw = N-1; nw >= 0; --nw ){

    if( nw != 0 )
      efficiency[nw-1]  = passEvent[nw]/totalEvent;
    else 
      efficiencyCentral = passEvent[nw]/totalEvent;

  }

  float uncertainty = (first != 0) ? TMath::RMS(N-1, efficiency)/efficiencyCentral : TMath::Max(fabs(efficiency[1]-efficiencyCentral), fabs(efficiency[0]-efficiencyCentral))/efficiencyCentral;

  delete [] passEvent;
  delete [] efficiency;

  return ( !onlyCentral ) ? uncertainty : efficiencyCentral;

}
