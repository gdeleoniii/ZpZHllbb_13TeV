R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagCalhead/BTagCalibrationStandalone_cpp.so)
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/untuplizer.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZee.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZmumu.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassJet.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagCalhead/BTagCalibrationStandalone.h"

void mZHLimit(string inputFile, string outputFile, string channel, int cat, int mzh=0){

  // setup calibration and reader

  BTagCalibration calib("csvv1", "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/CSVV1.csv");

  BTagCalibrationReader reader_l(BTagEntry::OP_LOOSE, "central");
  BTagCalibrationReader reader_c(BTagEntry::OP_LOOSE, "central");
  BTagCalibrationReader reader_b(BTagEntry::OP_LOOSE, "central");

  reader_l.load(calib, BTagEntry::FLAV_UDSG, "comb");
  reader_c.load(calib, BTagEntry::FLAV_C,    "mujets");
  reader_b.load(calib, BTagEntry::FLAV_B,    "mujets");

  // to read b-tag effinciency 

  TFile* f_l = TFile::Open(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/uncertainties/btagging.zjets/bTagEffroot/%s_udsgflavor_zjetsBtagEff.root", channel.data()));
  TFile* f_c = TFile::Open(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/uncertainties/btagging.zjets/bTagEffroot/%s_cflavor_zjetsBtagEff.root", channel.data()));
  TFile* f_b = mzh==0 ? 
    TFile::Open(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/uncertainties/btagging.zjets/bTagEffroot/%s_bflavor_zjetsBtagEff.root", channel.data())) :
    TFile::Open(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/uncertainties/btagging.signal/bTagEffroot/%s_bflavor_signalBtagEff.root", channel.data()));
  
  TH1F* h_l = (TH1F*)(f_l->Get(Form("%s_udsgflavor", channel.data())));
  TH1F* h_c = (TH1F*)(f_c->Get(Form("%s_cflavor",    channel.data())));
  TH1F* h_b = mzh==0 ?
    (TH1F*)(f_b->Get(Form("%s_bflavor",     channel.data()))) :
    (TH1F*)(f_b->Get(Form("%s_bflavor_m%i", channel.data(), mzh)));

  // read the ntuples (in pcncu)
  
  TreeReader data(inputFile.data());
  TFile f(inputFile.data());

  // Declare the histogram

  TH1F* h_totalEvents = (TH1F*)f.Get("h_totalEv");
     
  TH1F* h_mZprime = new TH1F("h_mZprime", "mZprime", 100, 400, 5000);

  h_mZprime->Sumw2();
  h_mZprime->GetXaxis()->SetTitle("mZprime");

  // begin of event loop

  fprintf(stdout, "Total events %lli\n", data.GetEntriesFast());

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    if( (unsigned)ev % 100000 == 0 )
      fprintf(stdout, "Still left events %lli\n", ev);

    data.GetEntry(ev);

    Bool_t         isData            = data.GetBool("isData");
    Float_t        eventWeight       = data.GetFloat("ev_weight");
    TClonesArray*  muP4              = (TClonesArray*) data.GetPtrTObject("muP4");
    TClonesArray*  eleP4             = (TClonesArray*) data.GetPtrTObject("eleP4");
    TClonesArray*  FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");

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
    if( fabs( (*thisLep+*thatLep).DeltaPhi(*thisJet) ) < 2.5 ) continue;
    if( fabs( (*thisLep+*thatLep).Eta() - (*thisJet).Eta() ) > 5 ) continue;

    // b-tag cut

    int nsubBjet = 0;

    float btagWeight = isData ? 1 : bTagWeight(data, goodFATJetID, &nsubBjet, h_l, h_c, h_b, reader_l, reader_c, reader_b);

    if( cat == 1 && nsubBjet != 1 ) continue;
    if( cat == 2 && nsubBjet != 2 ) continue;
    
    h_mZprime->Fill((*thisLep+*thatLep+*thisJet).M(), eventWeight*btagWeight);

  } // end of event loop

  fprintf(stderr, "Processed all events\n");

  TFile* outFile = new TFile(Form("%s_mZH%sLimit.root", outputFile.data(), channel.data()), "recreate");

  h_mZprime    ->Write("mZprime");
  h_totalEvents->Write("totalEvents");

  outFile->Write();

  delete outFile;
  
}
