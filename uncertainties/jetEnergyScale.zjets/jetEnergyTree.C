R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagCalhead/BTagCalibrationStandalone_cpp.so)
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/untuplizer.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZee.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZmumu.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassJet.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagCalhead/BTagCalibrationStandalone.h"

float crossSection(string thisPath){

  ifstream textFile("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/xSec.txt");
  string token;
  float crosssection = 0., thisNum = 0.;

  while( textFile >> token >> thisNum ){

    if( thisPath.find(token) != string::npos )
      crosssection = thisNum;

  }

  return crosssection;

}

void jetEnergyTree(string inputFile, string outputFile, string jes, string channel){

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
  TFile* f_b = TFile::Open(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagEffroot/%s_bflavor_zjetsBtagEff.root",    channel.data()));
  
  TH1F* h_l = (TH1F*)(f_l->Get(Form("%s_udsgflavor", channel.data())));
  TH1F* h_c = (TH1F*)(f_c->Get(Form("%s_cflavor",    channel.data())));
  TH1F* h_b = (TH1F*)(f_b->Get(Form("%s_bflavor",    channel.data())));

  // read the ntuples (in pcncu)

  TreeReader data(inputFile.data());

  TFile* f = new TFile(inputFile.data());
  TH1D* h_totalEvents = (TH1D*)f->Get("h_totalEv");

  // Create a tree to store variables

  TFile* outFile = new TFile(Form("%s_%s_%sMiniTree.root", outputFile.data(), jes.data(), channel.data()), "recreate");
  TTree* tree = new TTree("tree", "TreeForRooFit");

  Int_t   cat;
  Float_t prmass, evweight, mllbb;

  tree->Branch("cat",      &cat,      "cat/I");
  tree->Branch("mllbb",    &mllbb,    "mllbb/F");
  tree->Branch("prmass",   &prmass,   "prmass/F");
  tree->Branch("evweight", &evweight, "evweight/F");

  // Calculate the scale correspond to inputFile

  Float_t scale = 2512.*crossSection(outputFile.data())/h_totalEvents->Integral();

  // begin of event loop

  fprintf(stdout, "Total events %lli\n", data.GetEntriesFast());

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    if( (unsigned)ev % 100000 == 0 )
      fprintf(stdout, "Still left events %lli\n", ev);

    data.GetEntry(ev);

    Float_t       eventWeight = data.GetFloat("ev_weight");
    TClonesArray* muP4        = (TClonesArray*) data.GetPtrTObject("muP4");
    TClonesArray* eleP4       = (TClonesArray*) data.GetPtrTObject("eleP4");
    TClonesArray* FATjetP4    = (TClonesArray*) data.GetPtrTObject("FATjetP4");

    // select good leptons
      
    vector<int> goodLepID;

    if( channel == "ele" && !isPassZee(data,goodLepID)   ) continue;
    if( channel == "mu"  && !isPassZmumu(data,goodLepID) ) continue;

    TLorentzVector* thisLep = (channel=="ele") ? (TLorentzVector*)eleP4->At(goodLepID[0]) : (TLorentzVector*)muP4->At(goodLepID[0]);
    TLorentzVector* thatLep = (channel=="ele") ? (TLorentzVector*)eleP4->At(goodLepID[1]) : (TLorentzVector*)muP4->At(goodLepID[1]);

    // select good FATjet

    int goodFATJetID = -1;

    int JES = 0;

    if     ( jes == "up"   ) JES =  1;
    else if( jes == "down" ) JES = -1;

    if( !isPassJet(data, &goodFATJetID, thisLep, thatLep, false, JES) ) continue;

    TLorentzVector* thisJet = (TLorentzVector*)FATjetP4->At(goodFATJetID);

    if( (*thisLep+*thatLep+*thisJet).M() < 750 ) continue;

    // b-tag cut

    int nsubBjet = 0;

    float btagWeight = bTagWeight(data, goodFATJetID, &nsubBjet, h_l, h_c, h_b, reader_l, reader_c, reader_b);

    if     ( nsubBjet == 1 ) cat = 1;
    else if( nsubBjet == 2 ) cat = 2;      
    else                     cat = 0;
    
    mllbb    = (*thisLep+*thatLep+thisJet).M(); 
    prmass   = FATjetPRmassCorr[goodFATJetID];
    evweight = eventWeight * scale * btagWeight;

    tree->Fill();

  } // end of event loop

  fprintf(stdout, "Processed all events\n");

  tree->Write();  
  outFile->Write();

  delete f;
  delete outFile;

}
