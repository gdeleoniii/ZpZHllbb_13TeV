R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/uncertainties/btagging.zjets/BTagCalibrationStandalone_cpp.so)
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <TH1D.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/untuplizer.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZmumu.h"
#include "../BTagCalibrationStandalone.h"

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

void toyMC_mu(string inputFile, string outputFile){

  // setup calibration and reader

  string region[3] = {"central", "up", "down"};

  BTagCalibration calib("csvv1", "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/uncertainties/btagging.zjets/CSVV1.csv");
  BTagCalibrationReader* reader[3];

  // read the ntuples (in pcncu)

  TreeReader data(inputFile.data());

  TFile* f = new TFile(inputFile.data());
  TH1D* h_totalEvents = (TH1D*)f->Get("h_totalEv");

  // Create a tree to store variables

  TFile* outFile = new TFile(Form("%s_toyMC.root",outputFile.c_str()), "recreate");
  TTree* tree = new TTree("tree", "TreeForRooFit");

  Int_t   cat;
  Float_t mllbb, prmass;
  Float_t evweight[3] = {0};

  tree->Branch("cat",      &cat,     "cat/I");
  tree->Branch("mllbb",    &mllbb,   "mllbb/F");
  tree->Branch("prmass",   &prmass,  "prmass/F");

  for(int i = 0; i < 3; ++i){

    tree->Branch(Form("evweight%i",i), &(evweight[i]), Form("evweight%i/F",i));

    reader[i] = new BTagCalibrationReader(BTagEntry::OP_LOOSE, region[i].data());
    reader[i]->load(calib, BTagEntry::FLAV_B, "mujets");

  }

  // Calculate the scale correspond to inputFile

  Float_t scale = 2512.*crossSection(outputFile.data())/h_totalEvents->Integral();

  // begin of event loop

  fprintf(stdout, "Total events %lli\n", data.GetEntriesFast());

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    if( (unsigned)ev % 100000 == 0 )
      fprintf(stdout, "Still left events %lli\n", ev);

    data.GetEntry(ev);

    Float_t        eventWeight       = data.GetFloat("ev_weight"); 
    TClonesArray*  muP4              = (TClonesArray*) data.GetPtrTObject("muP4");
    Int_t          FATnJet           = data.GetInt("FATnJet");    
    Int_t*         FATnSubSDJet      = data.GetPtrInt("FATnSubSDJet");
    Float_t*       corrPRmass        = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    TClonesArray*  FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    vector<bool>&  FATjetPassIDLoose = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    vector<float>* FATsubjetSDCSV    = data.GetPtrVectorFloat("FATsubjetSDCSV", FATnJet);

    // select good leptons
      
    vector<int> goodLepID;
    if( !isPassZmumu(data,goodLepID) ) continue;

    TLorentzVector* thisLep = (TLorentzVector*)muP4->At(goodLepID[0]);   
    TLorentzVector* thatLep = (TLorentzVector*)muP4->At(goodLepID[1]);   

    // select good FATjet

    Int_t goodFATJetID = -1;
    TLorentzVector* thisJet = NULL;

    for( Int_t ij = 0; ij < FATnJet; ++ij ){

      thisJet = (TLorentzVector*)FATjetP4->At(ij);

      if( thisJet->Pt() < 200 ) continue;
      if( fabs(thisJet->Eta()) > 2.4 ) continue;
      if( !FATjetPassIDLoose[ij] ) continue;
      if( thisJet->DeltaR(*thisLep) < 0.8 || thisJet->DeltaR(*thatLep) < 0.8 ) continue;
      
      Int_t nsubBjet = 0;

      for( Int_t is = 0; is < FATnSubSDJet[ij]; ++is ){

	if( FATsubjetSDCSV[ij][is] > 0.605 ) nsubBjet++;

      }

      // b-tag cut

      if     ( nsubBjet == 1 ) cat = 1;
      else if( nsubBjet == 2 ) cat = 2;      
      else                     cat = 0;

      goodFATJetID = ij;
      break;

    } // end of FatnJet loop

    if( goodFATJetID < 0 ) continue;

    mllbb = (*thisLep+*thatLep+*thisJet).M();  

    if( mllbb < 750 ) continue;

    prmass = corrPRmass[goodFATJetID];

    for( int i = 0; i < 3; ++i ){
      evweight[i] = eventWeight * scale * reader[i]->eval_auto_bounds(region[i].data(), BTagEntry::FLAV_B, thisJet->Eta(), thisJet->Pt());
    }

    tree->Fill();

  } // end of event loop

  fprintf(stdout, "Processed all events\n");

  tree->Write();  
  outFile->Write();

  delete f;
  delete outFile;

}
