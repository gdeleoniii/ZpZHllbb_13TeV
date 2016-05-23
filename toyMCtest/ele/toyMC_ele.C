#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <TH1D.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "../../untuplizer.h"
#include "../../isPassZee.h"

Float_t CrossSection(string token){

  ifstream textFile("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/xSec.txt");
  string thisSample;

  Float_t crosssection = 0., thisNum = 0.;

  while( textFile >> thisSample >> thisNum ){

    if( token.find(thisSample) != string::npos )
      crosssection = thisNum;

  }

  return crosssection;
  
}

void toyMC_ele(string inputFile, string outputFile){

  // read the ntuples (in pcncu)

  TreeReader data(inputFile.data());

  TFile* f = new TFile(inputFile.data());
  TH1D* h_totalEvents = (TH1D*)f->Get("h_totalEv");

  // Create a tree to store variables

  TFile* outFile = new TFile(Form("%s_toyMC.root",outputFile.c_str()), "recreate");
  TTree* tree = new TTree("tree", "TreeForRooFit");

  Int_t   cat;
  Float_t mllbb, prmass, evweight;

  tree->Branch("cat",      &cat,      "cat/I");
  tree->Branch("mllbb",    &mllbb,    "mllbb/F");
  tree->Branch("prmass",   &prmass,   "prmass/F");
  tree->Branch("evweight", &evweight, "evweight/F");

  // Calculate the scale correspond to inputFile

  Float_t scale = 2512./(h_totalEvents->Integral()/CrossSection(inputFile.data()));

  // Mark minor backgounds

  Int_t minor = 1;

  if( (inputFile.find("WW") != string::npos) ||
      (inputFile.find("WZ") != string::npos) ||
      (inputFile.find("ZZ") != string::npos) ||
      (inputFile.find("ZH") != string::npos) ||
      (inputFile.find("TT") != string::npos)  ) minor = -1;

  // begin of event loop

  fprintf(stdout, "Total events %lli\n", data.GetEntriesFast());

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    if( (unsigned)ev % 100000 == 0 )
      fprintf(stdout, "Still left events %lli\n", ev);

    data.GetEntry(ev);

    Bool_t         isData            = data.GetBool("isData");
    Float_t        eventWeight       = data.GetFloat("ev_weight"); 
    TClonesArray*  eleP4             = (TClonesArray*) data.GetPtrTObject("eleP4");
    Int_t          FATnJet           = data.GetInt("FATnJet");    
    Int_t*         FATnSubSDJet      = data.GetPtrInt("FATnSubSDJet");
    Float_t*       corrPRmass        = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    TClonesArray*  FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    vector<bool>&  FATjetPassIDLoose = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    vector<float>* FATsubjetSDCSV    = data.GetPtrVectorFloat("FATsubjetSDCSV", FATnJet);

    // select good leptons
      
    vector<int> goodLepID;
    if( !isPassZee(data,goodLepID) ) continue;

    TLorentzVector* thisLep = (TLorentzVector*)eleP4->At(goodLepID[0]);   
    TLorentzVector* thatLep = (TLorentzVector*)eleP4->At(goodLepID[1]);   

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

    evweight = isData ? 1 : eventWeight * scale * minor;

    tree->Fill();

  } // end of event loop

  fprintf(stdout, "Processed all events\n");

  tree->Write();  
  outFile->Write();

  delete f;
  delete outFile;

}
