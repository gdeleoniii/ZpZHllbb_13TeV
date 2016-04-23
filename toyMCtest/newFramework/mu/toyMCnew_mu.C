#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <TH1D.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "../../../untuplizer.h"
#include "../../../isPassZmumu.h"

Float_t CrossSection(std::string token){

  std::ifstream textFile("textFile.txt");
  std::string thisSample;

  Float_t crosssection = 0.;
  Float_t thisNum = 0.;

  while( textFile >> thisSample >> thisNum ){

    if( token.find(thisSample) != std::string::npos )
      crosssection = thisNum;

  }

  return crosssection;
  
}

void toyMCnew_mu(std::string inputFile, std::string outputFile){

  // read the ntuples (in pcncu)

  TreeReader data(inputFile.data());

  TFile* f = new TFile(inputFile.data());
  TH1D* h_totalEvents = (TH1D*)f->Get("h_totalEv");

  // Create a tree to store variables

  TFile* outFile = new TFile(Form("%s_toyMCnew.root",outputFile.c_str()), "recreate");
  TTree* tree = new TTree("tree", "TreeForRooFit");

  Float_t mllbb, prmass, evweight;

  tree->Branch("mllbb",    &mllbb,    "mllbb/F");
  tree->Branch("prmass",   &prmass,   "prmass/F");
  tree->Branch("evweight", &evweight, "evweight/F");

  // Calculate the scale correspond to inputFile

  Float_t scale = 2080./(h_totalEvents->Integral()/CrossSection(inputFile.data())); // dataLumi = 2080/pb

  // Mark minor backgounds

  int minor = 1;

  if( (inputFile.find("WW") != std::string::npos) ||
      (inputFile.find("WZ") != std::string::npos) ||
      (inputFile.find("ZZ") != std::string::npos) ||
      (inputFile.find("TT") != std::string::npos)  ) minor = -1;

  // begin of event loop

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    if( (unsigned)ev % 500000 == 0 )
      fprintf(stdout, "Still left events %lli of %lli\n", ev, data.GetEntriesFast());

    data.GetEntry(ev);

    Bool_t         isData            = data.GetBool("isData");
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

    TLorentzVector *thisLep = NULL, *thatLep = NULL;

    if( !isPassZmumu(data,goodLepID) ) continue;

    thisLep = (TLorentzVector*)muP4->At(goodLepID[0]);   
    thatLep = (TLorentzVector*)muP4->At(goodLepID[1]);   

    // select good FATjet

    int goodFATJetID = -1;
    TLorentzVector* thisJet = NULL;

    for( int ij = 0; ij < FATnJet; ++ij ){

      thisJet = (TLorentzVector*)FATjetP4->At(ij);

      if( thisJet->Pt() < 200 ) continue;
      if( fabs(thisJet->Eta()) > 2.4 ) continue;
      if( !FATjetPassIDLoose[ij] ) continue;
      if( thisJet->DeltaR(*thisLep) < 0.8 || thisJet->DeltaR(*thatLep) < 0.8 ) continue;
      
      int nsubBjet = 0;

      for( int is = 0; is < FATnSubSDJet[ij]; ++is ){

	if( FATsubjetSDCSV[ij][is] > 0.605 ) ++nsubBjet;

      }

      // require at least 1 b-tag jet

      if( nsubBjet == 0 ) continue;
      
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
