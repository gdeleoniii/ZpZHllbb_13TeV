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

  // read the ntuples (in pcncu)

  TreeReader data(inputFile.data());

  TFile* f = new TFile(inputFile.data());
  TH1D* h_totalEvents = (TH1D*)f->Get("h_totalEv");

  // Create a tree to store variables

  TFile* outFile = new TFile(Form("%s_toyMC.root",outputFile.c_str()), "recreate");
  TTree* tree = new TTree("tree", "TreeForRooFit");

  Int_t   cat;
  Float_t prmass, evweight;
  Float_t mllbb[3] = {0};

  for(int i = 0; i < 3; ++i)
    tree->Branch(Form("mllbb%i",i),  &(mllbb[i]),  Form("mllbb%i/F",i));

  tree->Branch("cat", &cat, "cat/I");
  tree->Branch("prmass", &prmass, "prmass/F");
  tree->Branch("evweight", &evweight, "evweight/F");

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
    Float_t*       FATjetPRmassCorr  = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    Float_t*       FATjetCorrUncUp   = data.GetPtrFloat("FATjetCorrUncUp");
    Float_t*       FATjetCorrUncDown = data.GetPtrFloat("FATjetCorrUncDown");
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
    TLorentzVector thisJet(0,0,0,0);
    
    for( int js = 0; js < 3; ++js ){

      for( Int_t ij = 0; ij < FATnJet; ++ij ){

	TLorentzVector* myJet = (TLorentzVector*)FATjetP4->At(ij);

	*myJet *= (js==0) ? 1 : ( (js==1) ? (1+FATjetCorrUncUp[ij]) : (1-FATjetCorrUncDown[ij]) );

	if( myJet->Pt() < 200 ) continue;
	if( fabs(myJet->Eta()) > 2.4 ) continue;
	if( !FATjetPassIDLoose[ij] ) continue;
	if( myJet->DeltaR(*thisLep) < 0.8 || myJet->DeltaR(*thatLep) < 0.8 ) continue;
      
	Int_t nsubBjet = 0;

	for( Int_t is = 0; is < FATnSubSDJet[ij]; ++is ){

	  if( FATsubjetSDCSV[ij][is] > 0.605 ) nsubBjet++;

	}

	// b-tag cut

	if     ( nsubBjet == 1 ) cat = 1;
	else if( nsubBjet == 2 ) cat = 2;      
	else                     cat = 0;

	if( js == 0 )
	  goodFATJetID = ij;

	thisJet = *myJet;

	break;

      } // end of FatnJet loop

      mllbb[js] = (*thisLep+*thatLep+thisJet).M();
      
    }

    prmass = FATjetPRmassCorr[goodFATJetID];
    evweight = eventWeight * scale;

    if( goodFATJetID < 0 ) continue;
    if( mllbb[0] < 750 ) continue;

    tree->Fill();

  } // end of event loop

  fprintf(stdout, "Processed all events\n");

  tree->Write();  
  outFile->Write();

  delete f;
  delete outFile;

}
