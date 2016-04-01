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

Double_t CrossSection(std::string token){

  std::ifstream textFile("textFile.txt");
  std::string thisSample;

  double crosssection = 0.;
  double thisNum = 0.;

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

  Int_t cat = 0;
  Float_t mllbb, prmass, evweight;

  tree->Branch("cat",      &cat,      "cat/I");
  tree->Branch("mllbb",    &mllbb,    "mllbb/F");
  tree->Branch("prmass",   &prmass,   "prmass/F");
  tree->Branch("evweight", &evweight, "evweight/F");

  // Calculate the scale correspond to inputFile

  Double_t totalEvents  = h_totalEvents->Integral();
  Double_t crossSection = CrossSection(inputFile.data());
  Double_t scale        = 3000./(totalEvents/crossSection); // dataLumi = 3000/pb

  // Mark minor backgounds

  Int_t minor = 1;

  if( (inputFile.find("WW") != std::string::npos) ||
      (inputFile.find("WZ") != std::string::npos) ||
      (inputFile.find("ZZ") != std::string::npos) ||
      (inputFile.find("TT") != std::string::npos)  ) minor = -1;

  // begin of event loop

  Long64_t nentries = data.GetEntriesFast();

  for(Long64_t ev = 0; ev < nentries; ev++){

    if( ev % 1000000 == 0 )
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, nentries);

    data.GetEntry(ev);

    Bool_t         isData            = data.GetBool("isData");
    Float_t        eventWeight       = data.GetFloat("ev_weight"); 
    TClonesArray*  muP4              = (TClonesArray*) data.GetPtrTObject("muP4");
    Int_t          FATnJet           = data.GetInt("FATnJet");    
    Int_t*         FATnSubSDJet      = data.GetPtrInt("FATnSubSDJet");
    Float_t*       corrPRmass        = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    TClonesArray*  FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    vector<bool>&  FATjetPassIDLoose = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    vector<float>* FATsubjetSDPx     = data.GetPtrVectorFloat("FATsubjetSDPx", FATnJet);
    vector<float>* FATsubjetSDPy     = data.GetPtrVectorFloat("FATsubjetSDPy", FATnJet);
    vector<float>* FATsubjetSDPz     = data.GetPtrVectorFloat("FATsubjetSDPz", FATnJet);
    vector<float>* FATsubjetSDE      = data.GetPtrVectorFloat("FATsubjetSDE", FATnJet);
    vector<float>* FATsubjetSDCSV    = data.GetPtrVectorFloat("FATsubjetSDCSV", FATnJet);

    // select good leptons
      
    vector<Int_t> goodLepID;

    TLorentzVector* thisLep = NULL;
    TLorentzVector* thatLep = NULL;

    if( !isPassZmumu(data,goodLepID) ) continue;

    thisLep = (TLorentzVector*)muP4->At(goodLepID[0]);   
    thatLep = (TLorentzVector*)muP4->At(goodLepID[1]);   

    // select good FATjet

    Int_t goodFATJetID = -1;
    TLorentzVector* thisJet = NULL;

    for(Int_t ij = 0; ij < FATnJet; ij++){

      thisJet = (TLorentzVector*)FATjetP4->At(ij);

      if( thisJet->Pt() < 200 ) continue;
      if( fabs(thisJet->Eta()) > 2.4 ) continue;
      if( !FATjetPassIDLoose[ij] ) continue;
      if( thisJet->DeltaR(*thisLep) < 0.8 || thisJet->DeltaR(*thatLep) < 0.8 ) continue;
      
      Int_t nsubBjet = 0;

      for(Int_t is = 0; is < FATnSubSDJet[ij]; is++){

	if( FATsubjetSDCSV[ij][is] > 0.605 ) nsubBjet++;

      }

      Double_t subjetDeltaR = -1;

      if( nsubBjet == 2 ){

	TLorentzVector l4_subjet0(0,0,0,0);
	TLorentzVector l4_subjet1(0,0,0,0);

	l4_subjet0.SetPxPyPzE(FATsubjetSDPx[ij][0],
			      FATsubjetSDPy[ij][0],
			      FATsubjetSDPz[ij][0],
			      FATsubjetSDE [ij][0]);

	l4_subjet1.SetPxPyPzE(FATsubjetSDPx[ij][1],
			      FATsubjetSDPy[ij][1],
			      FATsubjetSDPz[ij][1],
			      FATsubjetSDE [ij][1]);

	subjetDeltaR = l4_subjet0.DeltaR(l4_subjet1);

      }

      // deltaR depends loose cut

      if     ( subjetDeltaR < 0.3 && nsubBjet > 0 ) cat = 1;
      else if( subjetDeltaR > 0.3 && nsubBjet > 1 ) cat = 2;
      else continue;
      
      goodFATJetID = ij;
      break;

    } // end of FatnJet loop

    if( goodFATJetID < 0 ) continue;

    mllbb = (*thisLep+*thatLep+*thisJet).M();  

    if( mllbb < 750 ) continue;

    prmass = corrPRmass[goodFATJetID];

    if( isData ) 
      evweight = 1;
    else
      evweight = eventWeight * scale * minor;

    tree->Fill();

  } // end of event loop

  fprintf(stderr, "Processed all events\n");

  tree->Write();  
  outFile->Write();

  delete f;
  delete outFile;

}
