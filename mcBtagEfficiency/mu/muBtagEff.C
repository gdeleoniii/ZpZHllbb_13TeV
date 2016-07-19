#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <TH1.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/untuplizer.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZmumu.h"

void muBtagEff(string inputFile, string outputFile){
  
  // read the ntuples (in pcncu)
  
  TreeReader data(inputFile.data());

  // Declare the histogram

  TFile* f = new TFile(inputFile.data());
  TH1D* h_totalEvents = (TH1D*)f->Get("h_totalEv");

  TH1F* h_lJetPtnoCSV = new TH1F("h_lJetPtnoCSV", "lJetPtnoCSV", 50, 0, 2000);
  TH1F* h_lJetPtwtCSV = new TH1F("h_lJetPtwtCSV", "lJetPtwtCSV", 50, 0, 2000);
  TH1F* h_cJetPtnoCSV = new TH1F("h_cJetPtnoCSV", "cJetPtnoCSV", 50, 0, 2000);
  TH1F* h_cJetPtwtCSV = new TH1F("h_cJetPtwtCSV", "cJetPtwtCSV", 50, 0, 2000);
  TH1F* h_bJetPtnoCSV = new TH1F("h_bJetPtnoCSV", "bJetPtnoCSV", 50, 0, 2000);
  TH1F* h_bJetPtwtCSV = new TH1F("h_bJetPtwtCSV", "bJetPtwtCSV", 50, 0, 2000);

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
    TClonesArray*  FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    vector<bool>&  FATjetPassIDLoose = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    vector<float>* FATsubjetSDCSV    = data.GetPtrVectorFloat("FATsubjetSDCSV", FATnJet);
    vector<float>* FATsubjetSDPx     = data.GetPtrVectorFloat("FATsubjetSDPx", FATnJet);
    vector<float>* FATsubjetSDPy     = data.GetPtrVectorFloat("FATsubjetSDPy", FATnJet);
    vector<float>* FATsubjetSDPz     = data.GetPtrVectorFloat("FATsubjetSDPz", FATnJet);
    vector<float>* FATsubjetSDE      = data.GetPtrVectorFloat("FATsubjetSDE", FATnJet);
    vector<int>*   FATsubjetFlavor   = data.GetPtrVectorInt("FATsubjetSDHadronFlavor", FATnJet);

    // select good reco level events     
    // select good leptons
      
    vector<int> goodLepID;

    if( !isPassZmumu(data,goodLepID) ) continue;

    TLorentzVector* thisLep = (TLorentzVector*)muP4->At(goodLepID[0]);
    TLorentzVector* thatLep = (TLorentzVector*)muP4->At(goodLepID[1]);

    // select good FATjet

    int goodFATJetID = -1;
    TLorentzVector thisJet(0,0,0,0);

    for( int ij = 0; ij < FATnJet; ++ij ){

      TLorentzVector* myJet = (TLorentzVector*)FATjetP4->At(ij);

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

    for( int is = 0; is < FATnSubSDJet[goodFATJetID]; ++is ){
      
      TLorentzVector thisSubJet;

      thisSubJet.SetPxPyPzE(FATsubjetSDPx[goodFATJetID][is],
			    FATsubjetSDPy[goodFATJetID][is],
			    FATsubjetSDPz[goodFATJetID][is],
			    FATsubjetSDE[goodFATJetID][is]);

      if( FATsubjetFlavor[goodFATJetID][is] == 1  ||
	  FATsubjetFlavor[goodFATJetID][is] == 2  ||
	  FATsubjetFlavor[goodFATJetID][is] == 3  ||
	  FATsubjetFlavor[goodFATJetID][is] == 21 ){

	h_lJetPtnoCSV->Fill(thisSubJet.Pt(),eventWeight);
	
	if( FATsubjetSDCSV[goodFATJetID][is] > 0.605 )
	  h_lJetPtwtCSV->Fill(thisSubJet.Pt(),eventWeight);

      }

      else if( FATsubjetFlavor[goodFATJetID][is] == 4 ){

	h_cJetPtnoCSV->Fill(thisSubJet.Pt(),eventWeight);

	if( FATsubjetSDCSV[goodFATJetID][is] > 0.605 )
	  h_cJetPtwtCSV->Fill(thisSubJet.Pt(),eventWeight);

      }
 
      else if( FATsubjetFlavor[goodFATJetID][is] == 5 ){

	h_bJetPtnoCSV->Fill(thisSubJet.Pt(),eventWeight);

	if( FATsubjetSDCSV[goodFATJetID][is] > 0.605 )
	  h_bJetPtwtCSV->Fill(thisSubJet.Pt(),eventWeight);

      } // end of if-else subjet flavor
 
    } // end of subjet loop                           

  } // end of event loop
  
  fprintf(stdout, "Processed all events\n");

  TFile* outFile = new TFile(Form("%s_muMCbtagEff.root",outputFile.c_str()), "recreate");

  h_lJetPtnoCSV->Write("lJetPtnoCSV");
  h_lJetPtwtCSV->Write("lJetPtwtCSV");
  h_cJetPtnoCSV->Write("cJetPtnoCSV");
  h_cJetPtwtCSV->Write("cJetPtwtCSV");
  h_bJetPtnoCSV->Write("bJetPtnoCSV");
  h_bJetPtwtCSV->Write("bJetPtwtCSV");
  h_totalEvents->Write("totalEvents");
  
  outFile->Write();

  delete f;
  delete outFile;

}
