#include <vector>
#include <string>
#include <iostream>
#include <TH1D.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "../../../../untuplizer.h"
#include "../../../../isPassZmumu.h"

void toyMCmu_cat1(std::string inputFile, std::string outputFile){

  // read the ntuples (in pcncu)

  TreeReader data(inputFile.data());

  // Declare the histogram

  TFile* f = new TFile(inputFile.data());
  TH1D* h_totalEvents = (TH1D*)f->Get("h_totalEv");

  TH1D* h_PRmassNoSIG = new TH1D("h_PRmassNoSIG", "PRmassNoSIG",    40,   40,  240);
  TH1D* h_PRmassAll   = new TH1D("h_PRmassAll",   "PRmassAll",      40,   40,  240);

  h_PRmassNoSIG->Sumw2();
  h_PRmassAll  ->Sumw2();

  h_PRmassNoSIG->GetXaxis()->SetTitle("corrPRmass");
  h_PRmassAll  ->GetXaxis()->SetTitle("corrPRmass");

  // begin of event loop
  int pass=0;
  Long64_t nentries = data.GetEntriesFast();

  for(Long64_t ev = 0; ev < nentries; ev++){

    if( ev % 1000000 == 0 )
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, nentries);

    data.GetEntry(ev);

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

      if( subjetDeltaR > 0.3 ) continue;
      if( nsubBjet < 1 ) continue;
            
      goodFATJetID = ij;
      break;

    }

    if( goodFATJetID < 0 ) continue;
    pass++;
    Float_t mllbb = (*thisLep+*thatLep+*thisJet).M();  

    if( mllbb < 750 ) continue;
    
    h_PRmassAll->Fill(corrPRmass[goodFATJetID],eventWeight);

    if( corrPRmass[goodFATJetID] > 30 && !(corrPRmass[goodFATJetID] > 65 && corrPRmass[goodFATJetID] < 135) && corrPRmass[goodFATJetID] < 300 )
      h_PRmassNoSIG->Fill(corrPRmass[goodFATJetID],eventWeight);

  } // end of event loop
  cout << pass << endl;
  TFile* outFile = new TFile(Form("%s_pseudoTest.root",outputFile.c_str()), "recreate");
  
  h_PRmassNoSIG->Write("corrPRmass");
  h_PRmassAll  ->Write("corrPRmassAll");
  h_totalEvents->Write("totalEvents");

  outFile->Write();

  delete f;
  delete outFile;

}
