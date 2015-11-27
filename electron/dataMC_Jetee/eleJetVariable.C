#include <vector>
#include <string>
#include <iostream>
#include <TH1D.h>
#include <TMath.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "../../untuplizer.h"
#include "../../readSample.h"
#include "../../isPassZee.h"

void eleJetVariable(std::string inputFile, std::string outputFile){

  // read the ntuples (in pcncu)

  std::vector<string> infiles;
    
  readSample(inputFile, infiles);
  
  TreeReader data(infiles);
  
  // Declare the histogram

  Int_t nBin = 20;
     
  TH1D* h_FATjetPt         = new TH1D("h_FATjetPt", "FATjetPt", nBin, 100, 1000);
  TH1D* h_FATjetEta        = new TH1D("h_FATjetEta", "FATjetEta", nBin, -4, 4);
  TH1D* h_FATjetCISVV2     = new TH1D("h_FATjetCISVV2", "FATjetCISVV2", nBin, 0, 1.2);
  TH1D* h_FATjetSDmass     = new TH1D("h_FATjetSDmass", "FATjetSDmass", nBin, 0, 200);
  TH1D* h_FATjetPRmass     = new TH1D("h_FATjetPRmass", "FATjetPRmass", nBin, 0, 200);
  TH1D* h_FATjetTau1       = new TH1D("h_FATjetTau1", "FATjetTau1", nBin, 0, 1);
  TH1D* h_FATjetTau2       = new TH1D("h_FATjetTau2", "FATjetTau2", nBin, 0, 1);
  TH1D* h_FATjetTau2dvTau1 = new TH1D("h_FATjetTau2dvTau1", "FATjetTau2dvTau1", nBin, 0, 1);
  TH1D* h_FATsubjetPt      = new TH1D("h_FATsubjetPt", "FATsubjetPt", nBin, 0, 800);
  TH1D* h_FATsubjetEta     = new TH1D("h_FATsubjetEta", "FATsubjetEta", nBin, -4, 4);  
  TH1D* h_FATsubjetSDCSV   = new TH1D("h_FATsubjetSDCSV", "FATsubjetSDCSV", nBin, 0, 1.2);
  TH1D* h_eventWeight      = new TH1D("h_eventWeight", "eventWeight", 100, -1, 1);

  h_FATjetPt        ->Sumw2();   
  h_FATjetEta       ->Sumw2();
  h_FATjetCISVV2    ->Sumw2();
  h_FATjetSDmass    ->Sumw2();
  h_FATjetPRmass    ->Sumw2();
  h_FATjetTau1      ->Sumw2();
  h_FATjetTau2      ->Sumw2();
  h_FATjetTau2dvTau1->Sumw2();
  h_FATsubjetPt     ->Sumw2();
  h_FATsubjetEta    ->Sumw2();
  h_FATsubjetSDCSV  ->Sumw2();

  h_FATjetPt        ->GetXaxis()->SetTitle("FATjetPt");
  h_FATjetEta       ->GetXaxis()->SetTitle("FATjetEta");
  h_FATjetCISVV2    ->GetXaxis()->SetTitle("FATjetCISVV2");
  h_FATjetSDmass    ->GetXaxis()->SetTitle("FATjetSDmass");
  h_FATjetPRmass    ->GetXaxis()->SetTitle("FATjetPRmass");
  h_FATjetTau1      ->GetXaxis()->SetTitle("FATjetTau1");
  h_FATjetTau2      ->GetXaxis()->SetTitle("FATjetTau2");
  h_FATjetTau2dvTau1->GetXaxis()->SetTitle("FATjetTau2dvTau1");
  h_FATsubjetPt     ->GetXaxis()->SetTitle("FATsubjetPt");
  h_FATsubjetEta    ->GetXaxis()->SetTitle("FATsubjetEta");
  h_FATsubjetSDCSV  ->GetXaxis()->SetTitle("FATsubjetSDCSV");
  h_eventWeight     ->GetXaxis()->SetTitle("eventWeight");
    
  // begin of event loop

  Int_t nPass[4] = {0};

  for( Long64_t ev = 0; ev < data.GetEntriesFast(); ev++ ){

    if( ev % 1000000 == 0 )
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());

    data.GetEntry(ev);

    Int_t          nVtx              = data.GetInt("nVtx");
    Float_t        mcWeight          = data.GetFloat("mcWeight");
    TClonesArray*  eleP4             = (TClonesArray*) data.GetPtrTObject("eleP4");
    Int_t          FATnJet           = data.GetInt("FATnJet");    
    Int_t*         FATnSubSDJet      = data.GetPtrInt("FATnSubSDJet");
    Float_t*       FATjetCISVV2      = data.GetPtrFloat("FATjetCISVV2");
    Float_t*       FATjetSDmass      = data.GetPtrFloat("FATjetSDmass");
    Float_t*       FATjetPRmass      = data.GetPtrFloat("FATjetPRmass");
    Float_t*       FATjetTau1        = data.GetPtrFloat("FATjetTau1");
    Float_t*       FATjetTau2        = data.GetPtrFloat("FATjetTau2");
    TClonesArray*  FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    vector<bool>&  FATjetPassIDLoose = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    vector<float>* FATsubjetSDPx     = data.GetPtrVectorFloat("FATsubjetSDPx", FATnJet);
    vector<float>* FATsubjetSDPy     = data.GetPtrVectorFloat("FATsubjetSDPy", FATnJet);
    vector<float>* FATsubjetSDPz     = data.GetPtrVectorFloat("FATsubjetSDPz", FATnJet);
    vector<float>* FATsubjetSDE      = data.GetPtrVectorFloat("FATsubjetSDE", FATnJet);
    vector<float>* FATsubjetSDCSV    = data.GetPtrVectorFloat("FATsubjetSDCSV", FATnJet);

    Double_t eventWeight = mcWeight;
    if( inputFile.find("DYJets") != std::string::npos ){
      if( eventWeight > 0 ) eventWeight = 1;
      else if( eventWeight < 0 ) eventWeight = -1;
    }
    else
      eventWeight = 1;
    
    h_eventWeight->Fill(0.,eventWeight);
    
    if( nVtx < 1 ) continue;

    // data trigger cut (electron channel)

    std::string* trigName = data.GetPtrString("hlt_trigName");
    vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
    const Int_t nsize = data.GetPtrStringSize();    
    bool passTrigger = false;
    
    for(Int_t it = 0; it < nsize; it++){
    
      std::string thisTrig = trigName[it];
      bool results = trigResult[it];
      
      if( thisTrig.find("HLT_Ele105") != std::string::npos && results==1 ){
	passTrigger = true;
	break;
      }
      
    }

    if( !passTrigger ) continue;

    nPass[0]++;

    // select good Electrons

    vector<Int_t> goodEleID;
    if( !isPassZee(data, goodEleID) ) continue;

    nPass[1]++;

    TLorentzVector* thisEle = (TLorentzVector*)eleP4->At(goodEleID[0]);
    TLorentzVector* thatEle = (TLorentzVector*)eleP4->At(goodEleID[1]);

    // select good FATjet

    Int_t goodFATJetID = -1;
    TLorentzVector* thisJet = NULL;

    for(Int_t ij = 0; ij < FATnJet; ij++){

      thisJet = (TLorentzVector*)FATjetP4->At(ij);

      if( thisJet->Pt() < 200 ) continue;
      if( fabs(thisJet->Eta()) > 2.4 ) continue;
      if( !FATjetPassIDLoose[ij] ) continue;
      if( FATnSubSDJet[ij] != 2 ) continue;
      if( thisJet->DeltaR(*thisEle) < 0.8 || thisJet->DeltaR(*thatEle) < 0.8 ) continue;
      
      goodFATJetID = ij;
      break;

    }

    if( goodFATJetID < 0 ) continue;
    nPass[2]++;

    h_FATjetPt        ->Fill(thisJet->Pt(),eventWeight);
    h_FATjetEta       ->Fill(thisJet->Eta(),eventWeight);
    h_FATjetCISVV2    ->Fill(FATjetCISVV2[goodFATJetID],eventWeight);
    h_FATjetSDmass    ->Fill(FATjetSDmass[goodFATJetID],eventWeight);
    h_FATjetPRmass    ->Fill(FATjetPRmass[goodFATJetID],eventWeight);
    h_FATjetTau1      ->Fill(FATjetTau1[goodFATJetID],eventWeight);
    h_FATjetTau2      ->Fill(FATjetTau2[goodFATJetID],eventWeight);
    h_FATjetTau2dvTau1->Fill(FATjetTau2[goodFATJetID]/FATjetTau1[goodFATJetID],eventWeight);

    TLorentzVector l4_subjet0(0,0,0,0);
    TLorentzVector l4_subjet1(0,0,0,0);

    l4_subjet0.SetPxPyPzE(FATsubjetSDPx[goodFATJetID][0],
			  FATsubjetSDPy[goodFATJetID][0],
			  FATsubjetSDPz[goodFATJetID][0],
			  FATsubjetSDE[goodFATJetID][0]);

    l4_subjet1.SetPxPyPzE(FATsubjetSDPx[goodFATJetID][1],
                          FATsubjetSDPy[goodFATJetID][1],
                          FATsubjetSDPz[goodFATJetID][1],
                          FATsubjetSDE[goodFATJetID][1]);

    h_FATsubjetPt   ->Fill(l4_subjet0.Pt(),eventWeight);
    h_FATsubjetPt   ->Fill(l4_subjet1.Pt(),eventWeight);
    h_FATsubjetEta  ->Fill(l4_subjet0.Eta(),eventWeight);
    h_FATsubjetEta  ->Fill(l4_subjet1.Eta(),eventWeight);
    h_FATsubjetSDCSV->Fill(FATsubjetSDCSV[goodFATJetID][0],eventWeight);
    h_FATsubjetSDCSV->Fill(FATsubjetSDCSV[goodFATJetID][1],eventWeight);

  } // end of event loop

  fprintf(stderr, "Processed all events\n");

  std::cout << "\nnPass[0] = " << nPass[0] 
	    << "\nnPass[1] = " << nPass[1] 
	    << "\nnPass[2] = " << nPass[2] 
	    << std::endl;

  std::string h_name[12] = {"FATjetPt","FATjetEta","FATjetCISVV2","FATjetSDmass",
			    "FATjetPRmass","FATjetTau1","FATjetTau2","FATjetTau2dvTau1",
			    "FATsubjetPt","FATsubjetEta","FATsubjetSDCSV","eventWeight"};

  TFile* outFile = new TFile(Form("%s_jeteeVariable.root",outputFile.c_str()), "recreate");
      
  h_FATjetPt        ->Write(h_name[0].data());
  h_FATjetEta       ->Write(h_name[1].data());
  h_FATjetCISVV2    ->Write(h_name[2].data());
  h_FATjetSDmass    ->Write(h_name[3].data());
  h_FATjetPRmass    ->Write(h_name[4].data());
  h_FATjetTau1      ->Write(h_name[5].data());
  h_FATjetTau2      ->Write(h_name[6].data());
  h_FATjetTau2dvTau1->Write(h_name[7].data());
  h_FATsubjetPt     ->Write(h_name[8].data());
  h_FATsubjetEta    ->Write(h_name[9].data());
  h_FATsubjetSDCSV  ->Write(h_name[10].data());
  h_eventWeight     ->Write(h_name[11].data());
  
  outFile->Write();
  
}
