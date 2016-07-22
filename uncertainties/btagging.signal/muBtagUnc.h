R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/uncertainties/btagging.signal/bTagCalhead/BTagCalibrationStandalone_cpp.so)

#include <vector>
#include <string>
#include <iostream>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/untuplizer.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZmumu.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/uncertainties/btagging.signal/bTagCalhead/BTagCalibrationStandalone.h"

float muBtagUnc(string inputFile, int cat, int jSF,
		TGraphAsymmErrors* g_l, TGraphAsymmErrors* g_c, TGraphAsymmErrors* g_b){

  // setup calibration and reader

  string region = (jSF==0) ? "central" : ((jSF==1) ? "up" : "down");

  BTagCalibration       calib ("csvv1", "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/uncertainties/btagging.signal/CSVV1.csv");
  BTagCalibrationReader reader_udsg(BTagEntry::OP_LOOSE, region.data());
  BTagCalibrationReader reader_c(BTagEntry::OP_LOOSE, region.data());
  BTagCalibrationReader reader_b(BTagEntry::OP_LOOSE, region.data());

  reader_udsg.load(calib, BTagEntry::FLAV_UDSG, "mujets");
  reader_c.load(calib, BTagEntry::FLAV_C, "mujets");
  reader_b.load(calib, BTagEntry::FLAV_B, "mujets");

  // read the ntuples (in pcncu)
  
  TreeReader data(inputFile.data());
  
  float passEvent = 0.;

  // begin of event loop

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    data.GetEntry(ev);

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

    int nsubBjet = 0;
    float pMC = 1., pData = 1.;

    for( int is = 0; is < FATnSubSDJet[goodFATJetID]; ++is ){
            
      TLorentzVector thisSubJet;
      
      thisSubJet.SetPxPyPzE(FATsubjetSDPx[goodFATJetID][is],
			    FATsubjetSDPy[goodFATJetID][is],
			    FATsubjetSDPz[goodFATJetID][is],
			    FATsubjetSDE[goodFATJetID][is]);

      float btagEff, scaleFactor;

      if( FATsubjetFlavor[goodFATJetID][is] != 4 && FATsubjetFlavor[goodFATJetID][is] != 5 ){

	btagEff = g_l->Eval(thisSubJet.Pt());
	scaleFactor = reader_udsg.eval_auto_bounds(region.data(), BTagEntry::FLAV_UDSG, thisSubJet.Eta(), thisSubJet.Pt());

      }

      else if( FATsubjetFlavor[goodFATJetID][is] == 4 ){

	btagEff = g_c->Eval(thisSubJet.Pt());
	scaleFactor = reader_c.eval_auto_bounds(region.data(), BTagEntry::FLAV_C, thisSubJet.Eta(), thisSubJet.Pt());

      }

      else if( FATsubjetFlavor[goodFATJetID][is] == 5 ){

	btagEff = g_b->Eval(thisSubJet.Pt());
	scaleFactor = reader_b.eval_auto_bounds(region.data(), BTagEntry::FLAV_B, thisSubJet.Eta(), thisSubJet.Pt());

      }

      if( FATsubjetSDCSV[goodFATJetID][is] > 0.605 ){

	pMC   *= btagEff;
	pData *= scaleFactor * btagEff;
	++nsubBjet;
	
      }
      
      else{
	
	pMC   *= (1 - btagEff);
	pData *= (1 - scaleFactor * btagEff);
	
      }
      
    } // end of subjet for loop
    
    // b-tag cut
    
    if( cat == 1 && nsubBjet != 1 ) continue;
    if( cat == 2 && nsubBjet != 2 ) continue;
        
    passEvent += (pData/pMC);

  } // end of event loop
  
  return passEvent/(float)data.GetEntriesFast();

}
