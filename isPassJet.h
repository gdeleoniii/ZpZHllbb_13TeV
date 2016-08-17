#include <vector>
#include <iostream>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/untuplizer.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagCalhead/BTagCalibrationStandalone.h"

bool isPassJet(TreeReader& data, int goodFATJetID, bool jetMassCut=true, 
	       TLorentzVector* thisLep=NULL, TLorentzVector* thatLep=NULL){

  Int_t          FATnJet           = data.GetInt("FATnJet");    
  Float_t*       FATjetPRmassCorr  = data.GetPtrFloat("FATjetPRmassL2L3Corr");
  TClonesArray*  FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");
  vector<bool>&  FATjetPassIDLoose = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));

  bool findJet = false;  

  TLorentzVector thisJet(0,0,0,0);

  for( int ij = 0; ij < FATnJet; ++ij ){

    TLorentzVector* myJet = (TLorentzVector*)FATjetP4->At(ij);

    if( myJet->Pt() < 200 ) continue;
    if( fabs(myJet->Eta()) > 2.4 ) continue;
    if( !FATjetPassIDLoose[ij] ) continue;
    if( myJet->DeltaR(*thisLep) < 0.8 || myJet->DeltaR(*thatLep) < 0.8 ) continue;
    if( jetMassCut && (FATjetPRmassCorr[ij] < 105 || FATjetPRmassCorr[ij] > 135) ) continue;
      
    goodFATJetID = ij;
    thisJet = *myJet;
    findJet = true;

    break;
 
  } // end of FatnJet loop

  return findJet ? true : false;

}

float bTagWeight(TreeReader& data, int goodFATJetID, int nsubBjet, 
		 BTagCalibrationReader& reader_l, BTagCalibrationReader& reader_c, BTagCalibrationReader& reader_b,
		 TGraphAsymmErrors* g_l, TGraphAsymmErrors* g_c, TGraphAsymmErrors* g_b){

  Int_t          FATnJet         = data.GetInt("FATnJet");
  Int_t*         FATnSubSDJet    = data.GetPtrInt("FATnSubSDJet");
  vector<float>* FATsubjetSDCSV  = data.GetPtrVectorFloat("FATsubjetSDCSV", FATnJet);
  vector<float>* FATsubjetSDPx   = data.GetPtrVectorFloat("FATsubjetSDPx", FATnJet);
  vector<float>* FATsubjetSDPy   = data.GetPtrVectorFloat("FATsubjetSDPy", FATnJet);
  vector<float>* FATsubjetSDPz   = data.GetPtrVectorFloat("FATsubjetSDPz", FATnJet);
  vector<float>* FATsubjetSDE    = data.GetPtrVectorFloat("FATsubjetSDE", FATnJet);
  vector<int>*   FATsubjetFlavor = data.GetPtrVectorInt("FATsubjetSDHadronFlavor", FATnJet);
 
  // Calculate b-tag weights

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
      scaleFactor = reader_l.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, thisSubJet.Eta(), thisSubJet.Pt());

    }

    else if( FATsubjetFlavor[goodFATJetID][is] == 4 ){

      btagEff = g_c->Eval(thisSubJet.Pt());
      scaleFactor = reader_c.eval_auto_bounds("central", BTagEntry::FLAV_C, thisSubJet.Eta(), thisSubJet.Pt());

    }
    else if( FATsubjetFlavor[goodFATJetID][is] == 5 ){

      btagEff = g_b->Eval(thisSubJet.Pt());
      scaleFactor = reader_b.eval_auto_bounds("central", BTagEntry::FLAV_B, thisSubJet.Eta(), thisSubJet.Pt());

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
      
  } // end of subjet loop

  return pData/pMC;

}
