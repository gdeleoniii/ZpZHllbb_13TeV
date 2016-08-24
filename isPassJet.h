#include <vector>
#include <iostream>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/untuplizer.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagCalhead/BTagCalibrationStandalone.h"

bool isPassJet(TreeReader& data, int *goodFATJetID, TLorentzVector* thisLep=NULL, TLorentzVector* thatLep=NULL, bool jetMassCut=true, int jetScale=0){

  Int_t         FATnJet           = data.GetInt("FATnJet");    
  Float_t*      FATjetPRmassCorr  = data.GetPtrFloat("FATjetPRmassL2L3Corr");
  Float_t*      FATjetCorrUncUp   = data.GetPtrFloat("FATjetCorrUncUp");
  Float_t*      FATjetCorrUncDown = data.GetPtrFloat("FATjetCorrUncDown");
  vector<bool>& FATjetPassIDLoose = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
  TClonesArray* FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");

  bool findJet = false;  

  for( int ij = 0; ij < FATnJet; ++ij ){

    TLorentzVector* myJet = (TLorentzVector*)FATjetP4->At(ij);

    if( jetScale == 1 ) 
      *myJet *= 1+FATjetCorrUncUp[ij];
    else if( jetScale == -1 )
      *myJet *= 1-FATjetCorrUncDown[ij];
    else
      *myJet *= 1;

    if( myJet->Pt() < 200 ) continue;
    if( fabs(myJet->Eta()) > 2.4 ) continue;
    if( !FATjetPassIDLoose[ij] ) continue;
    if( myJet->DeltaR(*thisLep) < 0.8 || myJet->DeltaR(*thatLep) < 0.8 ) continue;
    if( jetMassCut && (FATjetPRmassCorr[ij] < 105 || FATjetPRmassCorr[ij] > 135) ) continue;

    *goodFATJetID = ij;
    findJet = true;

    break;
 
  } // end of FatnJet loop

  return findJet ? true : false;

}

void noiseCleaning(TreeReader& data, string channel, int thisGoodLepID, int thatGoodLepID, int goodFATJetID, float *mZH){

  TClonesArray* muP4     = (TClonesArray*) data.GetPtrTObject("muP4");
  TClonesArray* eleP4    = (TClonesArray*) data.GetPtrTObject("eleP4");
  TClonesArray* FATjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");

  TLorentzVector* thisLep = (channel=="ele") ? (TLorentzVector*)eleP4->At(thisGoodLepID) : (TLorentzVector*)muP4->At(thisGoodLepID);
  TLorentzVector* thatLep = (channel=="ele") ? (TLorentzVector*)eleP4->At(thatGoodLepID) : (TLorentzVector*)muP4->At(thatGoodLepID);
  TLorentzVector* thisJet = (TLorentzVector*)FATjetP4->At(goodFATJetID);

  if( fabs( (*thisLep+*thatLep).DeltaPhi(*thisJet) ) < 2.5 ) return;
  if( fabs( (*thisLep+*thatLep).Eta() - (*thisJet).Eta() ) > 5 ) return;
  if( (*thisLep+*thatLep+*thisJet).M() < 750 ) return;

  *mZH = (*thisLep+*thatLep+*thisJet).M();

}

float bTagWeight(TreeReader& data, int goodFATJetID, int* nsubBjet, TH1F* h_l, TH1F* h_c, TH1F* h_b,
		 BTagCalibrationReader& reader_l, BTagCalibrationReader& reader_c, BTagCalibrationReader& reader_b,
		 string region="central"){

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
      
    thisSubJet.SetPxPyPzE(FATsubjetSDPx[goodFATJetID][is], FATsubjetSDPy[goodFATJetID][is], FATsubjetSDPz[goodFATJetID][is], FATsubjetSDE[goodFATJetID][is]);

    float thisSubJetPt;

    if( thisSubJet.Pt() < 30 ) 
      thisSubJetPt = 30.001;
    else if( thisSubJet.Pt() > 2000 )
      thisSubJetPt = 1999.999;
    else
      thisSubJetPt = thisSubJet.Pt();

    float btagEff, scaleFactor;

    if( FATsubjetFlavor[goodFATJetID][is] != 4 && FATsubjetFlavor[goodFATJetID][is] != 5 ){

      btagEff = h_l->GetBinContent(h_l->FindBin(thisSubJetPt));
      scaleFactor = reader_l.eval_auto_bounds(region.data(), BTagEntry::FLAV_UDSG, thisSubJet.Eta(), thisSubJet.Pt());

    }

    else if( FATsubjetFlavor[goodFATJetID][is] == 4 ){

      btagEff = h_c->GetBinContent(h_c->FindBin(thisSubJetPt));
      scaleFactor = reader_c.eval_auto_bounds(region.data(), BTagEntry::FLAV_C, thisSubJet.Eta(), thisSubJet.Pt());

    }

    else if( FATsubjetFlavor[goodFATJetID][is] == 5 ){

      btagEff = h_b->GetBinContent(h_b->FindBin(thisSubJetPt));
      scaleFactor = reader_b.eval_auto_bounds(region.data(), BTagEntry::FLAV_B, thisSubJet.Eta(), thisSubJet.Pt());

    }

    if( FATsubjetSDCSV[goodFATJetID][is] > 0.605 ){

      pMC   *= btagEff;
      pData *= scaleFactor * btagEff;
      ++(*nsubBjet);
      
    }
      
    else{
      
      pMC   *= (1 - btagEff);
      pData *= (1 - scaleFactor * btagEff);
      
    }
      
  } // end of subjet loop

  return pMC > 0 ? pData/pMC : 0;

}
