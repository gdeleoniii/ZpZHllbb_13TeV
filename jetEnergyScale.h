#include <iostream>
#include <fstream>
#include <TH1F.h>
#include <TFile.h>
#include <TString.h>
#include <TSystem.h>
#include <TLorentzVector.h>

Float_t jetPtEtaUnc(Double_t myPt, Double_t myEta){

  TString nCol_temp = gSystem->GetFromPipe("awk 'NR==2 {print NF}' Summer15_25nsV6_MC_Uncertainty_AK8PFchs.txt");
  TString nRow_temp = gSystem->GetFromPipe("awk '{print NR}' Summer15_25nsV6_MC_Uncertainty_AK8PFchs.txt | tail -n 1");

  Int_t nCol = atoi(nCol_temp.Data());
  Int_t nRow = atoi(nRow_temp.Data());

  const Int_t nPtBins = (nCol-3)/3;
  const Int_t nEtaBins = nRow; 

  Float_t etaDw[nEtaBins];
  Float_t etaUp[nEtaBins];
  Float_t ptDw[nPtBins];
  Float_t unc[nEtaBins][nPtBins];
  Float_t dummy;

  ifstream fin;
  fin.open("Summer15_25nsV6_MC_Uncertainty_AK8PFchs.txt");

  while(!fin.eof()){

    for(Int_t ieta = 0; ieta < nEtaBins; ieta++){
      fin >> etaDw[ieta] >> etaUp[ieta] >> dummy;

      for(Int_t ipt = 0; ipt < nPtBins; ipt++){
	fin >> ptDw[ipt] >> unc[ieta][ipt] >> dummy;

      }
    }
  }

  fin.close();

  TH1F* hPt = new TH1F("hPt", "", nPtBins-1, ptDw);
  TH1F* hEta = new TH1F("hEta", "", nEtaBins-1, etaDw);

  Int_t myPtId = hPt->FindBin(myPt)-1;
  Int_t myEtaId = hEta->FindBin(myEta)-1;

  // Doing interpolation
  Float_t a = (unc[myEtaId][myPtId+1]-unc[myEtaId][myPtId])/(ptDw[myPtId+1]-ptDw[myPtId]);
  Float_t b = (unc[myEtaId][myPtId]*ptDw[myPtId+1]-unc[myEtaId][myPtId+1]*ptDw[myPtId])/(ptDw[myPtId+1]-ptDw[myPtId]);
  Float_t myUnc = a*myPt+b;

  hPt->Delete();
  hEta->Delete();

  return myUnc;

}

TLorentzVector jetEnergyScale(Int_t mode, Float_t myPt, Float_t myEta, TLorentzVector jet_original){

  TLorentzVector falseVector(0,0,0,0);

  if( mode > 1 || mode < -1 ) return falseVector;

  Float_t myUnc = jetPtEtaUnc(myPt, myEta);
  Float_t corrUnc = (1 + mode*fabs(myUnc));

  TLorentzVector jet_corr = corrUnc * jet_original;

  return jet_corr;

}
