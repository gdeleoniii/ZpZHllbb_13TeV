#include <vector>
#include <string>
#include <iostream>
#include <TH1D.h>
#include <TGraph.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "../untuplizer.h"
#include "../isPassZmumu.h"

float getEfficiency(string inputFile){

  // read the ntuples (in pcncu)

  TreeReader data(inputFile.data());
  
  int genGoodMu  = 0;
  int recoGoodMu = 0;

  // begin of event loop

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    if( (unsigned)ev % 500000 == 0 )
      fprintf(stdout, "Still left events %lli of %lli\n", ev, data.GetEntriesFast());

    data.GetEntry(ev);

    Int_t    nGenPar     = data.GetInt("nGenPar");
    Int_t*   genParId    = data.GetPtrInt("genParId");
    Int_t*   genParSt    = data.GetPtrInt("genParSt");
    Int_t*   genMomParId = data.GetPtrInt("genMomParId");
    Int_t    nMu         = data.GetInt("nMu");
    Float_t* muMiniIsoEA = data.GetPtrFloat("muMiniIsoEA");
    TClonesArray* muP4   = (TClonesArray*) data.GetPtrTObject("muP4");
    vector<bool>& isHighPtMuon        = *((vector<bool>*) data.GetPtr("isHighPtMuon"));
    vector<bool>& isCustomTrackerMuon = *((vector<bool>*) data.GetPtr("isCustomTrackerMuon"));

    // select good generator level muons

    for(Int_t ig = nGenPar-1; ig >= 0; --ig){

      if( abs(genParId[ig]) != 13 || genParSt[ig] != 1 ) continue;
      if( genMomParId[ig] != 23 && genMomParId[ig] != genParId[ig] ) continue;

      ++genGoodMu;

    }

    // select good reco level muons

    for(Int_t im = nMu-1; im >= 0; --im){
      
      TLorentzVector* myMu = (TLorentzVector*)muP4->At(im);

      if( muMiniIsoEA[im] > 0.2 ) continue;
      if( !isHighPtMuon[im] && !isCustomTrackerMuon[im] ) continue;
      if( fabs(myMu->Eta()) > 2.4 ) continue;
      if( myMu->Pt() < 20 ) continue;

      ++recoGoodMu;

    }

  } // end of event loop

  fprintf(stdout, "Processed all events\n");
  
  return (float)recoGoodMu/(float)genGoodMu;

}

void signalEfficiency(){

  Float_t x_mzh[13] = {600,800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500};
  Float_t y_eff[13];

  for( int i = 0; i < 13; ++i )
    y_eff[i] = getEfficiency(Form("/data7/htong/skim_samples/mu/skim_mu_ZprimeToZhToZlephbb_narrow_M-%f_13TeV-madgraph.root",x_mzh[i]));

  TGraph g_eff(13, x_mzh, y_eff);

  TCanvas c("c", "", 0, 0, 800, 600);

  c.cd();
  g_eff.Draw();
  c.Print("muIDeff.pdf");

}
