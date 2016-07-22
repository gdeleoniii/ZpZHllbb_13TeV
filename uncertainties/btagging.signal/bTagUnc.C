#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include "eleBtagUnc.h"
#include "muBtagUnc.h"

void bTagUnc(){

  int mzh[11] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  TFile* fe_l = TFile::Open("bTagEffroot/ele_udsgflavor_zjetsBtagEff.root");
  TFile* fe_c = TFile::Open("bTagEffroot/ele_cflavor_zjetsBtagEff.root");
  TFile* fe_b = TFile::Open("bTagEffroot/ele_bflavor_signalBtagEff.root");

  TFile* fm_l = TFile::Open("bTagEffroot/mu_udsgflavor_zjetsBtagEff.root");
  TFile* fm_c = TFile::Open("bTagEffroot/mu_cflavor_zjetsBtagEff.root");
  TFile* fm_b = TFile::Open("bTagEffroot/mu_bflavor_signalBtagEff.root");

  for( int cat = 1; cat <= 2; ++cat ){

    FILE* fe = fopen(Form("eb%d_bTagUnc.txt",cat), "w");
    FILE* fm = fopen(Form("mb%d_bTagUnc.txt",cat), "w");

    for( int m = 0; m < 11; ++m ){

      TGraphAsymmErrors* ge_l = (TGraphAsymmErrors*)(fe_l->Get("ele_udsgflavor"));
      TGraphAsymmErrors* ge_c = (TGraphAsymmErrors*)(fe_c->Get("ele_cflavor"));
      TGraphAsymmErrors* ge_b = (TGraphAsymmErrors*)(fe_b->Get(Form("ele_bflavor_m%i",mzh[m])));

      TGraphAsymmErrors* gm_l = (TGraphAsymmErrors*)(fm_l->Get("mu_udsgflavor"));
      TGraphAsymmErrors* gm_c = (TGraphAsymmErrors*)(fm_c->Get("mu_cflavor"));
      TGraphAsymmErrors* gm_b = (TGraphAsymmErrors*)(fm_b->Get(Form("mu_bflavor_m%i",mzh[m])));

      fprintf(fe, "%d\t%g\t%g\t%g\n",
	      mzh[m],
	      eleBtagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[m]), cat, 0, ge_l, ge_c, ge_b),
	      eleBtagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[m]), cat, 1, ge_l, ge_c, ge_b),
	      eleBtagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[m]), cat, 2, ge_l, ge_c, ge_b));

      fprintf(fm, "%d\t%g\t%g\t%g\n",
	      mzh[m],
	      muBtagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[m]), cat, 0, gm_l, gm_c, gm_b),
	      muBtagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[m]), cat, 1, gm_l, gm_c, gm_b),
	      muBtagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[m]), cat, 2, gm_l, gm_c, gm_b));
    
    }

    fclose(fe);
    fclose(fm);

  }

}
