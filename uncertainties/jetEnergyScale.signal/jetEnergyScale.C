#include <iostream>
#include <cstdio>
#include "eleJetEnergyScale.h"
#include "muJetEnergyScale.h"

void pileUpWeight(){

  int mzh[11] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  for( int cat = 1; cat <= 2; ++cat ){

    FILE* fe = fopen(Form("eb%d_jetEnergyUnc.txt",cat), "w");
    FILE* fm = fopen(Form("mb%d_jetEnergyUnc.txt",cat), "w");

    // mass, central, up, down

    for( int i = 0; i < 11; ++i ){
    
      fprintf(fe, "%d\t%g\t%g\t%g\n",
	      mzh[i],
	      eleJetEnergyScale(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), cat, 0),
	      eleJetEnergyScale(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), cat, 1),
	      eleJetEnergyScale(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), cat, -1));

      fprintf(fm, "%d\t%g\t%g\t%g\n",
	      mzh[i],
	      muJetEnergyScale(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), cat, 0),
	      muJetEnergyScale(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), cat, 1),
	      muJetEnergyScale(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), cat, -1));
    
    }

    fclose(fe);
    fclose(fm);

  }

}
