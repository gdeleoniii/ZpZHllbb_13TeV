#include <iostream>
#include <cstdio>
#include "elePileUpWeight.h"
#include "muPileUpWeight.h"

void pileUpWeight(){

  int mzh[11] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  for( int cat = 1; cat <= 2; ++cat ){

    FILE* fe = fopen(Form("eb%d_pileUpUnc.txt",cat), "w");
    FILE* fm = fopen(Form("mb%d_pileUpUnc.txt",cat), "w");

    // mass, central, up, down

    for( int i = 0; i < 11; ++i ){
    
      fprintf(fe, "%d\t%g\t%g\t%g\n",
	      mzh[i],
	      elePileUpWeight(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",     mzh[i]), cat),
	      elePileUpWeight(Form("/data7/htong/skim_signalPileUpScaleUp/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), cat),
	      elePileUpWeight(Form("/data7/htong/skim_signalPileUpScaleDw/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), cat));

      fprintf(fm, "%d\t%g\t%g\t%g\n",
	      mzh[i],
	      muPileUpWeight(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",     mzh[i]), cat),
	      muPileUpWeight(Form("/data7/htong/skim_signalPileUpScaleUp/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), cat),
	      muPileUpWeight(Form("/data7/htong/skim_signalPileUpScaleDw/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), cat));
    
    }

    fclose(fe);
    fclose(fm);

  }

}
