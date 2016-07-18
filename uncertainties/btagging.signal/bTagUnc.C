#include <iostream>
#include <fstream>
#include <string>
#include "eleBtagUnc.h"
#include "muBtagUnc.h"

// do " .L BTagCalibrationStandalone.cpp+ " first before run this code

void bTagUnc(){

  int mzh[11] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  for( int cat = 1; cat <= 2; ++cat ){

    FILE* fe = fopen(Form("eb%d_bTagUnc.txt",cat), "w");
    FILE* fm = fopen(Form("mb%d_bTagUnc.txt",cat), "w");

    // mass, mur=1, pdf

    for( int i = 0; i < 11; ++i ){

      fprintf(fe, "%d\t%g\t%g\t%g\n",
	      mzh[i],
	      eleBtagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), cat, 0),
	      eleBtagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), cat, 1),
	      eleBtagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), cat, 2));

      fprintf(fm, "%d\t%g\t%g\t%g\n",
	      mzh[i],
	      muBtagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), cat, 0),
	      muBtagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), cat, 1),
	      muBtagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), cat, 2));
    
    }

    fclose(fe);
    fclose(fm);

  }

}
