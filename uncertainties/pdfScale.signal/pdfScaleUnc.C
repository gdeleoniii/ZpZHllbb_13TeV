#include <iostream>
#include <fstream>
#include <string>
#include "elePdfScaleUnc.h"
#include "muPdfScaleUnc.h"

void pdfScaleUnc(){

  int mzh[11] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  for( int cat = 1; cat <= 2; ++cat ){

    FILE* fe = fopen(Form("eb%d_pdfScaleUnc.txt",cat), "w");
    FILE* fm = fopen(Form("mb%d_pdfScaleUnc.txt",cat), "w");

    // mass, mur=1, pdf

    for( int i = 0; i < 11; ++i ){

      fprintf(fe, "%d\t%g\t%g\n",
	      mzh[i],
	      elePdfScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), cat, 0,2,1),
	      elePdfScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), cat, 10,109,1));

      fprintf(fm, "%d\t%g\t%g\n",
	      mzh[i],
	      muPdfScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), cat, 0,2,1),
	      muPdfScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), cat, 10,109,1));
    
    }

    fclose(fe);
    fclose(fm);

  }

}
