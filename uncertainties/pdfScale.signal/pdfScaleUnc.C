#include <fstream>
#include <TSystem.h>
#include "elePdfScaleUnc.h"
#include "muPdfScaleUnc.h"

void pdfScaleUnc(){

  int mzh[11] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  for( int cat = 1; cat <= 2; ++cat ){

    FILE* fe = fopen(Form("ele_%ibtag_pdfScaleUnc.txt", cat), "w");
    FILE* fm = fopen(Form("mu_%ibtag_pdfScaleUnc.txt", cat), "w");

    fprintf(fe, "mass\tpdfUnc\tscaleUnc\n");
    fprintf(fm, "mass\tpdfUnc\tscaleUnc\n");

    float pdfUnce[11], scaleUnce[11];
    float pdfUncm[11], scaleUncm[11];

    for( int i = 0; i < 11; ++i ){

      scaleUnce[i] = elePdfScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), cat, 0, 2);
      pdfUnce[i]   = elePdfScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), cat, 10, 109);

      scaleUncm[i] = muPdfScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), cat, 0, 2);
      pdfUncm[i]   = muPdfScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), cat, 10, 109);

      fprintf(fe, "%i\t%g\t%g\n", mzh[i], pdfUnce[i], scaleUnce[i]);      
      fprintf(fm, "%i\t%g\t%g\n", mzh[i], pdfUncm[i], scaleUncm[i]);
    
    }

    fclose(fe);
    fclose(fm);

  }

  gSystem->Exec("mkdir signalPdfScaleResults");
  gSystem->Exec("mv *txt signalPdfScaleResults");

}
