#include "pdfScaleUnc.h"

void pdfScaleUnc(){

  int mzh[11] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  for( int cat = 1; cat <= 2; ++cat ){

    FILE* fe = fopen(Form("ele_%ibtag_pdfScaleUnc.txt", cat), "w");
    FILE* fm = fopen(Form("mu_%ibtag_pdfScaleUnc.txt", cat), "w");

    fprintf(fe, "mass\tcentral\tpdfUnc\tcentral\tscaleUnc\n");
    fprintf(fm, "mass\tcentral\tpdfUnc\tcentral\tscaleUnc\n");

    float centraleP[11], centraleS[11], pdfUnce[11], scaleUnce[11];
    float centralmP[11], centralmS[11], pdfUncm[11], scaleUncm[11];

    for( int i = 0; i < 11; ++i ){

      centraleS[i] = pdfScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "ele", cat, mzh[i], 0, 2, true);
      scaleUnce[i] = pdfScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "ele", cat, mzh[i], 0, 2);
      centraleP[i] = pdfScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "ele", cat, mzh[i], 9, 109, true);
      pdfUnce[i]   = pdfScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "ele", cat, mzh[i], 9, 109);

      centralmS[i] = pdfScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "mu", cat, mzh[i], 0, 2, true);
      scaleUncm[i] = pdfScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "mu", cat, mzh[i], 0, 2);
      centralmP[i] = pdfScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "mu", cat, mzh[i], 9, 109, true);
      pdfUncm[i]   = pdfScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "mu", cat, mzh[i], 9, 109);

      fprintf(fe, "%i\t%.3f\t%.3f\t%.3f\t%.3f\n", mzh[i], centraleP[i], pdfUnce[i], centraleS[i], scaleUnce[i]);      
      fprintf(fm, "%i\t%.3f\t%.3f\t%.3f\t%.3f\n", mzh[i], centralmP[i], pdfUncm[i], centralmS[i], scaleUncm[i]);
    
    }

    fclose(fe);
    fclose(fm);

  }

  gSystem->Exec("mkdir signalPdfScaleResults");
  gSystem->Exec("mv *txt signalPdfScaleResults");

}
