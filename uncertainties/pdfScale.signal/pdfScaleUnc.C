#include "pdfScaleUnc.h"

void pdfScaleUnc(){

  int mzh[11] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  for( int cat = 1; cat <= 2; ++cat ){

    FILE* fpe = fopen(Form("ele_%ibtag_pdfUnc.txt", cat), "w");
    FILE* fpm = fopen(Form("mu_%ibtag_pdfUnc.txt", cat), "w");

    FILE* fse = fopen(Form("ele_%ibtag_scaleUnc.txt", cat), "w");
    FILE* fsm = fopen(Form("mu_%ibtag_scaleUnc.txt", cat), "w");

    fprintf(fpe, "mass\tcentral\tPDF\n");
    fprintf(fpm, "mass\tcentral\tPDF\n");

    fprintf(fse, "mass\tcentral\tQCD\n");
    fprintf(fsm, "mass\tcentral\tQCD\n");

    float centraleS[11], pdfUnce[11], scaleUnce[11];
    float centralmS[11], pdfUncm[11], scaleUncm[11];

    for( int i = 0; i < 11; ++i ){

      centraleS[i] = pdfScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "ele", cat, mzh[i], 0, 2, true);
      scaleUnce[i] = pdfScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "ele", cat, mzh[i], 0, 2);
      pdfUnce[i]   = pdfScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "ele", cat, mzh[i], 9, 109);

      centralmS[i] = pdfScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "mu", cat, mzh[i], 0, 2, true);
      scaleUncm[i] = pdfScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "mu", cat, mzh[i], 0, 2);
      pdfUncm[i]   = pdfScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "mu", cat, mzh[i], 9, 109);

      fprintf(fpe, "%i\t%.3f\t%.3f\n", mzh[i], centraleS[i], pdfUnce[i]+1);      
      fprintf(fpm, "%i\t%.3f\t%.3f\n", mzh[i], centralmS[i], pdfUncm[i]+1);
    
      fprintf(fse, "%i\t%.3f\t%.3f\n", mzh[i], centraleS[i], scaleUnce[i]+1);
      fprintf(fsm, "%i\t%.3f\t%.3f\n", mzh[i], centralmS[i], scaleUncm[i]+1);

    }

    fclose(fpe);
    fclose(fpm);

    fclose(fse);
    fclose(fsm);

  }

  gSystem->Exec("mkdir signalPdfResults; mkdir signalScaleResults");
  gSystem->Exec("mv *pdfUnc.txt signalPdfResults");
  gSystem->Exec("mv *scaleUnc.txt signalScaleResults");

}
