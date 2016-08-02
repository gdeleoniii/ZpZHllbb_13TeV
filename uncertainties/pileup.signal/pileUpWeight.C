#include "pileUpWeight.h"

void pileUpWeight(){

  int mzh[11] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  for( int cat = 1; cat <= 2; ++cat ){

    FILE* fe = fopen(Form("ele_%ibtag_pileUpUnc.txt", cat), "w");
    FILE* fm = fopen(Form("mu_%ibtag_pileUpUnc.txt", cat), "w");

    fprintf(fe, "mass\tcentral\tPUuncUp\tPUuncDw\n");
    fprintf(fm, "mass\tcentral\tPUuncUp\tPUuncDw\n");

    float pu0e[11], puUpe[11], puDwe[11];
    float pu0m[11], puUpm[11], puDwm[11];

    for( int i = 0; i < 11; ++i ){
    
      pu0e[i]  = pileUpWeight(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "ele", cat, mzh[i]);
      puUpe[i] = fabs( pileUpWeight(Form("/data7/htong/skim_signalPileUpScaleUp/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "ele", cat, mzh[i]) - pu0e[i] );
      puDwe[i] = fabs( pileUpWeight(Form("/data7/htong/skim_signalPileUpScaleDw/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "ele", cat, mzh[i]) - pu0e[i] );

      pu0m[i]  = pileUpWeight(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "mu", cat, mzh[i]);
      puUpm[i] = fabs( pileUpWeight(Form("/data7/htong/skim_signalPileUpScaleUp/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "mu", cat, mzh[i]) - pu0m[i] );
      puDwm[i] = fabs( pileUpWeight(Form("/data7/htong/skim_signalPileUpScaleDw/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "mu", cat, mzh[i]) - pu0m[i] );

      fprintf(fe, "%i\t%g\t%g\t%g\n", mzh[i], pu0e[i], puUpe[i], puDwe[i]);      
      fprintf(fm, "%i\t%g\t%g\t%g\n", mzh[i], pu0m[i], puUpm[i], puDwm[i]);
    
    }

    fclose(fe);
    fclose(fm);

  }

  gSystem->Exec("mkdir signalpuWeightResults");
  gSystem->Exec("mv *txt signalpuWeightResults");

}
