#include "jetEnergyScale.h"

void jetEnergyScale(){

  int mzh[11] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  for( int cat = 1; cat <= 2; ++cat ){

    FILE* fe = fopen(Form("ele_%ibtag_jesUnc.txt", cat), "w");
    FILE* fm = fopen(Form("mu_%ibtag_jesUnc.txt", cat), "w");

    fprintf(fe, "mass\tcentral\tJESuncUp\tJESuncDw\n");
    fprintf(fm, "mass\tcentral\tJESuncUp\tJESuncDw\n");

    float jes0e[11], jesUpe[11], jesDwe[11];
    float jes0m[11], jesUpm[11], jesDwm[11];

    for( int i = 0; i < 11; ++i ){

      jes0e[i]  = jetEnergyScale(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "central", "ele", cat, mzh[i]);
      jesUpe[i] = fabs( jetEnergyScale(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "up", "ele", cat, mzh[i]) - jes0e[i] );
      jesDwe[i] = fabs( jetEnergyScale(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "dw", "ele", cat, mzh[i]) - jes0e[i] );

      jes0m[i]  = jetEnergyScale(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "central", "mu", cat, mzh[i]);
      jesUpm[i] = fabs( jetEnergyScale(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "up", "mu", cat, mzh[i]) - jes0m[i] );
      jesDwm[i] = fabs( jetEnergyScale(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "dw", "mu", cat, mzh[i]) - jes0m[i] );

      fprintf(fe, "%i\t%g\t%g\t%g\n", mzh[i], jes0e[i], jesUpe[i], jesDwe[i]);      
      fprintf(fm, "%i\t%g\t%g\t%g\n", mzh[i], jes0m[i], jesUpm[i], jesDwm[i]);

    }

    fclose(fe);
    fclose(fm);

  }

  gSystem->Exec("mkdir signalJetEnScaleResults");
  gSystem->Exec("mv *txt signalJetEnScaleResults");

}
