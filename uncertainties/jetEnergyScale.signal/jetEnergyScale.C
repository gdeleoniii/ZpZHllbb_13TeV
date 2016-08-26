#include "jetEnergyScale.h"

void jetEnergyScale(){

  int mzh[11] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  for( int cat = 1; cat <= 2; ++cat ){

    FILE* fe = fopen(Form("ele_%ibtag_jesUnc.txt", cat), "w");
    FILE* fm = fopen(Form("mu_%ibtag_jesUnc.txt", cat), "w");

    fprintf(fe, "mass\tcentral\tJES_relativeUnc\n");
    fprintf(fm, "mass\tcentral\tJES_relativeUnc\n");

    float jes0e[11], jesUpe[11], jesDwe[11], jesUnce[11];
    float jes0m[11], jesUpm[11], jesDwm[11], jesUncm[11];

    for( int i = 0; i < 11; ++i ){

      jes0e[i]  = jetEnergyScale(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "ele", cat, mzh[i], 0);
      jesUpe[i] = fabs( jetEnergyScale(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "ele", cat, mzh[i], 1) - jes0e[i] );
      jesDwe[i] = fabs( jetEnergyScale(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "ele", cat, mzh[i], -1) - jes0e[i] );

      jesUnce[i] = TMath::Max(jesUpe[i],jesDwe[i])/jes0e[i];

      jes0m[i]  = jetEnergyScale(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "mu", cat, mzh[i], 0);
      jesUpm[i] = fabs( jetEnergyScale(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "mu", cat, mzh[i], 1) - jes0m[i] );
      jesDwm[i] = fabs( jetEnergyScale(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "mu", cat, mzh[i], -1) - jes0m[i] );

      jesUncm[i] = TMath::Max(jesUpm[i],jesDwm[i])/jes0m[i];

      fprintf(fe, "%i\t%.3f\t%.3f\n", mzh[i], jes0e[i], jesUnce[i]);      
      fprintf(fm, "%i\t%.3f\t%.3f\n", mzh[i], jes0m[i], jesUncm[i]);

    }

    fclose(fe);
    fclose(fm);

  }

  gSystem->Exec("mkdir signalJetEnScaleResults");
  gSystem->Exec("mv *txt signalJetEnScaleResults");

}
