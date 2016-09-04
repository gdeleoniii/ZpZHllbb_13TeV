#include "leptonTriggerUnc.h"

void leptonTriggerUnc(){

  int mzh[11] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  for( int cat = 1; cat <= 2; ++cat ){

    FILE* fe = fopen(Form("ele_%ibTag_lepTriggerUnc.txt", cat), "w");
    FILE* fm = fopen(Form("mu_%ibtag_lepTriggerUnc.txt", cat), "w");

    fprintf(fe, "mass\tcentral\tTrigger\n");
    fprintf(fm, "mass\tcentral\tTrigger\n");

    float lepTrigger0m[11], lepTriggerUpm[11], lepTriggerDwm[11], lepTriggerUncm[11];
    float lepTrigger0e[11], finalUncm[11];

    for( int i = 0; i < 11; ++i ){

      lepTrigger0e[i]  = leptonTriggerUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), "ele", cat, 0, mzh[i]);

      lepTrigger0m[i]  = leptonTriggerUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), "mu", cat, 0, mzh[i]);
      lepTriggerUpm[i] = fabs( leptonTriggerUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), "mu", cat, 1, mzh[i]) - lepTrigger0m[i]);
      lepTriggerDwm[i] = fabs( leptonTriggerUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), "mu", cat, -1, mzh[i]) - lepTrigger0m[i]);
    
      lepTriggerUncm[i] = 1+TMath::Max(lepTriggerUpm[i],lepTriggerDwm[i])/lepTrigger0m[i];

      // additional uncertainties (0.5% for trigger), only for muon channel in 2015 analysis

      finalUncm[i] = TMath::Sqrt(lepTriggerUncm[i]*lepTriggerUncm[i] + 0.005*0.005);

      fprintf(fe, "%i\t%.3f\t-\n",    (int)mzh[i], lepTrigger0e[i]);
      fprintf(fm, "%i\t%.3f\t%.3f\n", (int)mzh[i], lepTrigger0m[i], finalUncm[i]);

    }

    fclose(fe);
    fclose(fm);

  }

  gSystem->Exec("mkdir signalTriggerResults");
  gSystem->Exec("mv *txt signalTriggerResults");

}
