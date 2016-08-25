#include "leptonScaleUnc.h"

void leptonScaleUnc(){

  int mzh[11] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  for( int cat = 1; cat <= 2; ++cat ){

    FILE* fe = fopen(Form("ele_%ibTag_lepScaleUnc.txt", cat), "w");
    FILE* fm = fopen(Form("mu_%ibTag_lepScaleUnc.txt", cat), "w");

    fprintf(fe, "mass\tcentral\tlepScale_relativeUnc\n");
    fprintf(fm, "mass\tcentral\tlepScale_relativeUnc\n");

    float lepScale0e[11], lepScaleUpe[11], lepScaleDwe[11], lepScaleUnce[11];
    float lepScale0m[11], lepScaleUpm[11], lepScaleDwm[11], lepScaleUncm[11];

    for( int i = 0; i < 11; ++i ){

      lepScale0e[i]  = leptonScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "ele", cat, 0, mzh[i]);
      lepScaleUpe[i] = fabs( leptonScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "ele", cat, 1, mzh[i]) - lepScale0e[i]);
      lepScaleDwe[i] = fabs( leptonScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "ele", cat, -1, mzh[i]) - lepScale0e[i]);

      lepScaleUnce[i] = TMath::Max(lepScaleUpe[i],lepScaleDwe[i])/lepScale0e[i];

      lepScale0m[i]  = leptonScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), "mu", cat, 0, mzh[i]);
      lepScaleUpm[i] = fabs( leptonScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), "mu", cat, 1, mzh[i]) - lepScale0m[i]);
      lepScaleDwm[i] = fabs( leptonScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), "mu", cat, -1, mzh[i]) - lepScale0m[i]);
    
      lepScaleUncm[i] = TMath::Max(lepScaleUpm[i],lepScaleDwm[i])/lepScale0m[i];

      fprintf(fe, "%i\t%.3f\t%.3f\n", (int)mzh[i], lepScale0e[i], lepScaleUnce[i]);
      fprintf(fm, "%i\t%.3f\t%.3f\n", (int)mzh[i], lepScale0m[i], lepScaleUncm[i]);

    }

    fclose(fe);
    fclose(fm);

  }

  gSystem->Exec("mkdir signallepSFResults");
  gSystem->Exec("mv *txt signallepSFResults");

}
