#include "bTagUnc.h"

void bTagUnc(){

  int mzh[11] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  for( int cat = 1; cat <= 2; ++cat ){

    FILE* fe = fopen(Form("ele_%ibtag_bTagUnc.txt", cat), "w");
    FILE* fm = fopen(Form("mu_%ibtag_bTagUnc.txt", cat), "w");

    fprintf(fe, "mass\tcentral\tbTagSF_relativeUnc\n");
    fprintf(fm, "mass\tcentral\tbTagSF_relativeUnc\n");

    float btagsf0e[11], btagsfUpe[11], btagsfDwe[11], btagsfUnce[11];
    float btagsf0m[11], btagsfUpm[11], btagsfDwm[11], btagsfUncm[11];

    for( int i = 0; i < 11; ++i ){

      btagsf0e[i]  = btagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "ele", cat, "central", mzh[i]);
      btagsfUpe[i] = fabs( btagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "ele", cat, "up", mzh[i]) - btagsf0e[i]);
      btagsfDwe[i] = fabs( btagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), "ele", cat, "down", mzh[i]) - btagsf0e[i]);

      btagsfUnce[i] = TMath::Max(btagsfUpe[i],btagsfDwe[i])/btagsf0e[i];

      btagsf0m[i]  = btagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), "mu", cat, "central", mzh[i]);
      btagsfUpm[i] = fabs( btagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), "mu", cat, "up", mzh[i]) - btagsf0m[i]);
      btagsfDwm[i] = fabs( btagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), "mu", cat, "down", mzh[i]) - btagsf0m[i]);
    
      btagsfUncm[i] = TMath::Max(btagsfUpm[i],btagsfDwm[i])/btagsf0m[i];

      fprintf(fe, "%i\t%.3f\t%.3f\n", (int)mzh[i], btagsf0e[i], btagsfUnce[i]);
      fprintf(fm, "%i\t%.3f\t%.3f\n", (int)mzh[i], btagsf0m[i], btagsfUncm[i]);

    }

    fclose(fe);
    fclose(fm);

  }

  gSystem->Exec("mkdir signalbTagSFResults");
  gSystem->Exec("mv *txt signalbTagSFResults");

}
