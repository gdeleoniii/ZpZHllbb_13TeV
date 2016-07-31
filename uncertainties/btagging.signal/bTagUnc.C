#include <fstream>
#include <TSystem.h>
#include "eleBtagUnc.h"
#include "muBtagUnc.h"

void bTagUnc(){

  int mzh[11] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  for( int cat = 1; cat <= 2; ++cat ){

    FILE* fe = fopen(Form("ele_%ibtag_bTagUnc.txt", cat), "w");
    FILE* fm = fopen(Form("mu_%ibtag_bTagUnc.txt", cat), "w");

    fprintf(fe, "mass\tcentral\tbTagSFuncUp\tbTagSFuncDw\n");
    fprintf(fm, "mass\tcentral\tbTagSFuncUp\tbTagSFuncDw\n");

    float btagsf0e[11], btagsfUpe[11], btagsfDwe[11];
    float btagsf0m[11], btagsfUpm[11], btagsfDwm[11];

    for( int i = 0; i < 11; ++i ){

      btagsf0e[i]  = eleBtagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), cat, "central", mzh[i]);
      btagsfUpe[i] = eleBtagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), cat, "up", mzh[i]);
      btagsfDwe[i] = eleBtagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", mzh[i]), cat, "down", mzh[i]);

      btagsf0m[i]  = muBtagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), cat, "central", mzh[i]);
      btagsfUpm[i] = muBtagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), cat, "up", mzh[i]);
      btagsfDwm[i] = muBtagUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]), cat, "down", mzh[i]);
    
      fprintf(fe, "%i\t%g\t%g\t%g\n", (int)mzh[i], btagsf0e[i], btagsfUpe[i], btagsfDwe[i]);      
      fprintf(fm, "%i\t%g\t%g\t%g\n", (int)mzh[i], btagsf0m[i], btagsfUpm[i], btagsfDwm[i]);

    }

    fclose(fe);
    fclose(fm);

  }

  gSystem->Exec("mkdir signalbTagSFResults");
  gSystem->Exec("mv *txt signalbTagSFResults");

}
