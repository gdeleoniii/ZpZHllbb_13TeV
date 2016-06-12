#include <iostream>
#include <fstream>
#include <string>
#include "elePdfScaleUnc.h"
#include "muPdfScaleUnc.h"

void pdfScaleUnc(int cat, int first, int last, int iter, string type){

  int mzh[11] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  fstream fe,fm;

  fe.open(Form("eb%d_%sScale.txt",cat,type.c_str()), ios::out);
  fm.open(Form("mb%d_%sScale.txt",cat,type.c_str()), ios::out);

  for( int i = 0; i < 11; ++i ){
    
    fe << mzh[i] << "\t" << elePdfScaleUnc(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",mzh[i]),cat,first,last,iter) << endl;
    fm << mzh[i] << "\t" << muPdfScaleUnc (Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", mzh[i]),cat,first,last,iter) << endl;

  }

  fe.close();
  fm.close();

}
