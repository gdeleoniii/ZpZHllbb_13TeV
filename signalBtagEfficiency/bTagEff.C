#include <iostream>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include "eleBtagEff.h"
#include "muBtagEff.h"

const int N = 11;

void bTagEff(string channel, string Cat){

  int cat = std::stoi(Cat);
  int mzh[N] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  TCanvas c("c", "", 0, 0, 800, 600);

  for( int i = 0; i < N; ++i ){

    TGraphAsymmErrors *g = new TGraphAsymmErrors();

    if( channel == "ele" )
      g = eleBtagEff(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root",mzh[i]),cat);
    else
      g = muBtagEff(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root",mzh[i]),cat);

    TLegend leg(0.60, 0.70, 0.90, 0.87);

    leg.SetBorderSize(0);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.04);
    leg.AddEntry(g, Form("ele %i b-tag", cat), "lp");

    TLatex lar;

    lar.SetNDC(kTRUE);
    lar.SetTextSize(0.04);
    lar.SetLineWidth(5);

    c.cd();
    g->Draw("ap");
    leg.Draw();
    lar.DrawLatex(0.15, 0.83, "CMS");
    lar.DrawLatex(0.15, 0.79, "#it{#bf{Simulation}}");
    lar.DrawLatex(0.62, 0.83, Form("#bf{m_{ZH} = %i GeV}", mzh[i]));
    c.Print(Form((i==0) ? "%sBtagEff_%ibtag.pdf(" : ( (i==N-1) ? "%sBtagEff_%ibtag.pdf)" : "%sBtagEff_%ibtag.pdf" ), channel.data(), cat));

  }

}
