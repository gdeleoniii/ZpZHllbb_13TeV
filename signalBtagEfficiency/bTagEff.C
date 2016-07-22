#include <iostream>
#include <TFile.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include "eleBtagEff.h"
#include "muBtagEff.h"

void bTagEff(string channel, string flavor){

  int mzh[11] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  TFile f(Form("%s_%sflavor_signalBtagEff.root", channel.data(),flavor.data()),"recreate");

  for( int i = 0; i < 11; ++i ){

    TGraphAsymmErrors *g = new TGraphAsymmErrors();

    g = (channel=="ele") ? 
      eleBtagEff(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",mzh[i]), flavor.data()) :
      muBtagEff (Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",mzh[i]),  flavor.data());

    g->GetXaxis()->SetLimits(0,2000);
    g->Write(Form("%s_%sflavor_m%i",channel.data(),flavor.data(),mzh[i]));
    f.Write();

    TLegend leg(0.60, 0.70, 0.90, 0.87);

    leg.SetBorderSize(0);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.04);
    leg.AddEntry(g, Form("%s, %s flavor",channel.data(),flavor.data()), "lp");

    TLatex lar;

    lar.SetNDC(kTRUE);
    lar.SetTextSize(0.04);
    lar.SetLineWidth(5);

    TCanvas c("c", "", 0, 0, 800, 600);

    c.cd();
    g->Draw("ap");
    leg.Draw();
    lar.DrawLatex(0.15, 0.83, "CMS");
    lar.DrawLatex(0.15, 0.79, "#it{#bf{Simulation}}");
    lar.DrawLatex(0.62, 0.83, Form("#bf{m_{ZH} = %i GeV}",mzh[i]));
    c.Print(Form( (i==0) ? "%s_%sflavor_signalBtagEff.pdf(" : ( (i==10) ? "%s_%sflavor_signalBtagEff.pdf)" : "%s_%sflavor_signalBtagEff.pdf" ), channel.data(),flavor.data()));

  } // end of mass point

}
