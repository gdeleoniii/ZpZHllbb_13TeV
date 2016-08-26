#include <iostream>
#include <TFile.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include "bTaggingEff.h"

void bTaggingEff(string channel, string flavor="b"){

  float varBins[] = {30,50,70,100,140,200,300,670,2000};
  int   nvarBins  = sizeof(varBins)/sizeof(varBins[0])-1;

  int mzh[11] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  TGraphAsymmErrors *g[11];

  for( int i = 0; i < 11; ++i ){

    g[i] = btaggingEff(Form("/data7/htong/skim_NCUGlobalTuples/skim_%s_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root", channel.data(), mzh[i]), channel.data());
    
    TLegend leg(0.60, 0.70, 0.90, 0.87);

    leg.SetBorderSize(0);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.04);
    leg.AddEntry(g[i], Form("%s, %s flavor",channel.data(),flavor.data()), "lp");

    TLatex lar;

    lar.SetNDC(kTRUE);
    lar.SetTextSize(0.04);
    lar.SetLineWidth(5);

    TCanvas c("c", "", 0, 0, 800, 600);

    c.cd();
    g[i]->GetXaxis()->SetLimits(0,2000);
    g[i]->Draw("ap");
    leg.Draw();
    lar.DrawLatex(0.15, 0.83, "CMS");
    lar.DrawLatex(0.15, 0.79, "#it{#bf{Simulation}}");
    lar.DrawLatex(0.62, 0.83, Form("#bf{m_{ZH} = %i GeV}",mzh[i]));
    c.Print(Form( (i==0) ? "%s_%sflavor_signalBtagEff.pdf(" : ( (i==10) ? "%s_%sflavor_signalBtagEff.pdf)" : "%s_%sflavor_signalBtagEff.pdf" ), channel.data(),flavor.data()));

  } // end of mass point

  // Convert the TGraph to TH1F and write to a root file

  TFile f(Form("%s_%sflavor_signalBtagEff.root", channel.data(),flavor.data()),"recreate");
  
  for( int i = 0; i < 11; ++i ){
   
    TH1F h("h", "", nvarBins, varBins);

    for( int n = 0; n < g[i]->GetN(); ++n ){
      
      double x, y;
      g[i]->GetPoint(n,x,y);
      h.SetBinContent(h.FindBin(x), y);
      
    }
    
    h.Write(Form("%s_%sflavor_m%i",channel.data(),flavor.data(),mzh[i]));
    
  } // end of mass point 
  
}
