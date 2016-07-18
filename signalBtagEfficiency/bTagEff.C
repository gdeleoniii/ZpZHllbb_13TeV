#include <iostream>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include "eleBtagEff.h"
#include "muBtagEff.h"

void bTagEff(string channel, string flavor){

  TGraphAsymmErrors *g = new TGraphAsymmErrors();

  g = (channel=="ele") ? eleBtagEff(flavor.data()) : muBtagEff(flavor.data());
  
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
  lar.DrawLatex(0.62, 0.83, "#bf{All m_{ZH} mass points}");
  c.Print(Form("%s_%sflavor_signalBtagEff.pdf", channel.data(), flavor.data()));

}
