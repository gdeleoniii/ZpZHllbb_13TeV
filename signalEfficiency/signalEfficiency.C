#include "signalEfficiency.h"

void signalEfficiency(){

  float mzh[11] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  float y_eff_elc1[11], y_eff_elc2[11];
  float y_eff_muc1[11], y_eff_muc2[11];

  for( int i = 0; i < 11; ++i ){

    y_eff_elc1[i] = signalEfficiency(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", (int)mzh[i]), "ele", 1, (int)mzh[i]);
    y_eff_elc2[i] = signalEfficiency(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", (int)mzh[i]), "ele", 2, (int)mzh[i]);
    y_eff_muc1[i] = signalEfficiency(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", (int)mzh[i]), "mu", 1, (int)mzh[i]);
    y_eff_muc2[i] = signalEfficiency(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%i_13TeV-madgraph.root", (int)mzh[i]), "mu", 2, (int)mzh[i]);

  }

  TGraph *g_eff_elc1 = new TGraph(11, mzh, y_eff_elc1);
  TGraph *g_eff_elc2 = new TGraph(11, mzh, y_eff_elc2);
  TGraph *g_eff_muc1 = new TGraph(11, mzh, y_eff_muc1);
  TGraph *g_eff_muc2 = new TGraph(11, mzh, y_eff_muc2);

  g_eff_elc1->GetXaxis()->SetLimits(600,4200);
  g_eff_elc1->SetMinimum(0);
  g_eff_elc1->SetMaximum(0.4);
  g_eff_elc1->SetTitle("");
  g_eff_elc1->GetXaxis()->SetTitle("m_{ZH} (GeV)");
  g_eff_elc1->GetYaxis()->SetTitle("Acceptance #times efficiency");  
  g_eff_elc1->GetYaxis()->SetTitleOffset(1.3);
  g_eff_elc1->SetLineStyle(3);
  g_eff_elc1->SetLineWidth(2);
  g_eff_elc1->SetLineColor(kGreen+1);
  g_eff_elc1->SetMarkerColor(kGreen+1);
  g_eff_elc1->SetMarkerStyle(8);

  g_eff_elc2->SetLineWidth(2);
  g_eff_elc2->SetLineColor(kGreen+1);
  g_eff_elc2->SetMarkerColor(kGreen+1);
  g_eff_elc2->SetMarkerStyle(8);

  g_eff_muc1->SetLineStyle(3);
  g_eff_muc1->SetLineWidth(2);
  g_eff_muc1->SetLineColor(kBlue);
  g_eff_muc1->SetMarkerColor(kBlue);
  g_eff_muc1->SetMarkerStyle(8);

  g_eff_muc2->SetLineWidth(2);
  g_eff_muc2->SetLineColor(kBlue);
  g_eff_muc2->SetMarkerColor(kBlue);
  g_eff_muc2->SetMarkerStyle(8);
  
  TLegend leg(0.65, 0.70, 0.90, 0.87);

  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.04);
  leg.AddEntry(g_eff_elc1, "2e, 1 b-tag", "lp");
  leg.AddEntry(g_eff_elc2, "2e, 2 b-tag", "lp");
  leg.AddEntry(g_eff_muc1, "2#mu, 1 b-tag", "lp");
  leg.AddEntry(g_eff_muc2, "2#mu, 2 b-tag", "lp");

  TLatex lar;

  lar.SetTextSize(0.03);
  lar.SetLineWidth(5);

  TCanvas c("c", "", 0, 0, 800, 600);

  c.cd();
  g_eff_elc1->Draw();
  g_eff_elc1->Draw("3lpsame");
  g_eff_elc2->Draw("3lpsame");
  g_eff_muc1->Draw("3lpsame");
  g_eff_muc2->Draw("3lpsame");
  leg.Draw();
  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.60, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  c.Print("signalEfficiency.pdf");

  gSystem->Exec("mv signalEfficiency.pdf $HOME/www");

}
