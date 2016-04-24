#include <iostream>
#include <TAxis.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include "getEleEfficiency.h"
#include "getMuEfficiency.h"

const int N = 11;

void signalEfficiency(){

  Float_t x_mzh[N] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  Float_t y_eff_elc1[N], y_eyh_elc1[N], y_eyl_elc1[N];
  Float_t y_eff_elc2[N], y_eyh_elc2[N], y_eyl_elc2[N];

  Float_t y_eff_muc1[N], y_eyh_muc1[N], y_eyl_muc1[N];
  Float_t y_eff_muc2[N], y_eyh_muc2[N], y_eyl_muc2[N];

  for( int i = 0; i < N; ++i ){

    y_eff_elc1[i] = getEleEfficiency(Form("/data7/htong/skim_samples/ele/skim_ele_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",(Int_t)x_mzh[i]),1,0);
    cout << "aaaaa" << endl;    
y_eyh_elc1[i] = getEleEfficiency(Form("/data7/htong/skim_samples/ele/skim_ele_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",(Int_t)x_mzh[i]),1,1);
    y_eyl_elc1[i] = getEleEfficiency(Form("/data7/htong/skim_samples/ele/skim_ele_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",(Int_t)x_mzh[i]),1,-1);

    y_eff_elc2[i] = getEleEfficiency(Form("/data7/htong/skim_samples/ele/skim_ele_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",(Int_t)x_mzh[i]),2,0);
    y_eyh_elc2[i] = getEleEfficiency(Form("/data7/htong/skim_samples/ele/skim_ele_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",(Int_t)x_mzh[i]),2,1);
    y_eyl_elc2[i] = getEleEfficiency(Form("/data7/htong/skim_samples/ele/skim_ele_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",(Int_t)x_mzh[i]),2,-1);
    /*
    y_eff_muc1[i] = getEleEfficiency(Form("/data7/htong/skim_samples/mu/skim_mu_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",(Int_t)x_mzh[i]),1,0);
    y_eyh_muc1[i] = getEleEfficiency(Form("/data7/htong/skim_samples/mu/skim_mu_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",(Int_t)x_mzh[i]),1,1);
    y_eyl_muc1[i] = getEleEfficiency(Form("/data7/htong/skim_samples/mu/skim_mu_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",(Int_t)x_mzh[i]),1,-1);

    y_eff_muc2[i] = getEleEfficiency(Form("/data7/htong/skim_samples/mu/skim_mu_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",(Int_t)x_mzh[i]),2,0);
    y_eyh_muc2[i] = getEleEfficiency(Form("/data7/htong/skim_samples/mu/skim_mu_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",(Int_t)x_mzh[i]),2,1);
    y_eyl_muc2[i] = getEleEfficiency(Form("/data7/htong/skim_samples/mu/skim_mu_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",(Int_t)x_mzh[i]),2,-1);
    */
  }

  TGraphAsymmErrors *g_eff_elc1 = new TGraphAsymmErrors(N, x_mzh, y_eff_elc1, 0, 0, y_eyl_elc1, y_eyh_elc1);
  TGraphAsymmErrors *g_eff_elc2 = new TGraphAsymmErrors(N, x_mzh, y_eff_elc2, 0, 0, y_eyl_elc2, y_eyh_elc2);
  //TGraphAsymmErrors *g_eff_muc1 = new TGraphAsymmErrors(N, x_mzh, y_eff_muc1, 0, 0, y_eyl_muc1, y_eyh_muc1);
  //TGraphAsymmErrors *g_eff_muc2 = new TGraphAsymmErrors(N, x_mzh, y_eff_muc2, 0, 0, y_eyl_muc2, y_eyh_muc2);

  g_eff_elc1->SetMinimum(0);
  g_eff_elc1->SetMaximum(0.5);
  g_eff_elc1->SetTitle("");
  g_eff_elc1->GetXaxis()->SetTitle("m_{ZH} (GeV)");
  g_eff_elc1->GetYaxis()->SetTitle("Acceptance #times efficiency");  
  g_eff_elc1->GetYaxis()->SetTitleOffset(1.3);

  g_eff_elc1->SetLineWidth(2);
  g_eff_elc1->SetLineColor(kBlue);
  g_eff_elc1->SetMarkerColor(kBlue);
  g_eff_elc1->SetMarkerStyle(8);

  g_eff_elc2->SetLineStyle(3);
  g_eff_elc2->SetLineWidth(2);
  g_eff_elc2->SetLineColor(kBlue);
  g_eff_elc2->SetMarkerColor(kBlue);
  g_eff_elc2->SetMarkerStyle(8);
  /*
  g_eff_muc1->SetLineWidth(2);
  g_eff_muc1->SetLineColor(kRed);
  g_eff_muc1->SetMarkerColor(kRed);
  g_eff_muc1->SetMarkerStyle(8);

  g_eff_muc2->SetLineStyle(3);
  g_eff_muc2->SetLineWidth(2);
  g_eff_muc2->SetLineColor(kRed);
  g_eff_muc2->SetMarkerColor(kRed);
  g_eff_muc2->SetMarkerStyle(8);
  */
  TLegend leg(0.73, 0.80, 0.90, 0.87);
  
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.04);
  leg.AddEntry(g_eff_elc1, "1 b-tag", "lp");
  leg.AddEntry(g_eff_elc2, "2 b-tag", "lp");
  //leg.AddEntry(g_eff_muc1, "1 b-tag", "lp");
  //leg.AddEntry(g_eff_muc2, "2 b-tag", "lp");

  TLatex lar;

  lar.SetNDC(kTRUE);
  lar.SetTextSize(0.04);
  lar.SetLineWidth(5);

  TCanvas c("c", "", 0, 0, 800, 600);

  c.cd();
  g_eff_elc1->Draw();
  g_eff_elc2->Draw("lpsame");
  //g_eff_muc1->Draw("lpsame");
  //g_eff_muc2->Draw("lpsame");
  leg.Draw();
  lar.DrawLatex(0.15, 0.83, "CMS #it{#bf{Simulation}}");
  c.Print("signalEfficiency.pdf");

}
