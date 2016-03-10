#include <string>
#include <iostream>
#include <TF1.h>
#include <TH1.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include "../setNCUStyle.h"

void fitDraw(TH1D* h_prmass, TH1D* h_prmass_hollow,
	     TGraphAsymmErrors* g_errorBands, TGraphAsymmErrors* g_errorBands_hollow,
	     TH1D* h_bias, TH1D* h_upPull, TH1D* h_dwPull, 
	     std::string pdfName){

  setNCUStyle();
  gStyle->SetOptFit(0);
  gStyle->SetTitleSize(0.04,"XYZ");
  gStyle->SetLabelSize(0.03,"XYZ");

  h_prmass->SetMarkerStyle(8);
  h_prmass->SetMarkerSize(1.5);
  h_prmass->SetLineColor(kBlack);
  h_prmass->SetXTitle("Corrected pruned jet mass");
  h_prmass->SetYTitle("Event numbers");
  h_prmass->SetTitleFont(62); 

  h_prmass_hollow->SetMarkerStyle(8);
  h_prmass_hollow->SetMarkerSize(1.5);
  h_prmass_hollow->SetLineColor(kBlack);
  h_prmass_hollow->SetXTitle("Corrected pruned jet mass without signal region");
  h_prmass_hollow->SetYTitle("Event numbers");
  h_prmass_hollow->SetTitleFont(62);

  h_bias->SetLineWidth(1);
  h_bias->SetFillColor(kYellow);
  h_bias->GetXaxis()->SetTitle("Bias ((fit-true)/true)");

  h_upPull->SetLineWidth(1);
  h_upPull->SetLineColor(kRed);
  h_upPull->GetXaxis()->SetTitle("Pull ((fit-true)/unc)");

  h_dwPull->SetLineWidth(1);
  h_dwPull->SetLineColor(kBlue);
  h_dwPull->GetXaxis()->SetTitle("Pull ((fit-true)/unc)");

  TLegend* leg = new TLegend(0.35, 0.77, 0.87, 0.87);

  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);

  TCanvas* c = new TCanvas("c","",0,0,1000,800);
 
  TLatex* lar = new TLatex();

  lar->SetTextSize(0.035);
  lar->SetLineWidth(5);
  //lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2016");
  //lar->DrawLatexNDC(0.65, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
 
  c->cd();
  h_prmass->Draw();
  g_errorBands->Draw("3same");
  h_prmass->Draw("same");
  leg->AddEntry(h_prmass, "Error = #sqrt{N_{per bin}}", "lp");
  leg->AddEntry(g_errorBands, "Uncertainty based on fitting errors", "f");
  leg->Draw();
  lar->DrawLatexNDC(0.50, 0.65, Form("# of events : %f", h_prmass->Integral()));
  // lar->DrawLatexNDC(0.50, 0.65, Form("#chi^{2}/ndf: %f/%d", f_fitprmass->GetChisquare(), f_fitprmass->GetNDF()));
  c->Print(Form("%s.pdf(",pdfName.data()));
  
  leg->Clear();

  c->cd();
  h_prmass_hollow->Draw();
  g_errorBands_hollow->Draw("3same");
  h_prmass_hollow->Draw("same");
  leg->AddEntry(h_prmass_hollow, "Error = #sqrt{N_{per bin}}", "lp");
  leg->AddEntry(g_errorBands_hollow, "Uncertainty based on fitting errors", "f");
  leg->Draw();
  lar->DrawLatexNDC(0.50, 0.65, Form("# of events : %f", h_prmass_hollow->Integral()));
  //lar->DrawLatexNDC(0.50, 0.65, Form("#chi^{2}/ndf: %f/%d", f_fitmass_hollow->GetChisquare(), f_fitprmass_hollow->GetNDF()));
  c->Print(Form("%s.pdf",pdfName.data()));
  
  leg->Clear();

  c->cd();
  h_bias->Draw();
  c->Print(Form("%s.pdf",pdfName.data()));

  c->cd();
  h_upPull->Draw();
  h_dwPull->Draw("same");
  leg->AddEntry(h_upPull, "positive pull", "l");
  leg->AddEntry(h_dwPull, "negative pull", "l");
  leg->Draw();
  c->Print(Form("%s.pdf)",pdfName.data()));
  
}
