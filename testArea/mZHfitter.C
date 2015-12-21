#include <string>
#include <vector>
#include <iostream>
#include <TF1.h>
#include <TH1.h>
#include <TKey.h>
#include <TMath.h>
#include <TFile.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TCanvas.h>
#include <TLegend.h>
#include "../setNCUStyle.h"

double fitZHmass(double* v, double* p){

  double x = v[0];
  return p[0]*TMath::Exp(p[1]*x + p[2]/x);

}

void mZHfitter(std::string input, std::string output){

  setNCUStyle();
  gStyle->SetOptFit(0);
  gStyle->SetTitleSize(0.04,"XYZ");
  gStyle->SetLabelSize(0.03,"XYZ");

  // To get the name of histograms
  
  TFile *f_ = TFile::Open(input.data());
  f_->cd();
  
  TDirectory *current_sourcedir = gDirectory;
  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key;

  vector<std::string> h_name;

  while ( (key = (TKey*)nextkey()) ) {

    TObject *obj = key->ReadObj();

    if ( obj->IsA()->InheritsFrom("TH1") ) 
      h_name.push_back(obj->GetName());

  }

  TLegend* leg = new TLegend(0.35, 0.77, 0.87, 0.87);

  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);

  TLatex* lar = new TLatex();

  lar->SetTextSize(0.035);
  lar->SetLineWidth(5);

  TCanvas* c = new TCanvas("c","",0,0,1000,800);

  std::string region;

  if( input.find("side") != std::string::npos ) region = "side band";
  else if( input.find("signal") != std::string::npos) region = "signal region";

  int xmin = 700; 
  int xmax = 4500;

  for(unsigned int i = 0; i < h_name.size(); i++){

    TH1D *h_ZHmass = (TH1D*)(f_->Get(h_name[i].data()));

    h_ZHmass->SetAxisRange(xmin,xmax,"X");
    h_ZHmass->SetMarkerStyle(8);
    h_ZHmass->SetMarkerSize(1.5);
    h_ZHmass->SetLineColor(kBlack);
    h_ZHmass->SetXTitle(Form("ZH mass in %s",region.data()));
    h_ZHmass->SetYTitle("Event numbers");
    h_ZHmass->SetTitleFont(62);

    // Fit ZH mass

    TF1* f_fitZHmass = new TF1("f_fitZHmass", fitZHmass, xmin, xmax, 3);

    f_fitZHmass->SetLineWidth(2);
    f_fitZHmass->SetLineColor(kBlue);
    f_fitZHmass->SetParameters(0,0,0);
    h_ZHmass->Fit("f_fitZHmass", "Q", "", xmin, xmax);

    double chisqr = f_fitZHmass->GetChisquare();
    int ndf = f_fitZHmass->GetNDF();

    std::string cutvalue = h_name[i].substr(21,3);

    c->cd(); 
    h_ZHmass->Draw();
    leg->Clear();
    leg->AddEntry(h_ZHmass, Form("|#Delta #eta_{ZH}| < %s", cutvalue.data()), "lp");
    leg->Draw();
    lar->DrawLatexNDC(0.50, 0.65, Form("#chi^{2} / ndf: %f / %d", chisqr, ndf));
    lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2015"); 
    lar->DrawLatexNDC(0.65, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");

    std::string pdfname;

    if( i == 0 ) pdfname = output + "(";
    else if( i == h_name.size()-1 ) pdfname = output + ")";
    else pdfname = output;

    c->Print(pdfname.data());

  }

}
