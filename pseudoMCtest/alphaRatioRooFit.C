#include <string>
#include <vector>
#include <iostream>
#include <TF1.h>
#include <TH1.h>
#include <TPad.h>
#include <TMath.h>
#include <TFile.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMatrixD.h>
#include <TFitResult.h>
#include <TSystemDirectory.h>
#include <TGraphAsymmErrors.h>
#include "../setNCUStyle.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooExponential.h"
#include "RooEffProd.h"
#include "RooFormulaVar.h"
#include "RooPlot.h"
#include "RooTFnBinding.h" 
#include "RooFFTConvPdf.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooChi2Var.h"
#include "RooCBShape.h"
#include "RooVoigtian.h"

R__LOAD_LIBRARY(PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(PDFs/PdfDiagonalizer_cc.so)

using namespace RooFit;

const Double_t xmin  = 500;
const Double_t xmax  = 5000;
const Int_t    nBins = (xmax-xmin)/100;

Double_t dataLumi  = 3000; //pb-1
Double_t xSecDY100 = 139.4*1.23;
Double_t xSecDY200 = 42.75*1.23;
Double_t xSecDY400 = 5.497*1.23;
Double_t xSecDY600 = 2.21*1.23;

TFile* getFile(std::string infiles, std::string hname, 
	       Double_t crossSection, Double_t* scale){

  TFile* f = TFile::Open(infiles.data());
  TH1D*  h = NULL;

  if( hname.find("pMC") != std::string::npos ) 
    h = (TH1D*)(f->Get("eventWeight_pMC"));
  else if( hname.find("pDA") != std::string::npos )
    h = (TH1D*)(f->Get("eventWeight_pDA"));

  *scale = dataLumi/(h->Integral()/crossSection);

  return f;

}

TH1D* addSamples(std::vector<string>& infiles, std::string hname,
		 TFile* f_DY100, TFile* f_DY200, TFile* f_DY400, TFile* f_DY600){ 

  Double_t scaleDY100 = 0;
  Double_t scaleDY200 = 0;
  Double_t scaleDY400 = 0;
  Double_t scaleDY600 = 0;

  for(unsigned int i = 0; i < infiles.size(); i++){

    if( infiles[i].find("HT-100") != std::string::npos )
      f_DY100 = getFile(infiles[i].data(), hname.data(), xSecDY100, &scaleDY100);

    if( infiles[i].find("HT-200") != std::string::npos )
      f_DY200 = getFile(infiles[i].data(), hname.data(), xSecDY200, &scaleDY200);

    if( infiles[i].find("HT-400") != std::string::npos )
      f_DY400 = getFile(infiles[i].data(), hname.data(), xSecDY400, &scaleDY400);

    if( infiles[i].find("HT-600") != std::string::npos )
      f_DY600 = getFile(infiles[i].data(), hname.data(), xSecDY600, &scaleDY600);

  }

  TH1D* DY100_temp = (TH1D*)(f_DY100->Get(Form("%s",hname.c_str())));
  TH1D* DY200_temp = (TH1D*)(f_DY200->Get(Form("%s",hname.c_str())));
  TH1D* DY400_temp = (TH1D*)(f_DY400->Get(Form("%s",hname.c_str())));
  TH1D* DY600_temp = (TH1D*)(f_DY600->Get(Form("%s",hname.c_str())));

  TH1D* h_Total = (TH1D*)(f_DY100->Get(Form("%s",hname.c_str())))->Clone("h_Total");

  h_Total->Reset();
  h_Total->Add(DY100_temp,scaleDY100);
  h_Total->Add(DY200_temp,scaleDY200);
  h_Total->Add(DY400_temp,scaleDY400);
  h_Total->Add(DY600_temp,scaleDY600);

  return h_Total;

}

void myRatio(TH1D* h_numer, TH1D* h_denom){

  TH1D* h_ratio = (TH1D*)h_numer->Clone("h_ratio");

  h_ratio->Reset();

  Int_t nbin = h_ratio->GetNbinsX();
  Double_t ratio[nbin];
  Double_t error[nbin];
  Double_t numer_nbincontent[nbin];
  Double_t denom_nbincontent[nbin];
  Double_t numer_binerror[nbin];
  Double_t denom_binerror[nbin];

  for( Int_t i = 1; i <= nbin; i++ ){

    numer_nbincontent[i] = h_numer->GetBinContent(i);
    denom_nbincontent[i] = h_denom->GetBinContent(i);
    numer_binerror[i] = h_numer->GetBinError(i);
    denom_binerror[i] = h_denom->GetBinError(i);

    if( denom_nbincontent[i] <= 0 || numer_nbincontent[i] <= 0 ) continue;
    if( denom_binerror[i] <= 0 || numer_binerror[i] <= 0 ) continue;

    ratio[i] = (Double_t)numer_nbincontent[i]/denom_nbincontent[i];
    error[i] = (ratio[i])*sqrt(pow(numer_binerror[i]/numer_nbincontent[i],2)+pow(denom_binerror[i]/denom_nbincontent[i],2));

    h_ratio->SetBinContent(i,ratio[i]);
    h_ratio->SetBinError(i,error[i]);

  }

  h_ratio->SetLineColor(kBlack);
  h_ratio->SetTitle("");
  h_ratio->GetYaxis()->SetTitle("Predicted/Truth");
  h_ratio->GetYaxis()->SetTitleOffset(0.3);
  h_ratio->GetXaxis()->SetTitle("ZH mass in signal region of pseudo-data");
  h_ratio->GetXaxis()->SetLabelSize(0.1);
  h_ratio->GetXaxis()->SetTitleSize(0.125);
  h_ratio->GetYaxis()->SetLabelSize(0.1);
  h_ratio->GetYaxis()->SetTitleSize(0.1);
  h_ratio->GetYaxis()->SetNdivisions(505);
  h_ratio->GetYaxis()->SetRangeUser(0,2);
  h_ratio->Draw();

  Double_t x0 = h_denom->GetXaxis()->GetXmin();
  Double_t x1 = h_denom->GetXaxis()->GetXmax();
  Double_t y0 = 1.;
  Double_t y1 = 1.;

  TLine* one = new TLine(x0,y0,x1,y1);

  one->SetLineColor(2);
  one->SetLineStyle(1);
  one->SetLineWidth(2);
  one->Draw("same");

  h_ratio->Draw("same");

}

void alphaRatioRooFit(std::string outputFolder){

  setNCUStyle();
  gStyle->SetOptFit(0);
  gStyle->SetMarkerSize(0);
  gStyle->SetTitleSize(0.04,"XYZ");
  gStyle->SetLabelSize(0.03,"XYZ");
  gStyle->SetHistLineWidth(2);

  std::vector<string> infiles;
 
  TSystemDirectory *base = new TSystemDirectory("root","root");
  base->SetDirectory(outputFolder.data());
  TList *listOfFiles = base->GetListOfFiles();
  TIter fileIt(listOfFiles);
  TFile *fileH = new TFile();
  Long64_t nfiles = 0;

  while( (fileH = (TFile*)fileIt()) ){
    
    std::string fileN = fileH->GetName();
    std::string baseString = "root";
    if( fileN.find(baseString) == std::string::npos ) continue;
    infiles.push_back(Form("%s/%s",outputFolder.data(),fileN.data()));
    nfiles++;
    
  }

  TFile *f_DY100 = NULL;
  TFile *f_DY200 = NULL;
  TFile *f_DY400 = NULL;
  TFile *f_DY600 = NULL;

  // Declare prefer histogram and add them together

  TH1D* h_sideTotalBKG  = addSamples(infiles,"ZprimeSide_pMC",f_DY100,f_DY200,f_DY400,f_DY600);
  TH1D* h_signTotalBKG  = addSamples(infiles,"ZprimeSign_pMC",f_DY100,f_DY200,f_DY400,f_DY600);
  TH1D* h_sideDATA      = addSamples(infiles,"ZprimeSide_pDA",f_DY100,f_DY200,f_DY400,f_DY600);
  TH1D* h_signDATA      = addSamples(infiles,"ZprimeSign_pDA",f_DY100,f_DY200,f_DY400,f_DY600);
  TH1D* h_corrPRmass    = addSamples(infiles,"corrPRmass_pDA",f_DY100,f_DY200,f_DY400,f_DY600);
  TH1D* h_corrPRmassAll = addSamples(infiles,"corrPRmassAll_pDA",f_DY100,f_DY200,f_DY400,f_DY600);

  h_sideTotalBKG->SetLineWidth(2);
  h_sideTotalBKG->SetLineColor(kBlack);
  h_sideTotalBKG->SetXTitle("ZH mass in side band of pseudo-MC");
  h_sideTotalBKG->SetYTitle("Event numbers");
  h_sideTotalBKG->SetTitleFont(62);

  h_signTotalBKG->SetLineWidth(2);
  h_signTotalBKG->SetLineColor(kBlack);
  h_signTotalBKG->SetXTitle("ZH mass in signal region of pseudo-MC");
  h_signTotalBKG->SetYTitle("Event numbers");
  h_signTotalBKG->SetTitleFont(62);

  h_signDATA->SetLineWidth(2);
  h_signDATA->SetLineColor(kBlue);
  h_signDATA->SetXTitle("ZH mass in signal region of pseudo-data");
  h_signDATA->SetYTitle("Event numbers");
  h_signDATA->SetTitleFont(62);

  h_corrPRmass->SetLineWidth(2);
  h_corrPRmass->SetLineColor(kBlack);
  h_corrPRmass->SetXTitle("Side band corrected pruned mass in pseudo-data");
  h_corrPRmass->SetYTitle("Event numbers");
  h_corrPRmass->SetTitleFont(62);

  h_corrPRmassAll->SetLineWidth(2);
  h_corrPRmassAll->SetLineColor(kBlack);
  h_corrPRmassAll->SetXTitle("Corrected pruned mass in pseudo-data");
  h_corrPRmassAll->SetYTitle("Event numbers");
  h_corrPRmassAll->SetTitleFont(62);

  // Make the statistics error more like data

  for( Int_t i = 1; i <= nBins; i++ ){

    h_sideDATA->SetBinError(i,TMath::Sqrt(h_sideDATA->GetBinContent(i)));
    h_signDATA->SetBinError(i,TMath::Sqrt(h_signDATA->GetBinContent(i)));

  }

  for( Int_t i = 1; i <= h_corrPRmass->GetNbinsX(); i++ ){

    h_corrPRmass   ->SetBinError(i,TMath::Sqrt(h_corrPRmass   ->GetBinContent(i)));
    h_corrPRmassAll->SetBinError(i,TMath::Sqrt(h_corrPRmassAll->GetBinContent(i)));

  }

  // Calculate alpha ratio

  TH1D* h_alphaRatio = new TH1D("h_alphaRatio", "", nBins, xmin, xmax); 

  h_alphaRatio->Sumw2();
  h_alphaRatio->SetXTitle("ZH mass");
  h_alphaRatio->SetYTitle("Alpha ratio");
  h_alphaRatio->Divide(h_signTotalBKG,h_sideTotalBKG);
  h_alphaRatio->SetMinimum(0);

  // Calculate number of backgrounds in signal region

  TH1D* h_numbkgDATA = (TH1D*)h_alphaRatio->Clone("h_numbkgDATA");

  h_numbkgDATA->Reset();

  for( Int_t i = 1; i <= nBins; i++ ){

    Double_t alphaRatio      = h_alphaRatio->GetBinContent(i); 
    Double_t sideDATA        = h_sideDATA->GetBinContent(i);
    Double_t numbkgDATA      = alphaRatio*sideDATA;      
    Double_t alphaRatioError = h_alphaRatio->GetBinError(i);
    Double_t sideDATAError   = h_sideDATA->GetBinError(i);

    if( alphaRatio == 0 || sideDATA == 0 ) continue;

    Double_t numbkgDATAError = numbkgDATA*sqrt(pow((alphaRatioError/alphaRatio),2)+pow((sideDATAError/sideDATA),2));

    h_numbkgDATA->SetBinContent(i,numbkgDATA);
    h_numbkgDATA->SetBinError(i,numbkgDATAError);

  }

  RooRealVar mass("mass", "M_{jet}", 0.0, 240.0, "GeV");

  mass.setRange("lowerSB", 40.0, 65.0);
  mass.setRange("upperSB", 145.0, 9999.0);
  mass.setRange("allRange", 40.0, 9999.0);

  RooDataHist corrPRmassAll("corrPRmassAll", "corrPRmassAll", mass, Import(*h_corrPRmassAll));
  RooDataHist corrPRmass("corrPRmass", "corrPRmass", mass, Import(*h_corrPRmass));

  RooPlot* corrPRmassAllframe;
  RooPlot* corrPRmassframe;

  corrPRmassAllframe = mass.frame(Title(""));
  corrPRmassframe    = mass.frame(Title(""));

  corrPRmassAll.plotOn(corrPRmassAllframe);
  corrPRmass.plotOn(corrPRmassframe);

  RooRealVar slope ("slope",  "slope of the exp", -0.03, -100.0, 0.0);
  RooRealVar offset("offset", "offset of the erf", 39,  0.0, 5000.0);
  RooRealVar width ("width",  "width of the erf",  30,  0.0, 5000.0);
  RooRealVar yieldSideband("yieldSideband", "Zjets normalization in sideband region", 10., 1., 1.e4);

  yieldSideband.setVal( (Int_t)h_corrPRmassAll->Integral() );
  yieldSideband.setError( TMath::Sqrt((Int_t)(h_corrPRmassAll->Integral())) );

  RooErfExpPdf mj_pdf("mj_pdf","fiting mj spectrum",mass,slope,offset,width);
  RooExtendPdf mj_pdf_ext("mj_pdf_ext","extended p.d.f", mj_pdf, yieldSideband);

  RooFitResult* result = mj_pdf_ext.fitTo(corrPRmassAll, Extended(kTRUE), Minimizer("Minuit"), Range("allRange"), Save());

  mj_pdf_ext.plotOn(corrPRmassAllframe, VisualizeError(*result,1), FillColor(kYellow));
  mj_pdf_ext.plotOn(corrPRmassAllframe);

  RooRealVar slope_sb ("slope_sb",  "slope of the exp",  slope.getVal(), -100.0, 0.0);
  RooRealVar offset_sb("offset_sb", "offset of the erf", offset.getVal(), 0.0, 5000.0);
  RooRealVar width_sb ("width_sb",  "width of the erf",  width.getVal(),  0.0, 5000.0);
  RooRealVar yieldSideband_sb("yieldSideband_sb", "Zjets normalization in sideband region", 10.0, 1.0, 1.e4);

  yieldSideband_sb.setVal( (Int_t)(h_corrPRmass->Integral()) );
  yieldSideband_sb.setError( TMath::Sqrt((Int_t)(h_corrPRmass->Integral())) );

  RooErfExpPdf mj_pdf_sb("mj_pdf_sb", "fiting mj spectrum", mass, slope_sb, offset_sb, width_sb);
  RooExtendPdf mj_pdf_ext_sb("mj_pdf_ext_sb", "extended p.d.f", mj_pdf_sb, yieldSideband_sb);

  mj_pdf_ext_sb.fitTo(corrPRmassAll, Extended(kTRUE), Minimizer("Minuit"), Range("lowerSB,upperSB"));
  mj_pdf_ext_sb.plotOn(corrPRmassAllframe, LineStyle(kDashed), LineColor(kRed));

  //h_numbkgDATA->Scale(nBkgSig/h_numbkgDATA->Integral(0,h_numbkgDATA->GetNbinsX()+1));

  // Output results

  TLegend* leg = new TLegend(0.21, 0.77, 0.87, 0.87);

  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);

  TLatex* lar = new TLatex();

  lar->SetTextSize(0.035);
  lar->SetLineWidth(5);

  TCanvas* c = new TCanvas("c","",0,0,1000,800);

  c->cd()->SetLogy(0);
  corrPRmassAllframe->Draw();
  //  lar->DrawLatexNDC(0.50, 0.60, "#font[22]{#color[4]{f(x) = #frac{1}{2} p_{0} e^{p_{1}x} ( 1 + erf ( #frac{x - p_{2}}{p_{3}} ) )}}");
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2015");
  lar->DrawLatexNDC(0.62, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  c->Print("alphaRatioRooFit.pdf(");

  c->cd()->SetLogy(0);
  corrPRmassframe->Draw();
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2015");
  lar->DrawLatexNDC(0.62, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  c->Print("alphaRatioRooFit.pdf)");

  /*
  c->cd()->SetLogy(0);
  h_corrPRmass->SetMaximum(300);
  h_corrPRmass->Draw();
  g_errorBands->Draw("3same");
  h_corrPRmass->Draw("same");
  leg->Clear();
  leg->AddEntry(h_corrPRmass, "Error = #sqrt{N_{per bin}}", "le");
  leg->AddEntry(g_errorBands, "Uncertainty based on fitting errors", "f");
  leg->Draw();
  lar->DrawLatexNDC(0.50, 0.65, Form("#chi^{2} / ndf: %f / %d",chisqr_cpm,ndf_cpm));
  lar->DrawLatexNDC(0.50, 0.55, "#font[22]{#color[4]{f(x) = #frac{1}{2} p_{0} e^{p_{1}x} ( 1 + erf ( #frac{x - p_{2}}{p_{3}} ) )}}");
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2015");
  lar->DrawLatexNDC(0.62, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  c->Print("alphaRatio.pdf");

  c->cd()->SetLogy();
  h_corrPRmass->SetMaximum(3e3);
  h_corrPRmass->Draw();
  g_errorBands->Draw("3same");
  h_corrPRmass->Draw("same");
  leg->Clear();
  leg->AddEntry(h_corrPRmass, "Error = #sqrt{N_{per bin}}", "le");
  leg->AddEntry(g_errorBands, "Uncertainty based on fitting errors", "f");
  leg->Draw();
  lar->DrawLatexNDC(0.25, 0.40, Form("#chi^{2} / ndf: %f / %d",chisqr_cpm,ndf_cpm));
  lar->DrawLatexNDC(0.25, 0.30, "#font[22]{#color[4]{f(x) = #frac{1}{2} p_{0} e^{p_{1}x} ( 1 + erf ( #frac{x - p_{2}}{p_{3}} ) )}}");
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2015");
  lar->DrawLatexNDC(0.62, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  c->Print("alphaRatio.pdf");
  
  c->cd()->SetLogy(0);
  h_signTotalBKG->Draw();
  lar->DrawLatexNDC(0.50, 0.80, Form("#chi^{2} / ndf: %f / %d",chisqr_sgb,ndf_sgb));
  lar->DrawLatexNDC(0.50, 0.70, "#font[22]{#color[4]{f_{signal}(x) = p_{0} e^{p_{1}x + #frac{p_{2}}{x}}}}");
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2015");
  lar->DrawLatexNDC(0.62, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  c->Print("alphaRatio.pdf");

  c->cd()->SetLogy(0);
  h_sideTotalBKG->Draw();
  lar->DrawLatexNDC(0.50, 0.80, Form("#chi^{2} / ndf: %f / %d",chisqr_sdb,ndf_sdb));
  lar->DrawLatexNDC(0.50, 0.70, "#font[22]{#color[4]{f_{side}(x) = p_{0} e^{p_{1}x + #frac{p_{2}}{x}}}}");
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2015");
  lar->DrawLatexNDC(0.62, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  c->Print("alphaRatio.pdf");

  c->cd()->SetLogy(0);
  h_alphaRatio->Draw();
  f_fitAlphaR->Draw("same");
  leg->Clear();
  leg->AddEntry(f_fitAlphaR, "#frac{f_{signal}(x)}{f_{side}(x)}", "l");
  leg->Draw();
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2015");
  lar->DrawLatexNDC(0.62, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  c->Print("alphaRatio.pdf");  

  c->Divide(1,2);
  TPad* c_up = (TPad*) c->GetListOfPrimitives()->FindObject("c_1");
  TPad* c_dw = (TPad*) c->GetListOfPrimitives()->FindObject("c_2"); 
  Double_t up_height = 0.8;
  Double_t dw_correction = 1.455;
  Double_t dw_height = (1-up_height)*dw_correction;
  c_up->SetPad(0,1-up_height,1,1);
  c_dw->SetPad(0,0,1,dw_height);
  c_dw->SetBottomMargin(0.25);
  c_up->cd();
  h_signDATA->GetXaxis()->SetTitle("");
  h_signDATA->GetXaxis()->SetLabelOffset(999);
  h_signDATA->GetXaxis()->SetLabelSize(0);
  h_signDATA->Draw();
  h_numbkgDATA->Draw("same");
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2015");
  lar->DrawLatexNDC(0.70, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  leg->Clear();
  leg->AddEntry(h_signDATA, "Truth backgrounds", "le");
  leg->AddEntry(h_numbkgDATA, "Predicted backgrounds", "le");
  leg->Draw();
  c_up->RedrawAxis();
  c_dw->cd();
  myRatio(h_numbkgDATA,h_signDATA);
  c->Draw();
  c->Print("alphaRatio.pdf)");
  */
}
