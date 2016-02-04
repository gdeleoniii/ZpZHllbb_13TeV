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
#include "../../setNCUStyle.h"
#include "addBackgrounds.h"
#include "fitFunction.h"
#include "uncertainty.h"
#include "myRatio.h"

const double xmin  = 500;
const double xmax  = 5000;
const int    nBins = (xmax-xmin)/100;

/// Main function start ///

void anotherfit(std::string outputFolder){

  setNCUStyle();
  gStyle->SetOptFit(0);
  gStyle->SetTitleSize(0.04,"XYZ");
  gStyle->SetLabelSize(0.03,"XYZ");

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
  TFile *f_TT    = NULL;
  TFile *f_WW    = NULL;
  TFile *f_WZ    = NULL;
  TFile *f_ZZ    = NULL;

  // Declare prefer histogram and add them together
  addBackgrounds dyjets;
  
  TH1D* h_PRmassBKG     = dyjets.add(infiles,"corrPRmassAll_pMC",f_DY100,f_DY200,f_DY400,f_DY600);
  TH1D* h_hollow_PRmass = dyjets.add(infiles,"corrPRmass_pDA",   f_DY100,f_DY200,f_DY400,f_DY600);
  TH1D* h_PRmass        = dyjets.add(infiles,"corrPRmassAll_pDA",f_DY100,f_DY200,f_DY400,f_DY600);
  /*
  TH1D* h_PRmassBKG     = addMinorBkg(infiles,"corrPRmassAll_pMC",f_TT,f_WW,f_WZ,f_ZZ);
  TH1D* h_hollow_PRmass = addMinorBkg(infiles,"corrPRmass_pDA",   f_TT,f_WW,f_WZ,f_ZZ);
  TH1D* h_PRmass        = addMinorBkg(infiles,"corrPRmassAll_pDA",f_TT,f_WW,f_WZ,f_ZZ);
  */
  h_PRmassBKG->SetMarkerStyle(8);
  h_PRmassBKG->SetMarkerSize(1.5);
  h_PRmassBKG->SetLineColor(kBlack);
  h_PRmassBKG->SetXTitle("Side band corrected pruned mass in MC");
  h_PRmassBKG->SetYTitle("Event numbers");
  h_PRmassBKG->SetTitleFont(62);

  h_hollow_PRmass->SetMarkerStyle(8);
  h_hollow_PRmass->SetMarkerSize(1.5);
  h_hollow_PRmass->SetLineColor(kBlack);
  h_hollow_PRmass->SetXTitle("Side band corrected pruned mass in pseudo-data");
  h_hollow_PRmass->SetYTitle("Event numbers");
  h_hollow_PRmass->SetTitleFont(62);

  h_PRmass->SetMarkerStyle(8);
  h_PRmass->SetMarkerSize(1.5);
  h_PRmass->SetLineColor(kBlack);
  h_PRmass->SetXTitle("Corrected pruned mass in pseudo-data");
  h_PRmass->SetYTitle("Event numbers");
  h_PRmass->SetTitleFont(62);

  // Make the statistics error more like data

  for( int i = 1; i <= h_PRmass->GetNbinsX(); i++ ){

    h_PRmassBKG->SetBinError(i,TMath::Sqrt(h_PRmassBKG->GetBinContent(i)));
    h_hollow_PRmass->SetBinError(i,TMath::Sqrt(h_hollow_PRmass->GetBinContent(i)));
    h_PRmass->SetBinError(i,TMath::Sqrt(h_PRmass->GetBinContent(i)));

  }

  // Fit pruned mass with signal region in MC (BKG)

  TF1* f_fitPRmassBKG = new TF1("f_fitPRmassBKG", fitPRmass, 40, 240, 5);

  f_fitPRmassBKG->SetLineWidth(2);
  f_fitPRmassBKG->SetLineColor(kBlue);

  double parFitPRm[4] = {1224,-0.107,139.6,107.4};

  f_fitPRmassBKG->SetParameters(parFitPRm[0],parFitPRm[1],parFitPRm[2],parFitPRm[3],h_PRmassBKG->GetBinWidth(1));
  f_fitPRmassBKG->FixParameter(4,h_PRmassBKG->GetBinWidth(1));
  f_fitPRmassBKG->FixParameter(0,h_PRmassBKG->Integral());

  h_PRmassBKG->Fit("f_fitPRmassBKG", "Q", "", 40, 240);

  // Fit pruned mass with signal region  

  TF1* f_fitPRmass = new TF1("f_fitPRmass", fitPRmass, 40, 240, 5);

  f_fitPRmass->SetLineWidth(2);
  f_fitPRmass->SetLineColor(kBlue);

  f_fitPRmass->FixParameter(1,f_fitPRmassBKG->GetParameter(1));
  f_fitPRmass->FixParameter(2,f_fitPRmassBKG->GetParameter(2));
  f_fitPRmass->FixParameter(3,f_fitPRmassBKG->GetParameter(3));
  f_fitPRmass->FixParameter(4,h_PRmass->GetBinWidth(1));

  h_PRmass->Fit("f_fitPRmass", "Q", "", 40, 240);

  double chisqr_cpma = f_fitPRmass->GetChisquare();
  int ndf_cpma = f_fitPRmass->GetNDF();

  TFitResultPtr fitptr = h_PRmass->Fit(f_fitPRmass, "QS");
  TFitResult fitresult = (*fitptr);
  TMatrixD corrMatrix  = fitresult.GetCorrelationMatrix();

  double dummy = 0.0;

  TGraphAsymmErrors* g_errorBands = fitUncertainty(true, f_fitPRmass, &corrMatrix, fitPRmass, NULL, 0, &dummy, &dummy);

  g_errorBands->SetFillStyle(1001);
  g_errorBands->SetFillColor(kYellow);

  // Fit pruned mass without signal region 

  TF1* f_hollow_fitPRmass = new TF1("f_hollow_fitPRmass", hollow_fitPRmass, 40, 240, 5);

  f_hollow_fitPRmass->SetLineWidth(2);
  f_hollow_fitPRmass->SetLineColor(kBlue);
 
  f_hollow_fitPRmass->FixParameter(1,f_fitPRmassBKG->GetParameter(1));
  f_hollow_fitPRmass->FixParameter(2,f_fitPRmassBKG->GetParameter(2));
  f_hollow_fitPRmass->FixParameter(3,f_fitPRmassBKG->GetParameter(3));
  f_hollow_fitPRmass->FixParameter(4,h_hollow_PRmass->GetBinWidth(1));
 
  h_hollow_PRmass->Fit("f_hollow_fitPRmass", "Q", "", 40, 240);

  double chisqr_cpm = f_hollow_fitPRmass->GetChisquare();
  int ndf_cpm = f_hollow_fitPRmass->GetNDF();

  TFitResultPtr hollow_fitptr = h_hollow_PRmass->Fit(f_hollow_fitPRmass, "QS");
  TFitResult hollow_fitresult = (*hollow_fitptr);
  TMatrixD hollow_corrMatrix  = hollow_fitresult.GetCorrelationMatrix();

  TGraphAsymmErrors* g_hollow_errorBands = fitUncertainty(true, f_hollow_fitPRmass, &hollow_corrMatrix, hollow_fitPRmass, NULL, 0, &dummy, &dummy);

  g_hollow_errorBands->SetFillStyle(1001);
  g_hollow_errorBands->SetFillColor(kYellow);


  //// Fluctuate the pruned mass histogram to test the fitting results ////       

  TH1D* h_fluc = (TH1D*)h_PRmass->Clone("h_fluc");
  TH1D* h_hollow_fluc = (TH1D*)h_PRmass->Clone("h_hollow_fluc");
  TH1D* h_bias = new TH1D("h_bias", "", 100, -0.5, 0.5);

  TCanvas* ctemp = new TCanvas("ctemp","",0,0,1000,800);
  TFile* outemp = new TFile("anotherfit.root","recreate");

  h_bias->SetLineWidth(1);
  h_bias->SetFillColor(kYellow);
  h_bias->GetXaxis()->SetTitle("Bias ((fit-true)/true)");
  h_bias->GetYaxis()->SetTitle("Counts");

  for( int ntoy = 0; ntoy < 1000; ntoy++ ){
    if(ntoy%10==1)cout<<ntoy<<endl;
    h_fluc->Reset();
    h_hollow_fluc->Reset();

    for( int idata = 0; idata < (int)(h_PRmass->Integral()); idata++ ){

      h_fluc->Fill(f_fitPRmass->GetRandom(40,240));
     
    }

    int nd1 = (105 - h_fluc->GetBinLowEdge(1))/h_fluc->GetBinWidth(1)+1;
    int nd2 = (135 - h_fluc->GetBinLowEdge(1))/h_fluc->GetBinWidth(1);

    double nSigHist = h_fluc->Integral(nd1,nd2);

    for( int nbin = 1; nbin <= h_fluc->GetNbinsX(); nbin++ ){

      if( h_hollow_fluc->GetBinLowEdge(nbin) < 65 || h_hollow_fluc->GetBinLowEdge(nbin) > 140 ){

	h_hollow_fluc->SetBinContent(nbin,h_fluc->GetBinContent(nbin));
	h_hollow_fluc->SetBinError(nbin,h_fluc->GetBinError(nbin));

      }

    }

    TF1* f_fluc = new TF1("f_fluc", hollow_fitPRmass, 40, 240, 5);

    f_fluc->FixParameter(1,f_fitPRmassBKG->GetParameter(1));
    f_fluc->FixParameter(2,f_fitPRmassBKG->GetParameter(2));
    f_fluc->FixParameter(3,f_fitPRmassBKG->GetParameter(3));
    f_fluc->FixParameter(4,h_hollow_fluc->GetBinWidth(1));

    h_hollow_fluc->Fit("f_fluc", "Q", "", 40, 240);

    // temp region
    ctemp->cd();
    h_fluc->SetLineColor(kRed);
    h_fluc->Draw();
    h_hollow_fluc->Draw("same");
    ctemp->Write();
    //

    double nSigFit = f_fluc->Integral(105.00000000,135.00000000)/h_hollow_fluc->GetBinWidth(1);

    h_bias->Fill((nSigFit - nSigHist)/nSigHist);

  }
  outemp->Write();
  TF1* f_bias = new TF1("f_bias","gaus");
  f_bias->SetLineWidth(2);
  f_bias->SetLineColor(kBlue);
  h_bias->Fit("f_bias","Q");

  //// End of the test ////  

  // Output results

  TLegend* leg = new TLegend(0.35, 0.77, 0.87, 0.87);

  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);

  TLatex* lar = new TLatex();

  lar->SetTextSize(0.035);
  lar->SetLineWidth(5);

  TCanvas* c = new TCanvas("c","",0,0,1000,800);

  c->cd()->SetLogy(0);
  h_PRmass->Draw();
  g_errorBands->Draw("3same");
  h_PRmass->Draw("same");
  leg->Clear();
  leg->AddEntry(h_PRmass, "Error = #sqrt{N_{per bin}}", "lp");
  leg->AddEntry(g_errorBands, "Uncertainty based on fitting errors", "f");
  leg->Draw();
  lar->DrawLatexNDC(0.50, 0.65, Form("#chi^{2} / ndf: %f / %d",chisqr_cpma,ndf_cpma));
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2016");
  lar->DrawLatexNDC(0.65, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  c->Print("anotherfita.pdf(");

  c->cd()->SetLogy(0);
  h_PRmassBKG->Draw();
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2016");
  lar->DrawLatexNDC(0.65, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  c->Print("anotherfita.pdf");

  c->cd()->SetLogy(0);
  h_hollow_PRmass->Draw();
  g_hollow_errorBands->Draw("3same");
  h_hollow_PRmass->Draw("same");
  leg->Clear();
  leg->AddEntry(h_hollow_PRmass, "Error = #sqrt{N_{per bin}}", "lp");
  leg->AddEntry(g_errorBands, "Uncertainty based on fitting errors", "f");
  leg->Draw();
  lar->DrawLatexNDC(0.50, 0.65, Form("#chi^{2} / ndf: %f / %d",chisqr_cpm,ndf_cpm));
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2016");
  lar->DrawLatexNDC(0.65, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  c->Print("anotherfita.pdf");
  
  c->cd()->SetLogy();
  h_hollow_PRmass->SetMaximum(3e3);
  h_hollow_PRmass->Draw();
  g_hollow_errorBands->Draw("3same");
  h_hollow_PRmass->Draw("same");
  leg->Clear();
  leg->AddEntry(h_hollow_PRmass, "Error = #sqrt{N_{per bin}}", "lp");
  leg->AddEntry(g_hollow_errorBands, "Uncertainty based on fitting errors", "f");
  leg->Draw();
  lar->DrawLatexNDC(0.25, 0.40, Form("#chi^{2} / ndf: %f / %d",chisqr_cpm,ndf_cpm));
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2016");
  lar->DrawLatexNDC(0.65, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  c->Print("anotherfita.pdf");

  c->cd()->SetLogy(0);
  h_bias->Draw();
  lar->DrawLatexNDC(0.15, 0.94, "CMS work in progress");
  lar->DrawLatexNDC(0.65, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  c->Print("anotherfita.pdf)");
  
}
