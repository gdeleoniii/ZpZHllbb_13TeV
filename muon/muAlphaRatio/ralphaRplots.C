#include <string>
#include <vector>
#include <iostream>
#include <TF1.h>
#include <TH1.h>
#include <TMath.h>
#include <TFile.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystemDirectory.h>
#include "../../setNCUStyle.h"

const Double_t varBins[] = {600,800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500};
Int_t nvarBins = sizeof(varBins)/sizeof(varBins[0])-1;

Double_t dataLumi  = 831.7; //pb-1
Double_t xSecDY100 = 139.4*1.23;
Double_t xSecDY200 = 42.75*1.23;
Double_t xSecDY400 = 5.497*1.23;
Double_t xSecDY600 = 2.21*1.23;
Double_t xSecTTbar = 831.76;
Double_t xSecWW    = 118.7;
Double_t xSecWZ    = 47.13;
Double_t xSecZZ    = 16.523;

TFile* getFile(std::string infiles, Double_t crossSection, Double_t* scale){

  TFile* f = TFile::Open(infiles.data());
  TH1D*  h = (TH1D*)(f->Get("eventWeight"));

  *scale = dataLumi/(h->Integral()/crossSection);

  return f;

}

TH1D* fixDiscdBin(TH1D* h){

  for( Int_t nb = 0; nb < nvarBins; nb++ ){
    Int_t binRatio = (varBins[nb+1]-varBins[nb])/(varBins[1]-varBins[0]);
    h->SetBinContent(nb+1,(h->GetBinContent(nb+1)/binRatio));
    h->SetBinError(nb+1,(h->GetBinError(nb+1)/binRatio));
  }

  return h;

}

TH1D* addSamples(std::vector<string>& infiles, std::string hname,
		 TFile* f_DY100, TFile* f_DY200, TFile* f_DY400, TFile* f_DY600, 
		 TFile* f_TTbar, TFile* f_WW, TFile* f_WZ, TFile* f_ZZ){

  Double_t scaleDY100 = 0;
  Double_t scaleDY200 = 0;
  Double_t scaleDY400 = 0;
  Double_t scaleDY600 = 0;
  Double_t scaleTTbar = 0;
  Double_t scaleWW    = 0;
  Double_t scaleWZ    = 0;
  Double_t scaleZZ    = 0;

  for(unsigned int i = 0; i < infiles.size(); i++){

    if( infiles[i].find("HT-100") != std::string::npos )
      f_DY100 = getFile(infiles[i].data(), xSecDY100, &scaleDY100);

    if( infiles[i].find("HT-200") != std::string::npos )
      f_DY200 = getFile(infiles[i].data(), xSecDY200, &scaleDY200);

    if( infiles[i].find("HT-400") != std::string::npos )
      f_DY400 = getFile(infiles[i].data(), xSecDY400, &scaleDY400);

    if( infiles[i].find("HT-600") != std::string::npos )
      f_DY600 = getFile(infiles[i].data(), xSecDY600, &scaleDY600);

    if( infiles[i].find("TT_") != std::string::npos )
      f_TTbar = getFile(infiles[i].data(), xSecTTbar, &scaleTTbar);

    if( infiles[i].find("WW_") != std::string::npos )
      f_WW = getFile(infiles[i].data(), xSecWW, &scaleWW);

    if( infiles[i].find("WZ_") != std::string::npos )
      f_WZ = getFile(infiles[i].data(), xSecWZ, &scaleWZ);

    if( infiles[i].find("ZZ_") != std::string::npos )
      f_ZZ = getFile(infiles[i].data(), xSecZZ, &scaleZZ);

  }

  TH1D* DY100_temp = fixDiscdBin((TH1D*)(f_DY100->Get(Form("%s",hname.c_str()))));
  TH1D* DY200_temp = fixDiscdBin((TH1D*)(f_DY200->Get(Form("%s",hname.c_str()))));
  TH1D* DY400_temp = fixDiscdBin((TH1D*)(f_DY400->Get(Form("%s",hname.c_str()))));
  TH1D* DY600_temp = fixDiscdBin((TH1D*)(f_DY600->Get(Form("%s",hname.c_str()))));
  TH1D* TTbar_temp = fixDiscdBin((TH1D*)(f_TTbar->Get(Form("%s",hname.c_str()))));
  TH1D* WW_temp    = fixDiscdBin((TH1D*)(f_WW->Get(Form("%s",hname.c_str()))));
  TH1D* WZ_temp    = fixDiscdBin((TH1D*)(f_WZ->Get(Form("%s",hname.c_str()))));
  TH1D* ZZ_temp    = fixDiscdBin((TH1D*)(f_ZZ->Get(Form("%s",hname.c_str()))));

  TH1D* h_Total = (TH1D*)(f_DY100->Get(Form("%s",hname.c_str())))->Clone("h_Total");

  h_Total->Reset();
  h_Total->Add(DY100_temp,scaleDY100);
  h_Total->Add(DY200_temp,scaleDY200);
  h_Total->Add(DY400_temp,scaleDY400);
  h_Total->Add(DY600_temp,scaleDY600);
  h_Total->Add(TTbar_temp,scaleTTbar);
  h_Total->Add(WW_temp,scaleWW);
  h_Total->Add(WZ_temp,scaleWZ);
  h_Total->Add(ZZ_temp,scaleZZ);

  return h_Total;

}

TH1D* addData(std::vector<string>& infiles, std::string hname,
	      TFile* f_data0, TFile* f_data1){

  for(unsigned int i = 0; i < infiles.size(); i++){

    if( infiles[i].find("V3_2015") != std::string::npos )
      f_data0 = TFile::Open(infiles[i].data());

    if( infiles[i].find("V4_2015") != std::string::npos )
      f_data1 = TFile::Open(infiles[i].data());

  }

  TH1D* data0_temp = fixDiscdBin((TH1D*)(f_data0->Get(Form("%s",hname.c_str()))));
  TH1D* data1_temp = fixDiscdBin((TH1D*)(f_data1->Get(Form("%s",hname.c_str()))));

  TH1D* h_Total = (TH1D*)(f_data0->Get(Form("%s",hname.c_str())))->Clone("h_Total");

  h_Total->Reset();
  h_Total->Add(data0_temp);
  h_Total->Add(data1_temp);

  return h_Total;

}

Double_t fitZprime(Double_t* v, Double_t* p){

  Double_t x = v[0];
  return p[0]*TMath::Exp(p[1]*x + p[2]/x);

}

Double_t divFunc(Double_t* v, Double_t* p){

  Double_t x = v[0];
  return (p[0]*TMath::Exp(p[1]*x+p[2]/x))/(p[3]*TMath::Exp(p[4]*x+p[5]/x));

}

Double_t fitPRmass(Double_t* v, Double_t* p){
  
  Double_t x = v[0];
  return p[0]*TMath::Exp(p[1]*x)*0.5*(1+TMath::Erf((x-p[2])/p[3]));
  
}

void ralphaRplots(std::string outputFolder){

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
  TFile *f_TTbar = NULL;
  TFile *f_WW    = NULL;
  TFile *f_WZ    = NULL;
  TFile *f_ZZ    = NULL;
  TFile *f_data0 = NULL;
  TFile *f_data1 = NULL;

  // Declare prefer histogram and add them together

  TH1D* h_sideDATA      = addData(infiles,"ZprimeSide",f_data0,f_data1);
  TH1D* h_signDATA      = addData(infiles,"ZprimeSign",f_data0,f_data1);
  TH1D* h_corrPRmass    = addData(infiles,"corrPRmass",f_data0,f_data1);
  TH1D* h_corrPRmassAll = addData(infiles,"corrPRmassAll",f_data0,f_data1);
  TH1D* h_sideTotalBKG  = addSamples(infiles,"ZprimeSide",f_DY100,f_DY200,f_DY400,f_DY600,f_TTbar,f_WW,f_WZ,f_ZZ);
  TH1D* h_signTotalBKG  = addSamples(infiles,"ZprimeSign",f_DY100,f_DY200,f_DY400,f_DY600,f_TTbar,f_WW,f_WZ,f_ZZ);

  h_sideTotalBKG->SetLineWidth(2);
  h_sideTotalBKG->SetLineColor(kBlack);
  h_sideTotalBKG->SetXTitle("Side band ZH mass in MC");
  h_sideTotalBKG->SetYTitle("Event numbers");

  h_signTotalBKG->SetLineWidth(2);
  h_signTotalBKG->SetLineColor(kBlack);
  h_signTotalBKG->SetXTitle("Signal region ZH mass in MC");
  h_signTotalBKG->SetYTitle("Event numbers");

  h_signDATA->SetLineWidth(2);
  h_signDATA->SetLineColor(kBlue);
  h_signDATA->SetXTitle("ZH mass");
  h_signDATA->SetYTitle("Event numbers");

  h_corrPRmass->SetLineWidth(2);
  h_corrPRmass->SetLineColor(kBlack);
  h_corrPRmass->SetXTitle("Side band corrected pruned mass in data");
  h_corrPRmass->SetYTitle("Event numbers");

  h_corrPRmassAll->SetLineWidth(2);
  h_corrPRmassAll->SetLineColor(kBlack);
  h_corrPRmassAll->SetXTitle("Corrected pruned mass in data");
  h_corrPRmassAll->SetYTitle("Event numbers");

  // Calculate alpha ratio

  TH1D* h_alphaRatio   = new TH1D("h_alphaRatio", "", nvarBins, varBins);

  h_alphaRatio->Sumw2();
  h_alphaRatio->SetXTitle("ZH mass");
  h_alphaRatio->SetYTitle("Alpha ratio");
  h_alphaRatio->Divide(h_signTotalBKG,h_sideTotalBKG);
  h_alphaRatio->SetMinimum(0);

  // Calculate number of backgrounds in signal region

  TH1D* h_numbkgDATA = (TH1D*)h_alphaRatio->Clone("h_numbkgDATA");

  h_numbkgDATA->Reset();

  for( Int_t i = 1; i <= nvarBins; i++ ){

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

  // Fitting procedure

  TF1* f_fitPRmass    = new TF1("f_fitPRmass", fitPRmass, 40, 240, 4);
  TF1* f_fitPRmassAll = new TF1("f_fitPRmassAll", fitPRmass, 40, 240, 4);
  TF1* f_fitSideBkg   = new TF1("f_fitSideBkg", fitZprime, varBins[0], varBins[nvarBins], 3);
  TF1* f_fitSignBkg   = new TF1("f_fitSignBkg", fitZprime, varBins[0], varBins[nvarBins], 3);
  TF1* f_fitAlphaR    = new TF1("f_fitAlphaR", divFunc, varBins[0], varBins[nvarBins], 6);

  h_corrPRmass   ->Fit("f_fitPRmass", "", "", 40, 240);
  h_corrPRmassAll->Fit("f_fitPRmassAll", "", "", 40, 240);

  h_sideTotalBKG->Fit("f_fitSideBkg", "", "", varBins[0], varBins[nvarBins]);
  h_signTotalBKG->Fit("f_fitSignBkg", "", "", varBins[0], varBins[nvarBins]);

  f_fitAlphaR->SetParameter(0, f_fitSignBkg->GetParameter(0));
  f_fitAlphaR->SetParameter(1, f_fitSignBkg->GetParameter(1));
  f_fitAlphaR->SetParameter(2, f_fitSignBkg->GetParameter(2));
  f_fitAlphaR->SetParameter(3, f_fitSideBkg->GetParameter(0));
  f_fitAlphaR->SetParameter(4, f_fitSideBkg->GetParameter(1));
  f_fitAlphaR->SetParameter(5, f_fitSideBkg->GetParameter(2));

  f_fitPRmass ->SetLineWidth(2);
  f_fitSideBkg->SetLineWidth(2);
  f_fitSignBkg->SetLineWidth(2);
  f_fitAlphaR ->SetLineWidth(2);

  Double_t nBkgSig = f_fitPRmass->Integral(105,135)/h_corrPRmass->GetBinWidth(1);
  h_numbkgDATA->Scale(nBkgSig/h_numbkgDATA->Integral());
  h_numbkgDATA->SetMinimum(0);

  // Output results

  TLegend *leg = new TLegend(0.30, 0.70, 0.85, 0.80);

  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(h_signDATA, "signal region of data", "le");
  leg->AddEntry(h_numbkgDATA, "backgrounds in signal region of data", "le");

  TLatex *lar = new TLatex(0.50, 0.94, "CMS,  #sqrt{s} = 13 TeV, L = 831.7 pb^{-1}");

  lar->SetNDC(kTRUE);
  lar->SetTextSize(0.04);
  lar->SetLineWidth(5);

  TLatex *eqn = new TLatex(0.50, 0.80, "#font[22]{#color[4]{f(x) = #frac{1}{2} p_{0} e^{p_{1}x} ( 1 + erf ( #frac{x - p_{2}}{p_{3}} ) )}}");

  eqn->SetNDC(kTRUE);
  eqn->SetTextSize(0.04);

  TCanvas* c = new TCanvas("c","",0,0,1000,800);
  
  c->cd();
  h_corrPRmassAll->Draw();
  f_fitPRmassAll ->Draw("same");
  lar->Draw();
  eqn->Draw();
  c->Print("realAlphaRatio.pdf(");

  c->cd();
  h_corrPRmass->Draw();
  f_fitPRmass ->Draw("same");
  lar->Draw();
  eqn->Draw();
  c->Print("realAlphaRatio.pdf");
  
  c->cd();
  h_signTotalBKG->Draw();
  f_fitSignBkg  ->Draw("same");
  lar->Draw();
  c->Print("realAlphaRatio.pdf");

  c->cd();
  h_sideTotalBKG->Draw();
  f_fitSideBkg  ->Draw("same");
  lar->Draw();
  c->Print("realAlphaRatio.pdf");

  c->cd();
  h_alphaRatio->Draw();
  f_fitAlphaR ->Draw("same");
  lar->Draw();
  c->Print("realAlphaRatio.pdf");  

  c->cd();
  h_signDATA->Draw();
  h_numbkgDATA->Draw("same");
  leg->Draw();
  lar->Draw();
  c->Print("realAlphaRatio.pdf)");

}
