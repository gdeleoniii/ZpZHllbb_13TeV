#include "../setNCUStyle.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooPlot.h"

R__LOAD_LIBRARY(PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(PDFs/PdfDiagonalizer_cc.so)

using namespace RooFit;

const Double_t xmin  = 500;
const Double_t xmax  = 5000;
const Int_t    nBins = (xmax-xmin)/100;

Double_t dataLumi  = 3000; //pb-1
Double_t xSecDY100 = 147.4*1.23;
Double_t xSecDY200 = 40.99*1.23;
Double_t xSecDY400 = 5.678*1.23;
Double_t xSecDY600 = 2.198*1.23;

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

/// Main function start ///

void roofitCheck(std::string outputFolder){

  setNCUStyle();
  gStyle->SetMarkerSize(0);
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

  // Declare prefer histogram and add them together

  TH1D* h_corrPRmass    = addSamples(infiles,"corrPRmass_pDA",f_DY100,f_DY200,f_DY400,f_DY600);
  TH1D* h_corrPRmassAll = addSamples(infiles,"corrPRmassAll_pDA",f_DY100,f_DY200,f_DY400,f_DY600);

  // Make the statistics error more like data

  for( Int_t i = 1; i <= h_corrPRmass->GetNbinsX(); i++ ){

    h_corrPRmass   ->SetBinError(i,TMath::Sqrt(h_corrPRmass   ->GetBinContent(i)));
    h_corrPRmassAll->SetBinError(i,TMath::Sqrt(h_corrPRmassAll->GetBinContent(i)));

  }

  // Declare observable mass

  RooRealVar mass("mass", "M_{jet}", 0.0, 240.0, "GeV");

  mass.setRange("lowerSB", 40.0, 65.0);
  mass.setRange("upperSB", 145.0, 240);
  mass.setRange("allRange", 40.0, 240);

  // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'mass'

  RooDataHist corrPRmassAll("corrPRmassAll", "corrPRmassAll", mass, Import(*h_corrPRmassAll));
  RooDataHist corrPRmass("corrPRmass", "corrPRmass", mass, Import(*h_corrPRmass));

  // Plot and fit a RooDataHist

  RooPlot* corrPRmassAllframe;

  corrPRmassAllframe = mass.frame();
  corrPRmassAll.plotOn(corrPRmassAllframe);

  // Error Function * Exponential

  RooRealVar slope ("slope",  "slope of the exp", -0.03, -100.0, 0.0);
  RooRealVar offset("offset", "offset of the erf", 39,  0.0, 5000.0);
  RooRealVar width ("width",  "width of the erf",  30,  0.0, 5000.0);
  RooRealVar yieldSideband("yieldSideband", "Zjets normalization in sideband region", 10., 1., 1.e4);
  
  yieldSideband.setVal( (Int_t)h_corrPRmassAll->Integral() );
  yieldSideband.setError( TMath::Sqrt((Int_t)(h_corrPRmassAll->Integral())) );
  yieldSideband.setConstant(kTRUE);

  RooErfExpPdf mj_pdf("mj_pdf","fiting mj spectrum",mass,slope,offset,width);
  RooExtendPdf mj_pdf_ext("mj_pdf_ext","extended pdf", mj_pdf, yieldSideband);

  RooFitResult* result = mj_pdf_ext.fitTo(corrPRmassAll, SumW2Error(kFALSE), Extended(kTRUE), Minimizer("Minuit"), Range("allRange"), Save());

  mj_pdf_ext.plotOn(corrPRmassAllframe, VisualizeError(*result,1), FillColor(kYellow));
  mj_pdf_ext.plotOn(corrPRmassAllframe);

  RooRealVar slope_sb ("slope_sb",  "slope of the exp",  slope.getVal(), -1.0, 0.0);
  RooRealVar offset_sb("offset_sb", "offset of the erf", offset.getVal(), 10.0, 100.0);
  RooRealVar width_sb ("width_sb",  "width of the erf",  width.getVal(),  10.0, 100.0);
  RooRealVar yieldSideband_sb("yieldSideband_sb", "Zjets normalization in sideband region", 10.0, 1.0, 1.e4);

  yieldSideband_sb.setVal( (Int_t)(0.5*h_corrPRmass->Integral()) );
  yieldSideband_sb.setError( TMath::Sqrt((Int_t)(0.5*h_corrPRmass->Integral())) );
  yieldSideband_sb.setConstant(kTRUE);

  RooErfExpPdf mj_pdf_sb("mj_pdf_sb", "fiting mj spectrum", mass, slope_sb, offset_sb, width_sb);
  RooExtendPdf mj_pdf_ext_sb("mj_pdf_ext_sb", "extended pdf", mj_pdf_sb, yieldSideband_sb);

  mj_pdf_ext_sb.fitTo(corrPRmassAll, SumW2Error(kFALSE), Extended(kTRUE), Minimizer("Minuit"), Range("lowerSB,upperSB"));
  mj_pdf_ext_sb.plotOn(corrPRmassAllframe, LineStyle(kDashed), LineColor(kRed));
  corrPRmassAll.plotOn(corrPRmassAllframe);

  // Output results

  TLatex* lar = new TLatex();

  lar->SetTextSize(0.035);
  lar->SetLineWidth(5);

  TCanvas* c = new TCanvas("c","",0,0,1000,800);

  c->cd();
  corrPRmassAllframe->Draw();
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2015");
  lar->DrawLatexNDC(0.65, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  c->Print("roofitCheck.pdf");

}
