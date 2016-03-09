#include "../setNCUStyle.h"
#include "readHists.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooPlot.h"

R__LOAD_LIBRARY(PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(PDFs/PdfDiagonalizer_cc.so)

using namespace RooFit;

void roofitCheck(std::string path){

  setNCUStyle();
  gStyle->SetTitleSize(0.04,"XYZ");
  gStyle->SetLabelSize(0.03,"XYZ");

  readHist dy100(Form("%s/DYjets/DYJetsToLL_M-50_HT-100to200_13TeV_pseudoTest.root",path.data()));
  readHist dy200(Form("%s/DYjets/DYJetsToLL_M-50_HT-200to400_13TeV_pseudoTest.root",path.data()));
  readHist dy400(Form("%s/DYjets/DYJetsToLL_M-50_HT-400to600_13TeV_pseudoTest.root",path.data()));
  readHist dy600(Form("%s/DYjets/DYJetsToLL_M-50_HT-600toInf_13TeV_pseudoTest.root",path.data()));
 
  TH1D* h_prmass = (TH1D*)(dy100.getHist("corrPRmassAll"))->Clone("h_prmass");

  h_prmass->Reset();
  h_prmass->Add(dy100.getHist("corrPRmassAll"));
  h_prmass->Add(dy200.getHist("corrPRmassAll"));
  h_prmass->Add(dy400.getHist("corrPRmassAll"));
  h_prmass->Add(dy600.getHist("corrPRmassAll"));

  TH1D* h_prmass_hollow = (TH1D*)(dy100.getHist("corrPRmass"))->Clone("h_prmass_hollow");

  h_prmass_hollow->Reset();
  h_prmass_hollow->Add(dy100.getHist("corrPRmass"));
  h_prmass_hollow->Add(dy200.getHist("corrPRmass"));
  h_prmass_hollow->Add(dy400.getHist("corrPRmass"));
  h_prmass_hollow->Add(dy600.getHist("corrPRmass"));

  // Declare observable mass

  RooRealVar mass("mass", "M_{jet}", 0.0, 240.0, "GeV");

  mass.setRange("lowerSB", 40.0, 65.0);
  mass.setRange("upperSB", 145.0, 240);
  mass.setRange("allRange", 40.0, 240);

  // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'mass'

  RooDataHist prmass("prmass", "prmass", mass, Import(*h_prmass));
  RooDataHist prmass_hollow("prmass_hollow", "prmass_hollow", mass, Import(*h_prmass_hollow));

  // Plot and fit a RooDataHist

  RooPlot* prmassframe;

  prmassframe = mass.frame();
  prmass.plotOn(prmassframe);

  // Error Function * Exponential

  RooRealVar slope ("slope",  "slope of the exp", -0.03, -100.0, 0.0);
  RooRealVar offset("offset", "offset of the erf", 39,  0.0, 5000.0);
  RooRealVar width ("width",  "width of the erf",  30,  0.0, 5000.0);
  RooRealVar yieldSideband("yieldSideband", "Zjets normalization in sideband region", 10., 1., 1.e4);
  
  yieldSideband.setVal((Int_t)h_prmass->Integral());
  yieldSideband.setError(TMath::Sqrt((Int_t)(h_prmass->Integral())));
  yieldSideband.setConstant(kTRUE);

  RooErfExpPdf mj_pdf("mj_pdf","fiting mj spectrum",mass,slope,offset,width);
  RooExtendPdf mj_pdf_ext("mj_pdf_ext","extended pdf", mj_pdf, yieldSideband);

  RooFitResult* result = mj_pdf_ext.fitTo(prmass, SumW2Error(kFALSE), Extended(kTRUE), Minimizer("Minuit"), Range("allRange"), Save());

  mj_pdf_ext.plotOn(prmassframe, VisualizeError(*result,1), FillColor(kYellow));
  mj_pdf_ext.plotOn(prmassframe);

  RooRealVar slope_sb ("slope_sb",  "slope of the exp",  slope.getVal(), -1.0, 0.0);
  RooRealVar offset_sb("offset_sb", "offset of the erf", offset.getVal(), 10.0, 100.0);
  RooRealVar width_sb ("width_sb",  "width of the erf",  width.getVal(),  10.0, 100.0);
  RooRealVar yieldSideband_sb("yieldSideband_sb", "Zjets normalization in sideband region", 10.0, 1.0, 1.e4);

  yieldSideband_sb.setVal((Int_t)(0.5*h_prmass_hollow->Integral()));
  yieldSideband_sb.setError(TMath::Sqrt((Int_t)(0.5*h_prmass_hollow->Integral())));
  yieldSideband_sb.setConstant(kTRUE);

  RooErfExpPdf mj_pdf_sb("mj_pdf_sb", "fiting mj spectrum", mass, slope_sb, offset_sb, width_sb);
  RooExtendPdf mj_pdf_ext_sb("mj_pdf_ext_sb", "extended pdf", mj_pdf_sb, yieldSideband_sb);

  //mj_pdf_ext_sb.fitTo(prmass, SumW2Error(kFALSE), Extended(kTRUE), Minimizer("Minuit"), Range("lowerSB,upperSB"));
  //mj_pdf_ext_sb.plotOn(prmassframe, LineStyle(kDashed), LineColor(kRed));
  prmass.plotOn(prmassframe);

  TCanvas* c = new TCanvas("c","",0,0,1000,800);

  c->cd();
  prmassframe->Draw();
  c->Print("roofitCheck.pdf");

}
