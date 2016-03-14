#include "../../setNCUStyle.h"
#include "readHists.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooPlot.h"

R__LOAD_LIBRARY(PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(PDFs/PdfDiagonalizer_cc.so)

using namespace RooFit;

void rooFitTest(std::string path){

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

  // Define variables
  RooRealVar jetMass("jetMass", "M_{jet}", 0, 300, "GeV"); 
  jetMass.setRange("lowerSB",   40,  65);
  jetMass.setRange("upperSB",  135, 300);
  jetMass.setRange("allRange",  40, 300);

  RooRealVar constZjet ("constZjet",   "slope of the exp",  -0.020,  -1.,   0.);
  RooRealVar offsetZjet("offsetZjet", "offset of the erf",     30., -50., 200.);
  RooRealVar widthZjet ("widthZjet",   "width of the erf",    100.,   1., 200.);
  
  offsetZjet.setConstant(true);

  RooErfExpPdf modelZjet("modelZjet", "Error function for Z+jets mass", jetMass, constZjet, offsetZjet, widthZjet);


  RooDataHist h_jetMass("h_jetMass", "", jetMass, Import(*h_prmass));

  RooFitResult* Zjet = modelZjet.fitTo(h_jetMass, 
				       SumW2Error(true), 
				       Range("allRange"), 
				       Strategy(2),
				       Minimizer("Minuit2"), 
				       Save(1));

  // Plot and fit a RooDataHist

  RooPlot* jetMassFrame = jetMass.frame();
  h_jetMass.plotOn(jetMassFrame);
  modelZjet.plotOn(jetMassFrame, VisualizeError(*Zjet,1), FillColor(kYellow));
  modelZjet.plotOn(jetMassFrame);





  // Generate some toy data

  RooAbsData* data = modelZjet.generate(mass, 1000);

  modelZjet.fitTo(*data);



  RooPlot* biasFrame;




  TCanvas* c = new TCanvas("c","",0,0,1000,800);
  c->cd();
  jetMassFrame->Draw();
  c->Print("roofitCheck.pdf(");

  c->cd();
  biasFrame->Draw();
  c->Print("roofitCheck.pdf)");

  delete data;

}
