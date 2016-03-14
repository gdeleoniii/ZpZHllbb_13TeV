#include <string>
#include <iostream>
#include <TF1.h>
#include <TH1.h>
#include <TGraphAsymmErrors.h>
#include "readHists.h"
#include "fitTest.h"
#include "fitDraw.h"

void forDYjets(std::string path, std::string pdfName){

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

  TGraphAsymmErrors* g_errorBands = new TGraphAsymmErrors(); 
  TGraphAsymmErrors* g_errorBands_hollow = new TGraphAsymmErrors();

  TH1D* h_bias   = new TH1D("h_bias",   "", 50, -1, 1);
  TH1D* h_upPull = new TH1D("h_upPull", "", 50, -1, 1);
  TH1D* h_dwPull = new TH1D("h_dwPull", "", 50, -1, 1);

  fitTest(h_prmass, h_prmass_hollow, 
	  &g_errorBands, &g_errorBands_hollow,
	  h_bias, h_upPull, h_dwPull);

  fitDraw(h_prmass, h_prmass_hollow,
	  g_errorBands, g_errorBands_hollow,
	  h_bias, h_upPull, h_dwPull,
	  pdfName);

}
