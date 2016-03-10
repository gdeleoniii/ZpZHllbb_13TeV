#include <string>
#include <iostream>
#include <TF1.h>
#include <TH1.h>
#include <TGraphAsymmErrors.h>
#include "readHists.h"
#include "fitTest.h"
#include "fitDraw.h"

void forDibosons(std::string path, std::string pdfName){

  readHist ww(Form("%s/diBosons/WW_TuneCUETP8M1_13TeV_pseudoTest.root",path.data()));
  readHist wz(Form("%s/diBosons/WZ_TuneCUETP8M1_13TeV_pseudoTest.root",path.data()));
  readHist zz(Form("%s/diBosons/ZZ_TuneCUETP8M1_13TeV_pseudoTest.root",path.data()));
 
  TH1D* h_prmass = (TH1D*)(ww.getHist("corrPRmassAll"))->Clone("h_prmass");

  h_prmass->Reset();
  h_prmass->Add(ww.getHist("corrPRmassAll"));
  h_prmass->Add(wz.getHist("corrPRmassAll"));
  h_prmass->Add(zz.getHist("corrPRmassAll"));

  TH1D* h_prmass_hollow = (TH1D*)(ww.getHist("corrPRmass"))->Clone("h_prmass_hollow");

  h_prmass_hollow->Reset();
  h_prmass_hollow->Add(ww.getHist("corrPRmass"));
  h_prmass_hollow->Add(wz.getHist("corrPRmass"));
  h_prmass_hollow->Add(zz.getHist("corrPRmass"));

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
