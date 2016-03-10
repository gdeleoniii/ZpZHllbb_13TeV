#include <string>
#include <iostream>
#include <TF1.h>
#include <TH1.h>
#include <TGraphAsymmErrors.h>
#include "readHists.h"
#include "fitTest.h"
#include "fitDraw.h"

void forTTbar(std::string path, std::string pdfName){

  readHist ttbar(Form("%s/ttbar/TT_TuneCUETP8M1_13TeV_pseudoTest.root",path.data()));
 
  TH1D* h_prmass = (TH1D*)(ttbar.getHist("corrPRmassAll"))->Clone("h_prmass");

  h_prmass->Reset();
  h_prmass->Add(ttbar.getHist("corrPRmassAll"));

  TH1D* h_prmass_hollow = (TH1D*)(ttbar.getHist("corrPRmass"))->Clone("h_prmass_hollow");

  h_prmass_hollow->Reset();
  h_prmass_hollow->Add(ttbar.getHist("corrPRmass"));

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
