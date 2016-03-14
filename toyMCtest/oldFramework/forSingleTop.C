#include <string>
#include <iostream>
#include <TF1.h>
#include <TH1.h>
#include <TGraphAsymmErrors.h>
#include "readHists.h"
#include "fitTest.h"
#include "fitDraw.h"

void forSingleTop(std::string path, std::string pdfName){

  readHist st1(Form("%s/singleTop/ST_s_4f_leptonDecays_13TeV_pseudoTest.root",path.data()));
  readHist st2(Form("%s/singleTop/ST_t_antitop_4f_leptonDecays_13TeV_pseudoTest.root",path.data()));
  readHist st3(Form("%s/singleTop/ST_t_top_4f_leptonDecays_13TeV_pseudoTest.root",path.data()));
  readHist st4(Form("%s/singleTop/ST_tW_antitop_5f_inclusiveDecays_13TeV_pseudoTest.root",path.data()));
  readHist st5(Form("%s/singleTop/ST_tW_top_5f_inclusiveDecays_13TeV_pseudoTest.root",path.data())); 

  TH1D* h_prmass = (TH1D*)(st1.getHist("corrPRmassAll"))->Clone("h_prmass");

  h_prmass->Reset();
  h_prmass->Add(st1.getHist("corrPRmassAll"));
  h_prmass->Add(st2.getHist("corrPRmassAll"));
  h_prmass->Add(st3.getHist("corrPRmassAll"));
  h_prmass->Add(st4.getHist("corrPRmassAll"));
  h_prmass->Add(st5.getHist("corrPRmassAll"));

  TH1D* h_prmass_hollow = (TH1D*)(st1.getHist("corrPRmass"))->Clone("h_prmass_hollow");

  h_prmass_hollow->Reset();
  h_prmass_hollow->Add(st1.getHist("corrPRmass"));
  h_prmass_hollow->Add(st2.getHist("corrPRmass"));
  h_prmass_hollow->Add(st3.getHist("corrPRmass"));
  h_prmass_hollow->Add(st4.getHist("corrPRmass"));
  h_prmass_hollow->Add(st5.getHist("corrPRmass"));  

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
