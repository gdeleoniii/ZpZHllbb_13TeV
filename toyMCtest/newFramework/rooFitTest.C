#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
//#include "RooRealVar.h"
//#include "RooPlot.h"
R__LOAD_LIBRARY(PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(PDFs/PdfDiagonalizer_cc.so)
using namespace RooFit;

void rooFitTest(std::string path){

  // Suppress all the INFO message

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  // Input files and sum all backgrounds

  TChain* tree = new TChain("tree");

  tree->Add(Form("%s/DYjets/DYJetsToLL_M-50_HT-100to200_13TeV_toyMCnew.root", path.data()));
  tree->Add(Form("%s/DYjets/DYJetsToLL_M-50_HT-200to400_13TeV_toyMCnew.root", path.data()));
  tree->Add(Form("%s/DYjets/DYJetsToLL_M-50_HT-400to600_13TeV_toyMCnew.root", path.data()));
  tree->Add(Form("%s/DYjets/DYJetsToLL_M-50_HT-600toInf_13TeV_toyMCnew.root", path.data()));




  TH1D* h = new TH1D("h","",40,40,240);
  h->Sumw2();
  int ccat;
  float mjet,ew;
  tree->SetBranchAddress("cat",&ccat);
  tree->SetBranchAddress("prmass",&mjet);
  tree->SetBranchAddress("evweight",&ew);
  for(long i=0; i<tree->GetEntries(); i++ ){
    tree->GetEntry(i);
    if(ccat != 1) continue;
    h->Fill(mjet,ew);
  }




  // Define all the variables from the trees

  RooRealVar cat ("cat", "", 0, 2);
  RooRealVar mJet("prmass", "M_{jet}", 40.,  240., "GeV");
  RooRealVar mZH ("mllbb",   "M_{ZH}",  0., 2000., "GeV");
  RooRealVar evWeight("evweight", "", 0, 1.e3);

  // Set the range in jet mass

  mJet.setRange("allRange", 40., 240.);
  mJet.setRange("lowSB",    30.,  65.);
  mJet.setRange("highSB",  135., 300.);

  RooBinning binsmJet(40, 240);
  binsmJet.addUniform(40, 40, 240);

  // Define the RooArgSet which will include all the variables defined before

  RooArgSet variables(cat, mJet, mZH, evWeight);

  // Create dataset 

  RooDataSet setDYjets("setDYjets", "setDYjets", variables, Cut("cat==1"), WeightVar(evWeight), Import(*tree));

  // Define variables of ErfExp function

  RooRealVar constZjet ("constZjet",  "slope of the exp",   -0.02,  -1.,   0.);
  RooRealVar offsetZjet("offsetZjet", "offset of the erf",    30., -50., 400.);
  RooRealVar widthZjet ("widthZjet",  "width of the erf",    100.,   0., 400.);
  
  offsetZjet.setConstant(true);

  // Define DYjets model

  RooErfExpPdf modelZjet("modelZjet", "Error function for Z+jets mass", mJet, constZjet, offsetZjet, widthZjet);

  // Fit to whole range of jet mass

  RooFitResult* mJet_result = modelZjet.fitTo(setDYjets, 
					      SumW2Error(true), 
					      Range("allRange"), 
					      Strategy(2),
					      Minimizer("Minuit2"), 
					      Save(1));




  // Plot on a frame

  RooPlot* mJetFrame = mJet.frame();

  //modelZjet.plotOn(mJetFrame, VisualizeError(*mJet_result,1,false), FillColor(kYellow));
  modelZjet.plotOn(mJetFrame);
  setDYjets.plotOn(mJetFrame, Binning(binsmJet));


  // Generate some toy data

  //RooAbsData* data = modelZjet.generate(mass, 1000);

  //modelZjet.fitTo(*data);


  TCanvas* c = new TCanvas("c","",0,0,1000,800);
  c->cd();
  mJetFrame->Draw();
  c->Print("roofitCheck.pdf(");

  c->cd();
  h->Draw();
  c->Print("roofitCheck.pdf)");
  //  delete data;

}
