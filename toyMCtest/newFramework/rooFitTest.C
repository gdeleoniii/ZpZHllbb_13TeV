#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
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
  mJet.setRange("lowSB",    40.,  65.);
  mJet.setRange("highSB",  135., 240.);
  mJet.setRange("signal",  105., 135.);

  RooBinning binsmJet(40, 40, 240);

  // Define the RooArgSet which will include all the variables defined before

  RooArgSet variables(cat, mJet, mZH, evWeight);

  // Create dataset 

  TCut catCut = "cat==1";
  TCut sbCut  = "prmass>40 && !(prmass>65 && prmass<135) && prmass<240";

  RooDataSet setDYjets("setDYjets", 
		       "setDYjets",
		       variables,
		       Cut(catCut),
		       WeightVar(evWeight), 
		       Import(*tree));

  RooDataSet setDYjetsSB("setDYjetsSB",
			 "setDYjetsSB",
			 variables,
			 Cut(catCut && sbCut), 
			 WeightVar(evWeight),
			 Import(*tree));

  // Define variables of ErfExp function

  RooRealVar constant("constant",  "slope of the exp",   -0.02,  -1.,   0.);
  RooRealVar offset  ("offset", "offset of the erf",    30., -50., 200.);
  RooRealVar width   ("width",  "width of the erf",    100.,   0., 200.);
  
  offset.setConstant(true);

  // Define DYjets model

  RooErfExpPdf model("model", "Error function for Z+jets mass", mJet, constant, offset, width);

  // Fit to whole range of jet mass

  RooFitResult* mJet_result = model.fitTo(setDYjets, 
					  SumW2Error(true), 
					  Range("allRange"),
					  Strategy(2),
					  Minimizer("Minuit2"), 
					  Save(1));

  RooFitResult* mJetSB_result = model.fitTo(setDYjetsSB,
					    SumW2Error(true),
					    Range("allRange"),
					    Strategy(2),
					    Minimizer("Minuit2"),
					    Save(1));

  // Produce n toyMCs to study fit bias and pull
  

  TH1D* h_bias   = new TH1D("h_bias",   "", 50, -1, 1);
  TH1D* h_upPull = new TH1D("h_upPull", "", 50, -1, 1);
  TH1D* h_dwPull = new TH1D("h_dwPull", "", 50, -1, 1);

  for( int ntoy = 0; ntoy < 1000; ntoy++ ){

    RooDataSet* setToyMC = model.generate(RooArgSet(mJet), setDYjets.sumEntries());
    RooDataSet* thisToyMC("thisToyMC", "thisToyMC", RooArgSet(mJet), Cut(sbCut), Import(*setToyMC));

    RooRealVar constant_toyMC ("constant_toyMC",  "slope of the exp", -0.02, -1., 0.);
    RooRealVar offset_toyMC("offset_toyMC", "offset of the erf", offset.getVal());
    RooRealVar width_toyMC ("width_toyMC",  "width of the erf",  width.getVal());

    RooErfExpPdf model_toyMC("model_toyMC", "Error function for Z+jets mass", mJet, constant_toyMC, offset_toyMC, width_toyMC);
    
    RooFitResult* toyMC_result = model_toyMC.fitTo(thisToyMC,
						   SumW2Error(true),
						   Range("allRange"),
						   Strategy(2),
						   Minimizer("Minuit2"),
						   Save(1));

    double nSigReal = thisToyMC->sumEntries("signal");
    double nSigFit  = model_toyMC->createIntegral("signal");

    // RooFormulaVar formula("formula", "formula", );
    //double upfitUnc = formula.getPropagatedError(toyMC_result);
    //double dwfitUnc = formula.getPropagatedError(toyMC_result);

    h_bias->Fill((nSigFit - nSigReal)/nSigReal);
    //h_upPull->Fill((nSigFit - nSigReal)/upfitUnc);
    //h_dwPull->Fill((nSigFit - nSigReal)/dwfitUnc);

  } // End of ntoy loop

  // Plot on a frame

  RooPlot* mJetFrame = mJet.frame();

  setDYjets.plotOn(mJetFrame, Binning(binsmJet));
  model.plotOn(mJetFrame, VisualizeError(*mJet_result,1), FillColor(kYellow));
  setDYjets.plotOn(mJetFrame, Binning(binsmJet));
  model.plotOn(mJetFrame);

  RooPlot* mJetFrameSB = mJet.frame();

  setDYjetsSB.plotOn(mJetFrameSB, Binning(binsmJet));
  model.plotOn(mJetFrameSB, VisualizeError(*mJetSB_result,1), FillColor(kYellow));
  setDYjetsSB.plotOn(mJetFrameSB, Binning(binsmJet));
  model.plotOn(mJetFrameSB);

  TCanvas* c = new TCanvas("c","",0,0,1000,800);
  c->cd();
  mJetFrame->Draw();
  c->Print("roofitCheck.pdf(");

  c->cd();
  mJetFrameSB->Draw();
  c->Print("roofitCheck.pdf");

  c->cd();
  h->Draw();
  c->Print("roofitCheck.pdf)");

}
