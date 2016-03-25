#ifndef __CINT__
#endif
R__LOAD_LIBRARY(PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(PDFs/PdfDiagonalizer_cc.so)
using namespace RooFit;

void rooFitTest(std::string path){

  // Suppress all the INFO message

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooMsgService::instance().setSilentMode(true);

  // Input files and sum all backgrounds

  TChain* tree = new TChain("tree");

  tree->Add(Form("%s/DYjets/DYJetsToLL_M-50_HT-100to200_13TeV_toyMCnew.root", path.data()));
  tree->Add(Form("%s/DYjets/DYJetsToLL_M-50_HT-200to400_13TeV_toyMCnew.root", path.data()));
  tree->Add(Form("%s/DYjets/DYJetsToLL_M-50_HT-400to600_13TeV_toyMCnew.root", path.data()));
  tree->Add(Form("%s/DYjets/DYJetsToLL_M-50_HT-600toInf_13TeV_toyMCnew.root", path.data()));

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

  RooArgSet variables(cat, mJet, mZH, evWeight);

  TCut catCut = "cat==1";
  TCut sbCut  = "prmass>40 && !(prmass>65 && prmass<135) && prmass<240";
  TCut sigCut = "prmass>105 && prmass<135";


  /*******************************************/
  /*                 ALL RANGE               */
  /*******************************************/

  RooDataSet setDYjets("setDYjets", 
		       "setDYjets",
		       variables,
		       Cut(catCut),
		       WeightVar(evWeight), 
		       Import(*tree));

  RooRealVar constant("constant",  "slope of the exp",   -0.02,  -1.,   0.);
  RooRealVar offset  ("offset",   "offset of the erf",     30., -50., 200.);
  RooRealVar width   ("width",     "width of the erf",    100.,   0., 200.);

  offset.setConstant(true);

  RooErfExpPdf model("model", "Error function for Z+jets mass", mJet, constant, offset, width);

  RooFitResult* mJet_result = model.fitTo(setDYjets, 
					  SumW2Error(true), 
					  Range("allRange"),
					  Strategy(2),
					  Minimizer("Minuit2"), 
					  Save(1));


  RooPlot* mJetFrame = mJet.frame();

  setDYjets.plotOn(mJetFrame, Binning(binsmJet));
  model.plotOn(mJetFrame, VisualizeError(*mJet_result,1), FillColor(kYellow));
  setDYjets.plotOn(mJetFrame, Binning(binsmJet));
  model.plotOn(mJetFrame);

  // Produce n toyMCs to study fit bias and pull                                                                                                                                                        

  RooMCStudy toyMC(model, model, mJet);
  toyMC.generateAndFit(1000, setDYjets.sumEntries(), "", "");
  RooPlot* pullconstFrame = toyMC.plotPull(constant, -5., 5., 25, true);
  RooPlot* pullwidthFrame = toyMC.plotPull(width, -5., 5., 25, true);
  

  /*******************************************/
  /*                SIDE BAND                */
  /*******************************************/

  RooDataSet setDYjetsSB("setDYjetsSB",
			 "setDYjetsSB",
			 variables,
			 Cut(catCut && sbCut), 
			 WeightVar(evWeight),
			 Import(*tree));


  RooRealVar constantSB("constantSB",  "slope of the exp", -0.02,  -1.,   0.);
  RooRealVar offsetSB  ("offsetSB",   "offset of the erf",   30., -50., 200.);
  RooRealVar widthSB   ("widthSB",     "width of the erf",   70.,   0., 200.);

  offsetSB.setConstant(true);
  width.setConstant(true);

  RooErfExpPdf modelSB("modelSB", "Error function for Z+jets mass", mJet, constant, offset, width);

  RooFitResult* mJetSB_result = modelSB.fitTo(setDYjetsSB,
					      SumW2Error(true),
					      Range("lowSB,highSB"),
					      Strategy(2),
					      Minimizer("Minuit2"),
					      Save(1));

  RooPlot* mJetFrameSB = mJet.frame();

  setDYjetsSB.plotOn(mJetFrameSB, Binning(binsmJet));
  modelSB.plotOn(mJetFrameSB, Range("allRange"), VisualizeError(*mJetSB_result,1), FillColor(kYellow));
  setDYjetsSB.plotOn(mJetFrameSB, Binning(binsmJet));
  modelSB.plotOn(mJetFrameSB, Range("allRange"));


  /*******************************************/
  /*            BIAS AND PULL                */
  /*******************************************/

  // Produce n toyMCs to study fit bias and pull  
  // Properties of pull: mean is 0 if there is no bias; width is 1 if error is correct

  TH1D* h_bias = new TH1D("h_bias", "", 50, -3, 3);
  TH1D* h_pull = new TH1D("h_pull", "", 50, -3, 3);

  for( int ntoy = 0; ntoy < 1000; ntoy++ ){

    RooArgSet mjet(mJet);

    RooDataSet* setToyMC = model.generate(mjet, setDYjets.sumEntries());
    RooDataSet  thisToyMC("thisToyMC", "thisToyMC", mjet, Cut(sbCut), Import(*setToyMC));

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

    RooAbsReal* getnFit = model_toyMC.createIntegral(mjet, Range("signal"));

    double nSigFit  = getnFit->getVal();
    double nSigReal = setDYjets.sumEntries(sigCut);

    RooFormulaVar formula("formula", "formula", mjet);
    double fitUnc = formula.getPropagatedError(toyMC_result);

    h_bias->Fill((nSigFit - nSigReal)/nSigReal);
    h_pull->Fill((nSigFit - nSigReal)/fitUnc);

  } // End of ntoy loop


  /*******************************************/
  /*                 OUTPUT                  */
  /*******************************************/

  TCanvas* c = new TCanvas("c","",0,0,1000,800);
  c->cd();
  mJetFrame->Draw();
  c->Print("roofitCheck.pdf(");

  c->cd();
  mJetFrameSB->Draw();
  c->Print("roofitCheck.pdf");

  c->cd();
  pullconstFrame->Draw();
  c->Print("roofitCheck.pdf");

  c->cd();
  pullwidthFrame->Draw();
  c->Print("roofitCheck.pdf)");
  
}
