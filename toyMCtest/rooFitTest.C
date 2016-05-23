R__LOAD_LIBRARY(PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(PDFs/PdfDiagonalizer_cc.so)
using namespace RooFit;

void rooFitTest(string channel, string catcut, bool pullTest=true){

  // Suppress all the INFO message

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  // Input files and sum all backgrounds

  TChain* tree = new TChain("tree");

  if( channel != "ele" && channel != "mu" ) return;
  
  tree->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-100to200_13TeV_toyMC.root", channel.data()));
  tree->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-200to400_13TeV_toyMC.root", channel.data()));
  tree->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-400to600_13TeV_toyMC.root", channel.data()));
  tree->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-600toInf_13TeV_toyMC.root", channel.data()));

  // Define all the variables from the trees

  RooRealVar cat ("cat", "", 0, 2);
  RooRealVar mJet("prmass", "M_{jet}", 30.,  300., "GeV");
  RooRealVar mZH ("mllbb",   "M_{ZH}",  0., 2000., "GeV");
  RooRealVar evWeight("evweight", "", 0, 1.e3);

  // Set the range in jet mass

  mJet.setRange("allRange", 30., 300.);
  mJet.setRange("lowSB",    30.,  65.);
  mJet.setRange("highSB",  135., 300.);
  mJet.setRange("signal",  105., 135.);

  RooBinning binsmJet(54, 30, 300);

  RooArgSet variables(cat, mJet, mZH, evWeight);

  TCut catCut = Form("cat==%s", catcut.c_str());
  TCut sbCut  = "prmass>30 && !(prmass>65 && prmass<135) && prmass<300";
  TCut sigCut = "prmass>105 && prmass<135";

  /*******************************************/
  /*                 ALL RANGE               */
  /*******************************************/

  // Create a dataset from a tree -> to process an unbinned likelihood fitting

  RooDataSet dataSet("dataSet", 
		     "dataSet",
		     variables,
		     Cut(catCut),
		     WeightVar(evWeight), 
		     Import(*tree));

  // Define the variables for pdf

  RooRealVar constant("constant",  "slope of the exp", -0.02,  -1.,   0.);
  RooRealVar offset  ("offset",   "offset of the erf",   30., -50., 200.);
  RooRealVar width   ("width",     "width of the erf",  100.,   0., 200.);

  // Set the parameter fixed when fitting

  offset.setConstant(true);

  // Define the pdf and fitting

  RooErfExpPdf model("model", "Error function for Z+jets mass", mJet, constant, offset, width);
  
  RooFitResult* mJet_result = model.fitTo(dataSet, 
					  SumW2Error(true), 
					  Range("allRange"),
					  Strategy(2),
					  Minimizer("Minuit2"), 
					  Save(1));

  // Plot results on a frame

  RooPlot* mJetFrame = mJet.frame();

  dataSet.plotOn(mJetFrame,
		 Binning(binsmJet));
 
  model.plotOn(mJetFrame, 
	       Normalization(dataSet.sumEntries(),RooAbsReal::NumEvent),
	       VisualizeError(*mJet_result),
	       FillColor(kYellow));

  dataSet.plotOn(mJetFrame,
		 Binning(binsmJet));

  model.plotOn(mJetFrame, 
	       Normalization(dataSet.sumEntries(),RooAbsReal::NumEvent));

  /*******************************************/
  /*                SIDE BAND                */
  /*******************************************/

  RooDataSet dataSetSB("dataSetSB",
		       "dataSetSB",
		       variables,
		       Cut(catCut && sbCut), 
		       WeightVar(evWeight),
		       Import(*tree));


  RooRealVar constantSB("constantSB",  "slope of the exp", -0.02,  -1.,   0.);
  RooRealVar offsetSB  ("offsetSB",   "offset of the erf", offset.getVal(), -50., 200.);
  RooRealVar widthSB   ("widthSB",     "width of the erf", 70.,   0., 200.);

  offsetSB.setConstant(true);

  RooErfExpPdf modelSB("modelSB", "Error function for Z+jets mass", mJet, constantSB, offsetSB, widthSB);

  RooFitResult* mJetSB_result = modelSB.fitTo(dataSetSB,
					      SumW2Error(true),
					      Range("lowSB,highSB"),
					      Strategy(2),
					      Minimizer("Minuit2"),
					      Save(1));

  RooPlot* mJetFrameSB = mJet.frame();

  dataSetSB.plotOn(mJetFrameSB,
		   Binning(binsmJet));

  modelSB.plotOn(mJetFrameSB, 
		 Normalization(dataSetSB.sumEntries(),RooAbsReal::NumEvent),
		 Range("allRange"),
		 VisualizeError(*mJetSB_result),
		 FillColor(kYellow));

  dataSetSB.plotOn(mJetFrameSB,
		   Binning(binsmJet));

  modelSB.plotOn(mJetFrameSB, 
		 Normalization(dataSetSB.sumEntries(),RooAbsReal::NumEvent),
		 Range("allRange"));

  model.plotOn(mJetFrameSB,
	       Normalization(dataSet.sumEntries(),RooAbsReal::NumEvent),
	       Range("allRange"),
	       LineStyle(7),
	       LineColor(kRed));

  /*******************************************/
  /*            BIAS AND PULL                */
  /*******************************************/

  // Produce n toyMCs to study fit bias and pull  
  // Properties of pull: mean is 0 if there is no bias; width is 1 if error is correct
  // Fit is converge: the fit really finds a set of parameter values that minimizes -log likelihood instead of finding a local minima

  TH1D* h_bias = new TH1D("h_bias", "", 16, -2, 2);
  TH1D* h_pull = new TH1D("h_pull", "", 40, -5, 5);

  RooMsgService::instance().setSilentMode(true);

  for( int ntoy = 1999; ntoy >= 0; --ntoy ){

    if( !pullTest ) break;

    RooArgSet mjet(mJet);

    RooDataSet* setToyMC = modelSB.generate(mjet, dataSet.sumEntries());
    RooDataSet  thisToyMC("thisToyMC", "thisToyMC", mjet, Cut(sbCut), Import(*setToyMC));

    RooRealVar constant_toyMC("constant_toyMC",  "slope of the exp", -0.02,  -1.,   0.);
    RooRealVar offset_toyMC  ("offset_toyMC",   "offset of the erf", offset.getVal(), -50., 200.);
    RooRealVar width_toyMC   ("width_toyMC",     "width of the erf", 70.,   0., 200.);

    offset_toyMC.setConstant(true);

    RooErfExpPdf model_toyMC("model_toyMC", "Error function for Z+jets mass", mJet, constant_toyMC, offset_toyMC, width_toyMC);
    
    RooFitResult* toyMC_result = model_toyMC.fitTo(thisToyMC,
						   SumW2Error(true),
						   Range("lowSB,highSB"),
						   Strategy(2),
						   Minimizer("Minuit2"),
						   Save(1));
    if( toyMC_result->status() != 0 ) continue;
   
    double nsbreal = setToyMC->sumEntries(sbCut)/setToyMC->sumEntries();

    RooRealVar nSBReal("nSBReal", "", nsbreal, 0., 1.);
   
    RooAbsReal* nSIGFit = model_toyMC.createIntegral(mjet, NormSet(mjet), Range("signal"));
    RooAbsReal* nSBFit  = model_toyMC.createIntegral(mjet, NormSet(mjet), Range("lowSB,highSB"));

    RooFormulaVar formula("formula", "ev in signal region of toyMC", "@0*@1/@2", RooArgList(nSBReal, *nSIGFit, *nSBFit));
    
    double nSigFit  = nSIGFit->getVal();
    double nSigReal = setToyMC->sumEntries(sigCut)/setToyMC->sumEntries();
    double fitUnc   = formula.getPropagatedError(*toyMC_result);

    h_bias->Fill((nSigFit - nSigReal)/nSigReal);
    h_pull->Fill((nSigFit - nSigReal)/fitUnc);

  } // End of ntoy loop

  RooMsgService::instance().setSilentMode(false);

  RooRealVar bias("bias", "bias", -2, 2);
  RooRealVar pull("pull", "pull", -5, 5);

  RooDataHist hbias("hbias", "", bias, Import(*h_bias));
  RooDataHist hpull("hpull", "", pull, Import(*h_pull));

  RooRealVar  m("m",  "mean", 0, -10, 10);
  RooRealVar  s("s", "sigma", 3, 0.1, 10);
  RooGaussian g("g", "gauss", pull, m, s);

  g.fitTo(hpull);

  RooPlot* biasFrame = bias.frame();
  hbias.plotOn(biasFrame);

  RooPlot* pullFrame = pull.frame();
  hpull.plotOn(pullFrame);
  g.plotOn(pullFrame);
  g.paramOn(pullFrame);

  /*******************************************/
  /*                 OUTPUT                  */
  /*******************************************/

  TCanvas* c = new TCanvas("c","",0,0,1000,800);

  c->cd();
  mJetFrame->Draw();
  c->Print(Form("rooFit_toyMC_%s_cat%s.pdf(", channel.data(), catcut.data()));

  c->cd();
  mJetFrameSB->Draw();
  c->Print(Form("rooFit_toyMC_%s_cat%s.pdf",  channel.data(), catcut.data()));
  
  c->cd();
  biasFrame->Draw();
  c->Print(Form("rooFit_toyMC_%s_cat%s.pdf",  channel.data(), catcut.data()));

  c->cd();
  pullFrame->Draw();
  c->Print(Form("rooFit_toyMC_%s_cat%s.pdf)", channel.data(), catcut.data()));

}
