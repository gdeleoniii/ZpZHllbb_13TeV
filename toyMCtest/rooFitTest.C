R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/PdfDiagonalizer_cc.so)
using namespace RooFit;

void rooFitTest(string channel, string catcut, bool pullTest=true){

  // Suppress all the INFO message

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooMsgService::instance().setSilentMode(true);

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
  RooRealVar evWeight("evweight", "", 0., 1.e3);

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

  // Define the pdf and fitting

  RooRealVar lamda("lamda", "lamda", -0.02, -5., 5.);
  RooExponential model("model", "Exponential function for Z+jets mass", mJet, lamda);

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

  TLegend* leg1 = new TLegend(0.60,0.67,0.85,0.80);

  leg1->AddEntry(mJetFrame->findObject(mJetFrame->nameOf(0)), "MC with statistical errors", "lep");
  leg1->AddEntry(mJetFrame->findObject(mJetFrame->nameOf(3)), "Fit curve with errors", "l");
  leg1->Draw();
  
  mJetFrame->addObject(leg1);
  mJetFrame->SetTitle("");

  /*******************************************/
  /*                SIDE BAND                */
  /*******************************************/

  RooDataSet dataSetSB("dataSetSB",
		       "dataSetSB",
		       variables,
		       Cut(catCut && sbCut), 
		       WeightVar(evWeight),
		       Import(*tree));

  RooRealVar lamdaSB("lamdaSB", "lamda", -0.02, -5., 5.);
  RooExponential modelSB("modelSB", "Exponential function for Z+jets mass", mJet, lamdaSB);

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

  TLegend* leg2 = new TLegend(0.60,0.62,0.85,0.80);
 
  leg2->AddEntry(mJetFrameSB->findObject(mJetFrameSB->nameOf(0)), "MC with statistical errors", "lep");
  leg2->AddEntry(mJetFrameSB->findObject(mJetFrameSB->nameOf(3)), "Fit curve with errors", "l");
  leg2->AddEntry(mJetFrameSB->findObject(mJetFrameSB->nameOf(4)), "Fit curve of all range", "l");
  leg2->Draw();
  
  mJetFrameSB->addObject(leg2);
  mJetFrameSB->SetTitle("");

  /*******************************************/
  /*            BIAS AND PULL                */
  /*******************************************/

  // Produce n toyMCs to study fit bias and pull  
  // Properties of pull: mean is 0 if there is no bias; width is 1 if error is correct
  // Fit is converge: the fit really finds a set of parameter values that minimizes -log likelihood instead of finding a local minima

  TH1F* h_bias = new TH1F("h_bias", "", 16, -2, 2);
  TH1F* h_pull = new TH1F("h_pull", "", 40, -5, 5);

  RooMsgService::instance().setSilentMode(true);

  for( int ntoy = 1000; ntoy > 0; --ntoy ){

    if( !pullTest ) break;

    RooArgSet mjet(mJet);

    RooDataSet* setToyMC = modelSB.generate(mjet, dataSet.sumEntries());
    RooDataSet  thisToyMC("thisToyMC", "thisToyMC", mjet, Cut(sbCut), Import(*setToyMC));

    RooRealVar lamda_toyMC("lamda_toyMC", "lamda", -0.02, -5., 5.);
    RooExponential model_toyMC("model_toyMC", "Exponential function for Z+jets mass", mJet, lamda_toyMC);

    RooFitResult* toyMC_result = model_toyMC.fitTo(thisToyMC,
						   SumW2Error(true),
						   Range("lowSB,highSB"),
						   Strategy(2),
						   Minimizer("Minuit2"),
						   Save(1));
    if( toyMC_result->status() != 0 ) continue;
   
    float nsbreal = setToyMC->sumEntries(sbCut)/setToyMC->sumEntries();

    RooRealVar nSBReal("nSBReal", "", nsbreal, 0., 1.);
   
    RooAbsReal* nSIGFit = model_toyMC.createIntegral(mjet, NormSet(mjet), Range("signal"));
    RooAbsReal* nSBFit  = model_toyMC.createIntegral(mjet, NormSet(mjet), Range("lowSB,highSB"));

    RooFormulaVar formula("formula", "ev in signal region of toyMC", "@0*@1/@2", RooArgList(nSBReal, *nSIGFit, *nSBFit));
    
    float nSigFit  = nSIGFit->getVal();
    float nSigReal = setToyMC->sumEntries(sigCut)/setToyMC->sumEntries();
    float fitUnc   = formula.getPropagatedError(*toyMC_result);

    h_bias->Fill((nSigFit - nSigReal)/nSigReal);
    h_pull->Fill((nSigFit - nSigReal)/fitUnc);

  } // End of ntoy loop

  RooMsgService::instance().setSilentMode(false);

  RooRealVar bias("bias", "Bias", -2, 2);
  RooRealVar pull("pull", "Pull", -5, 5);

  RooDataHist hbias("hbias", "", bias, Import(*h_bias));
  RooDataHist hpull("hpull", "", pull, Import(*h_pull));

  RooRealVar  mean("mean", "mean", 0, -10, 10);
  RooRealVar  sigma("sigma", "sigma", 3, 0.1, 10);
  RooGaussian gb("gb", "gauss", bias, mean, sigma);
  RooGaussian gp("gp", "gauss", pull, mean, sigma);

  gb.fitTo(hbias);

  RooPlot* biasFrame = bias.frame();
  hbias.plotOn(biasFrame);
  gb.plotOn(biasFrame);
  gb.paramOn(biasFrame,Layout(0.675,0.9,0.8));
  biasFrame->getAttText()->SetTextSize(0.025);
  biasFrame->SetTitle("");

  gp.fitTo(hpull);

  RooPlot* pullFrame = pull.frame();
  hpull.plotOn(pullFrame);
  gp.plotOn(pullFrame);
  gp.paramOn(pullFrame,Layout(0.675,0.9,0.8));
  pullFrame->getAttText()->SetTextSize(0.025);
  pullFrame->SetTitle("");

  /*******************************************/
  /*                 OUTPUT                  */
  /*******************************************/

  TLatex* lar = new TLatex();

  lar->SetTextSize(0.035);
  lar->SetLineWidth(5);
  
  TCanvas* c = new TCanvas("c","",0,0,1000,800);

  c->cd();
  mJetFrame->Draw();
  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.75, 0.85, Form("%s  %s btag", channel.data(), catcut.data()));
  c->Print(Form("rooFit_toyMC_%s_cat%s.pdf(", channel.data(), catcut.data()));

  c->cd();
  mJetFrameSB->Draw();

  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.75, 0.85, Form("%s  %s btag", channel.data(), catcut.data()));
  c->Print(Form("rooFit_toyMC_%s_cat%s.pdf",  channel.data(), catcut.data()));
  
  c->cd();
  biasFrame->Draw();
  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.75, 0.85, Form("%s  %s btag", channel.data(), catcut.data()));
  c->Print(Form("rooFit_toyMC_%s_cat%s.pdf",  channel.data(), catcut.data()));

  c->cd();
  pullFrame->Draw();
  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.75, 0.85, Form("%s  %s btag", channel.data(), catcut.data()));
  c->Print(Form("rooFit_toyMC_%s_cat%s.pdf)", channel.data(), catcut.data()));

}
