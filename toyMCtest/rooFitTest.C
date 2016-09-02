R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/PdfDiagonalizer_cc.so)
using namespace RooFit;

void rooFitTest(string channel, string catcut, bool pullTest=true){

  // Suppress all the INFO message

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooMsgService::instance().setSilentMode(true);
  gROOT->ProcessLine("gErrorIgnoreLevel=kWarning;");

  // Input files and sum all backgrounds

  TChain* tree = new TChain("tree");

  tree->Add(Form("Zjets/DYJetsToLL_M-50_HT-100to200_13TeV_%sMiniTree.root", channel.data()));
  tree->Add(Form("Zjets/DYJetsToLL_M-50_HT-200to400_13TeV_%sMiniTree.root", channel.data()));
  tree->Add(Form("Zjets/DYJetsToLL_M-50_HT-400to600_13TeV_%sMiniTree.root", channel.data()));
  tree->Add(Form("Zjets/DYJetsToLL_M-50_HT-600toInf_13TeV_%sMiniTree.root", channel.data()));

  // Define all the variables from the trees

  RooRealVar cat ("cat", "", 0, 2);
  RooRealVar mJet("prmass", "M_{jet}", 30., 300., "GeV");
  RooRealVar evWeight("evweight", "", 0., 1.e3);

  // Set the range in jet mass

  mJet.setRange("allRange", 30., 300.);
  mJet.setRange("lowSB",    30.,  65.);
  mJet.setRange("highSB",  135., 300.);
  mJet.setRange("signal",  105., 135.);

  RooBinning binsmJet(54, 30, 300);

  RooArgSet variables(cat, mJet, evWeight);

  TCut catCut = Form("cat==%s", catcut.data());
  TCut sbCut  = "prmass>30 && !(prmass>65 && prmass<135) && prmass<300";
  TCut sigCut = "prmass>105 && prmass<135";

  // Create a dataset from a tree -> to process an unbinned likelihood fitting

  RooDataSet dataSet  ("dataSet",   "dataSet",   variables, Cut(catCut),          WeightVar(evWeight), Import(*tree));
  RooDataSet dataSetSB("dataSetSB", "dataSetSB", variables, Cut(catCut && sbCut), WeightVar(evWeight), Import(*tree));

  // Total events number

  RooRealVar nMcEvents  ("nMcEvents",   "nMcEvents",   0., 1.e10);
  RooRealVar nSBMcEvents("nSBMcEvents", "nSBMcEvents", 0., 1.e10);

  nMcEvents.setVal(dataSet.sumEntries());
  nSBMcEvents.setVal(dataSetSB.sumEntries());

  // ALL RANGE

  RooRealVar     lamda("lamda", "lamda", -0.025, -0.04, -0.01);
  RooExponential model("model", "Exponential function for Z+jets mass", mJet, lamda);
  RooExtendPdf   ext_model("ext_model", "ext_model", model, nMcEvents);
  RooFitResult*  mJet_result = ext_model.fitTo(dataSet, SumW2Error(true), Extended(true), Range("allRange"), Strategy(2), Minimizer("Minuit2"), Save(1));

  // SIDE BAND

  RooRealVar     lamdaSB("lamdaSB", "lamda", -0.025, -0.04, -0.01);
  RooExponential modelSB("modelSB", "Exponential function for Z+jets mass", mJet, lamdaSB);
  RooExtendPdf   ext_modelSB("ext_modelSB", "ext_modelSB", modelSB, nSBMcEvents);
  RooFitResult*  mJetSB_result = ext_modelSB.fitTo(dataSetSB, SumW2Error(true), Extended(true), Range("lowSB,highSB"), Strategy(2), Minimizer("Minuit2"), Save(1));

  fprintf(stdout, "lamda=%f\tlamdaSB=%f\n", lamda.getVal(), lamdaSB.getVal());

  // Produce n toyMCs to study fit bias and pull  
  // Properties of pull: mean is 0 if there is no bias; width is 1 if error is correct
  // Fit is converge: the fit really finds a set of parameter values that minimizes -log likelihood instead of finding a local minima

  TH1F* h_bias = new TH1F("h_bias", "", 19, -9.5, 9.5);
  TH1F* h_pull = new TH1F("h_pull", "", 19, -9.5, 9.5);

  for( int ntoy = 1000; ntoy > 0; --ntoy ){

    if( !pullTest ) break;

    RooArgSet mjet(mJet);

    RooDataSet* setToyMC = ext_model.generate(mjet);
    RooDataSet  thisToyMC("thisToyMC", "thisToyMC", mjet, Cut(sbCut), Import(*setToyMC));

    RooRealVar  nToyMcEvents("nToyMcEvents", "nToyMcEvents", 0., 1.e10);
    nToyMcEvents.setVal(thisToyMC.sumEntries());

    RooRealVar     lamda_toyMC("lamda_toyMC", "lamda", -0.025, -0.04, -0.01);
    RooExponential model_toyMC("model_toyMC", "Exponential function for Z+jets mass", mJet, lamda_toyMC);
    RooExtendPdf   ext_model_toyMC("ext_model_toyMC", "ext_model_toyMC", model_toyMC, nToyMcEvents);
    RooFitResult*  toyMC_result = ext_model_toyMC.fitTo(thisToyMC, SumW2Error(true), Extended(true), Range("lowSB,highSB"), Strategy(2), Minimizer("Minuit2"), Save(1));

    // fprintf(stdout, "nToy=%i\tlamdaToy=%f\tstatus=%i\n", ntoy, lamda_toyMC.getVal(), toyMC_result->status());

    if( toyMC_result->status() != 0 ) continue;
   
    // calulate normalize factor

    RooAbsReal* nSIGFit = ext_model_toyMC.createIntegral(mjet, NormSet(mjet), Range("signal"));
    RooAbsReal* nSBFit  = ext_model_toyMC.createIntegral(mjet, NormSet(mjet), Range("lowSB,highSB"));

    RooRealVar nSBHist("nSBHist", "nSBHist", 0., 1.e10);

    nSBHist.setVal(setToyMC->sumEntries(sbCut));
    nSBHist.setConstant(true);

    float toyNormFactor = nSBHist.getVal()*(nSIGFit->getVal()/nSBFit->getVal());
    float toyNormHiste  = setToyMC->sumEntries(sigCut);

    RooFormulaVar formula("formula", "events in signal region of toyMC", "@0*@1/@2", RooArgList(nSBHist, *nSIGFit, *nSBFit));

    h_bias->Fill((toyNormFactor - toyNormHiste)/toyNormHiste);
    h_pull->Fill((toyNormFactor - toyNormHiste)/formula.getPropagatedError(*toyMC_result));

  } // End of ntoy loop

  RooRealVar bias("bias", "Bias", -9.5, 9.5);
  RooRealVar pull("pull", "Pull", -9.5, 9.5);

  RooDataHist hbias("hbias", "", bias, Import(*h_bias));
  RooDataHist hpull("hpull", "", pull, Import(*h_pull));

  RooRealVar gbmean("gbmean", "mean", 0.0, -5.0, 5.0);
  RooRealVar gbsigma("gbsigma", "sigma", 1.0, 0.1, 4.0);

  RooRealVar gpmean("gpmean", "mean", 0.0, -5.0, 5.0);
  RooRealVar gpsigma("gpsigma", "sigma", 1.0, 0.1, 4.0);

  RooGaussian gb("gb", "gauss", bias, gbmean, gbsigma);
  RooGaussian gp("gp", "gauss", pull, gpmean, gpsigma);

  gb.fitTo(hbias);
  gp.fitTo(hpull);

  // Another toy MC study using RooMCStudy 

  RooMCStudy* mcstudy = new RooMCStudy(model, mJet, Binned(false), Silence(true), Extended(true), FitOptions(Save(1), PrintEvalErrors(false)));
  
  mcstudy->generateAndFit(1000, dataSet.sumEntries());

  // Plot the results on frame 

  RooPlot* mJetFrame        = mJet.frame();
  RooPlot* mJetSBFrame      = mJet.frame();
  RooPlot* biasFrame        = bias.frame();
  RooPlot* pullFrame        = pull.frame();
  RooPlot* mJetPullFrame    = mJet.frame();
  RooPlot* mJetSBPullFrame  = mJet.frame();
  RooPlot* mcstudyPullFrame = mcstudy->plotPull(lamda, FrameRange(-9.5,9.5), Bins(19), FitGauss(true));

  dataSet  .plotOn(mJetFrame, Binning(binsmJet)); 
  ext_model.plotOn(mJetFrame, VisualizeError(*mJet_result,1,false), FillStyle(3002));
  dataSet  .plotOn(mJetFrame, Binning(binsmJet));
  ext_model.plotOn(mJetFrame);

  mJetPullFrame->addObject(mJetFrame->pullHist(), "P");

  dataSetSB  .plotOn(mJetSBFrame, Binning(binsmJet));
  ext_modelSB.plotOn(mJetSBFrame, Range("allRange"), VisualizeError(*mJetSB_result,1,false), FillStyle(3002));
  dataSetSB  .plotOn(mJetSBFrame, Binning(binsmJet));
  ext_modelSB.plotOn(mJetSBFrame, Range("allRange"));

  mJetSBPullFrame->addObject(mJetSBFrame->pullHist(), "P");

  ext_model.plotOn(mJetSBFrame, Normalization(dataSet.sumEntries(),RooAbsReal::NumEvent), Range("allRange"), LineStyle(7), LineColor(kRed));

  hbias.plotOn(biasFrame);
  gb.plotOn(biasFrame);
  gb.paramOn(biasFrame,Layout(0.65,0.9,0.8));

  hpull.plotOn(pullFrame);
  gp.plotOn(pullFrame);
  gp.paramOn(pullFrame,Layout(0.65,0.9,0.8));

  // Output results

  TLatex lar;

  lar.SetTextSize(0.035);
  lar.SetLineWidth(5);
  
  float up_height = 0.82;
  float dw_height = (1-up_height)*1.445;

  TCanvas c0("c0","",0,0,1000,800);
  TLegend leg0(0.60,0.67,0.85,0.80);
  
  c0.Divide(1,2);

  TPad* c0_up = (TPad*)c0.GetListOfPrimitives()->FindObject("c0_1");
  TPad* c0_dw = (TPad*)c0.GetListOfPrimitives()->FindObject("c0_2"); 

  c0_up->SetPad(0,1-up_height,1,1);
  c0_dw->SetPad(0,0,1,dw_height);
  c0_dw->SetBottomMargin(0.25);
  c0_up->cd()->SetLogy(1);

  leg0.AddEntry(mJetFrame->findObject(mJetFrame->nameOf(0)), "MC with statistical errors", "lep");
  leg0.AddEntry(mJetFrame->findObject(mJetFrame->nameOf(3)), "Fit curve with errors", "l");
  leg0.Draw();

  mJetFrame->addObject(&leg0);
  mJetFrame->SetTitle("");
  mJetFrame->SetMinimum(1e-3);
  mJetFrame->SetMaximum(catcut=="1"?100:10);
  mJetFrame->GetXaxis()->SetTitle("");
  mJetFrame->GetXaxis()->SetLabelOffset(999);
  mJetFrame->Draw();

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.62, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.83, Form("%s, %s btag", channel.data(), catcut.data()));

  c0_up->RedrawAxis();

  c0_dw->cd()->SetLogy(0);

  mJetPullFrame->SetTitle("");
  mJetPullFrame->GetYaxis()->SetTitle("Pulls");
  mJetPullFrame->GetYaxis()->SetTitleOffset(0.25);
  mJetPullFrame->GetXaxis()->SetLabelSize(0.125);
  mJetPullFrame->GetXaxis()->SetTitleSize(0.125);
  mJetPullFrame->GetYaxis()->SetLabelSize(0.125);
  mJetPullFrame->GetYaxis()->SetTitleSize(0.125);
  mJetPullFrame->GetYaxis()->SetNdivisions(505);
  mJetPullFrame->SetMinimum(-4);
  mJetPullFrame->SetMaximum(4);
  mJetPullFrame->Draw();

  c0.Draw();
  c0.Print(Form("rooFit_toyMC_%s_cat%s.pdf(", channel.data(), catcut.data()));

  TCanvas c1("c1","",0,0,1000,800);
  TLegend leg1(0.60,0.62,0.85,0.80);

  c1.Divide(1,2);

  TPad* c1_up = (TPad*)c1.GetListOfPrimitives()->FindObject("c1_1");
  TPad* c1_dw = (TPad*)c1.GetListOfPrimitives()->FindObject("c1_2");

  c1_up->SetPad(0,1-up_height,1,1);
  c1_dw->SetPad(0,0,1,dw_height);
  c1_dw->SetBottomMargin(0.25);
  c1_up->cd()->SetLogy(1);

  leg1.AddEntry(mJetSBFrame->findObject(mJetSBFrame->nameOf(0)), "MC with statistical errors", "lep");
  leg1.AddEntry(mJetSBFrame->findObject(mJetSBFrame->nameOf(3)), "Fit curve with errors", "l");
  leg1.AddEntry(mJetSBFrame->findObject(mJetSBFrame->nameOf(4)), "Fit curve of all range", "l");
  leg1.Draw();

  mJetSBFrame->addObject(&leg1);
  mJetSBFrame->SetTitle("");
  mJetSBFrame->SetMinimum(1e-3);
  mJetSBFrame->SetMaximum(catcut=="1"?100:10);
  mJetSBFrame->GetXaxis()->SetTitle("");
  mJetSBFrame->GetXaxis()->SetLabelOffset(999);
  mJetSBFrame->Draw();

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.62, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.83, Form("%s, %s btag", channel.data(), catcut.data()));

  c1_up->RedrawAxis();
  c1_dw->cd()->SetLogy(0);

  mJetSBPullFrame->SetTitle("");
  mJetSBPullFrame->GetYaxis()->SetTitle("Pulls");
  mJetSBPullFrame->GetYaxis()->SetTitleOffset(0.25);
  mJetSBPullFrame->GetXaxis()->SetLabelSize(0.125);
  mJetSBPullFrame->GetXaxis()->SetTitleSize(0.125);
  mJetSBPullFrame->GetYaxis()->SetLabelSize(0.125);
  mJetSBPullFrame->GetYaxis()->SetTitleSize(0.125);
  mJetSBPullFrame->GetYaxis()->SetNdivisions(505);
  mJetSBPullFrame->SetMinimum(-4);
  mJetSBPullFrame->SetMaximum(4);
  mJetSBPullFrame->Draw();

  c1.Draw();
  c1.Print(Form("rooFit_toyMC_%s_cat%s.pdf",  channel.data(), catcut.data()));
  
  TCanvas c("c","",0,0,1000,800);

  c.cd();
  mcstudyPullFrame->getAttText()->SetTextSize(0.025);
  mcstudyPullFrame->SetTitle("");
  mcstudyPullFrame->Draw();
  c.Print(Form("rooFit_toyMC_%s_cat%s.pdf",  channel.data(), catcut.data()));

  c.Clear();
  c.cd();
  biasFrame->getAttText()->SetTextSize(0.025);
  biasFrame->SetTitle("");  
  biasFrame->Draw();
  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.72, 0.83, Form("%s, %s btag", channel.data(), catcut.data()));
  c.Print(Form("rooFit_toyMC_%s_cat%s.pdf",  channel.data(), catcut.data()));

  c.Clear();
  c.cd();
  pullFrame->getAttText()->SetTextSize(0.025);
  pullFrame->SetTitle("");
  pullFrame->Draw();
  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.72, 0.83, Form("%s, %s btag", channel.data(), catcut.data()));
  c.Print(Form("rooFit_toyMC_%s_cat%s.pdf)", channel.data(), catcut.data()));

}
