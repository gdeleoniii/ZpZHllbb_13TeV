R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/PdfDiagonalizer_cc.so)
using namespace RooFit;

void rooFitData(string channel, string catcut, bool removeMinor=true){

  // Suppress all the INFO message

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooMsgService::instance().setSilentMode(true);

  // Input files and sum all backgrounds

  TChain* treeData  = new TChain("tree");
  TChain* treeZjets = new TChain("tree");

  if( channel == "ele" ){

    treeData->Add(Form("data/SingleElectron-Run2015D-v1_%sMiniTree.root", channel.data()));
    treeData->Add(Form("data/SingleElectron-Run2015D-v4_%sMiniTree.root", channel.data()));

  }

  else if( channel == "mu" ){

    treeData->Add(Form("data/SingleMuon-Run2015D-v1_%sMiniTree.root", channel.data()));
    treeData->Add(Form("data/SingleMuon-Run2015D-v4_%sMiniTree.root", channel.data()));

  }

  else return;

  treeZjets->Add(Form("Zjets/DYJetsToLL_M-50_HT-100to200_13TeV_%sMiniTree.root", channel.data()));
  treeZjets->Add(Form("Zjets/DYJetsToLL_M-50_HT-200to400_13TeV_%sMiniTree.root", channel.data()));
  treeZjets->Add(Form("Zjets/DYJetsToLL_M-50_HT-400to600_13TeV_%sMiniTree.root", channel.data()));
  treeZjets->Add(Form("Zjets/DYJetsToLL_M-50_HT-600toInf_13TeV_%sMiniTree.root", channel.data()));

  // To remove minor background contribution in data set (weight is -1)

  if( removeMinor ){

    treeData->Add(Form("minor/WW_TuneCUETP8M1_13TeV_%sMiniTree.root",     channel.data()));
    treeData->Add(Form("minor/WZ_TuneCUETP8M1_13TeV_%sMiniTree.root",     channel.data()));
    treeData->Add(Form("minor/ZZ_TuneCUETP8M1_13TeV_%sMiniTree.root",     channel.data()));
    treeData->Add(Form("minor/TT_TuneCUETP8M1_13TeV_%sMiniTree.root",     channel.data()));
    treeData->Add(Form("minor/ZH_HToBB_ZToLL_M125_13TeV_%sMiniTree.root", channel.data()));

  }

  // Define all the variables from the trees

  RooRealVar cat ("cat", "", 0, 2);
  RooRealVar mJet("prmass", "M_{jet}", 30., 300., "GeV");
  RooRealVar mZH ("mllbb", "M_{ZH}", 800., 4000., "GeV");
  RooRealVar evWeight("evweight", "", 0., 1.e3);

  // Set the range in zh mass and in jet mass

  mZH.setRange("fullRange", 800., 4000.);

  mJet.setRange("allRange", 30., 300.);
  mJet.setRange("lowSB",    30.,  65.);
  mJet.setRange("highSB",  135., 300.);
  mJet.setRange("signal",  105., 135.);

  RooBinning binsmZH(32, 800, 4000);
  RooBinning binsmJet(27, 30, 300);

  RooArgSet variables(cat, mJet, mZH, evWeight);

  TCut catCut = Form("cat==%s", catcut.c_str());
  TCut sbCut  = "prmass>30 && !(prmass>65 && prmass<135) && prmass<300";
  TCut sigCut = "prmass>105 && prmass<135";
  TCut specialCut = "(mllbb<1600 || mllbb>1700)";

  // Create a dataset from a tree -> to process an unbinned likelihood fitting

  RooDataSet dataSetData   ("dataSetData",    "dataSetData",    variables, Cut(catCut),           WeightVar(evWeight), Import(*treeData));
  RooDataSet dataSetDataSB ("dataSetDataSB",  "dataSetDataSB",  variables, Cut(catCut && sbCut),  WeightVar(evWeight), Import(*treeData));
  RooDataSet dataSetDataSG ("dataSetDataSG",  "dataSetDataSG",  variables, Cut(catCut && sigCut), WeightVar(evWeight), Import(*treeData));
  RooDataSet dataSetZjets  ("dataSetZjets",   "dataSetZjets",   variables, Cut(catCut),           WeightVar(evWeight), Import(*treeZjets));
  RooDataSet dataSetZjetsSB("dataSetZjetsSB", "dataSetZjetsSB", variables, Cut(catCut && sbCut && specialCut),  WeightVar(evWeight), Import(*treeZjets));  
  RooDataSet dataSetZjetsSG("dataSetZjetsSG", "dataSetZjetsSG", variables, Cut(catCut && sigCut), WeightVar(evWeight), Import(*treeZjets));
  
  // Total events number

  RooRealVar nMcEvents    ("nMcEvents",     "nMcEvents",     0., 1.e10);
  RooRealVar nSBMcEvents  ("nSBMcEvents",   "nSBMcEvents",   0., 1.e10);
  RooRealVar nSGMcEvents  ("nSGMcEvents",   "nSGMcEvents",   0., 1.e10);
  RooRealVar nDataEvents  ("nDataEvents",   "nDataEvents",   0., 1.e10);
  RooRealVar nSBDataEvents("nSBDataEvents", "nSBDataEvents", 0., 1.e10);

  nMcEvents.setVal(dataSetZjets.sumEntries());
  nMcEvents.setConstant(true);
  
  nSBMcEvents.setVal(dataSetZjetsSB.sumEntries());
  nSBMcEvents.setConstant(true);
  
  nSGMcEvents.setVal(dataSetZjetsSG.sumEntries());
  nSGMcEvents.setConstant(true);
  
  nDataEvents.setVal(dataSetData.sumEntries());
  nDataEvents.setConstant(true);

  nSBDataEvents.setVal(dataSetDataSB.sumEntries());
  nSBDataEvents.setConstant(true);

  // Side band jet mass in data

  RooRealVar lamda("lamda", "lamda", -0.025, -0.030, -0.005);

  RooExponential model_mJetSB("model_mJetSB", "model_mJetSB", mJet, lamda);
  RooExtendPdf ext_model_mJetSB("ext_model_mJetSB", "ext_model_mJetSB", model_mJetSB, nSBDataEvents);

  RooFitResult* mJetSB_result = ext_model_mJetSB.fitTo(dataSetDataSB, SumW2Error(true), Extended(true), Range("lowSB,highSB"), Strategy(2), Minimizer("Minuit2"), Save(1));

  fprintf(stdout, "lamda value is %g\n", lamda.getVal());

  // Normalize factor to normalize the background in signal region of data

  RooAbsReal* nSIGFit = ext_model_mJetSB.createIntegral(RooArgSet(mJet), NormSet(mJet), Range("signal"));
  RooAbsReal* nSBFit  = ext_model_mJetSB.createIntegral(RooArgSet(mJet), NormSet(mJet), Range("lowSB,highSB"));

  float normFactor = nSBDataEvents.getVal()*(nSIGFit->getVal()/nSBFit->getVal());
 
  fprintf(stdout, "The normalization factor is %g\n", normFactor);
 
  // Alpha ratio part
  
  float bmin, bmax;

  if( channel == "ele" ){
    bmin = (catcut=="1") ?  600. : 1400.;
    bmax = (catcut=="1") ? 1500. : 2600.;
  }

  else if( channel == "mu" ){
    bmin = (catcut=="1") ? 100. : 2100.; 
    bmax = (catcut=="1") ? 900. : 2900.;
  }
    
  RooRealVar p0("p0", "p0", -0.002, -0.005, 0.);
  RooRealVar p1("p1", "p1", (bmin+bmax)*0.5, bmin, bmax);

  RooGenericPdf model_ZHSB("model_ZHSB", "model_ZHSB", "TMath::Exp(@1*@0+@2/@0)", RooArgSet(mZH,p0,p1));
  RooExtendPdf ext_model_ZHSB("ext_model_ZHSB", "ext_model_ZHSB", model_ZHSB, nSBMcEvents);

  RooFitResult* mZHSB_result = ext_model_ZHSB.fitTo(dataSetZjetsSB, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));
  RooAbsReal* nZHSBFit = ext_model_ZHSB.createIntegral(RooArgSet(mZH), NormSet(mZH), Range("fullRange"));

  // Fit ZH mass in signal region

  float dmin, dmax;
   
  if( channel == "ele" ){
    dmin = (catcut=="1") ? 3100. : 0.;
    dmax = (catcut=="1") ? 3900. : 1.;
  }
    
  else if( channel == "mu" ){
    dmin = (catcut=="1") ? 0. : 8.;
    dmax = (catcut=="1") ? 1. : 18.;
  }
    
  RooRealVar p2("p2", "p2", -0.002, -0.005, 0.);
  RooRealVar p3("p3", "p3", (dmin+dmax)*0.5, dmin, dmax);

  RooGenericPdf model_ZHSG("model_ZHSG", "model_ZHSG", "TMath::Exp(@1*@0+@2/@0)", RooArgSet(mZH,p2,p3));
  RooExtendPdf ext_model_ZHSG("ext_model_ZHSG", "ext_model_ZHSG", model_ZHSG, nSGMcEvents);

  RooFitResult* mZHSG_result = ext_model_ZHSG.fitTo(dataSetZjetsSG, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));
  RooAbsReal* nZHSGFit = ext_model_ZHSG.createIntegral(RooArgSet(mZH), NormSet(mZH), Range("fullRange"));

  // Fit ZH mass in side band region (data)

  float pdmin, pdmax;
    
  if( channel == "ele" ){
    pdmin = 0.;
    pdmax = 2.;
  }
    
  else if( channel == "mu" ){
    pdmin = (catcut=="1") ? 5000. : 1300.;
    pdmax = (catcut=="1") ? 5500. : 1500.;
  }

  RooRealVar pd0("pd0", "pd0", -0.002, -0.005, 0.);
  RooRealVar pd1("pd1", "pd1", (pdmin+pdmax)*0.5, pdmin, pdmax);

  RooGenericPdf model_ZH("model_ZH", "model_ZH", "TMath::Exp(@1*@0+@2/@0)", RooArgSet(mZH,pd0,pd1));
  RooExtendPdf  ext_model_ZH("ext_model_ZH", "ext_model_ZH", model_ZH, nDataEvents);

  RooFitResult* mZH_result = ext_model_ZH.fitTo(dataSetDataSB, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));

  fprintf(stdout, "p0=%f\tp1=%f\tp2=%f\tp3=%f\tpd0=%f\tpd1=%f\n\n", p0.getVal(),p1.getVal(),p2.getVal(),p3.getVal(),pd0.getVal(),pd1.getVal());

  // Multiply the model of background in data side band with the model of alpha ratio to the a model of background in data signal region

  float normConst = ((TF1*)ext_model_ZHSB.asTF(mZH, RooArgList(p0,p1)))->Integral(800,4000) / ((TF1*)ext_model_ZHSG.asTF(mZH, RooArgList(p2,p3)))->Integral(800,4000);

  RooGenericPdf model_alpha  ("model_alpha", "model_alpha", Form("%f*TMath::Exp(%f*@0+%f/@0)/TMath::Exp(%f*@0+%f/@0)", normConst, p2.getVal(),p3.getVal(),p0.getVal(),p1.getVal()), RooArgSet(mZH));
  RooProdPdf    model_sigData("model_sigData", "ext_model_ZH*model_alpha", RooArgList(ext_model_ZH,model_alpha));

  // Plot the results to a frame 

  /*-*-*-*-*-*-*/

  RooPlot* mJetFrame = mJet.frame();

  dataSetDataSB   .plotOn(mJetFrame, Binning(binsmJet));
  ext_model_mJetSB.plotOn(mJetFrame, Range("allRange"), VisualizeError(*mJetSB_result), FillColor(kYellow));
  dataSetDataSB   .plotOn(mJetFrame, Binning(binsmJet));
  ext_model_mJetSB.plotOn(mJetFrame, Range("allRange"));

  TLegend* leg = new TLegend(0.60,0.72,0.85,0.85);

  leg->AddEntry(mJetFrame->findObject(mJetFrame->nameOf(2)), "Data side band", "lep");
  leg->AddEntry(mJetFrame->findObject(mJetFrame->nameOf(3)), "Fit curve with errors", "l");
  leg->Draw();

  mJetFrame->addObject(leg);
  mJetFrame->SetMinimum(0);
  mJetFrame->SetTitle("");

  /*-*-*-*-*-*-*/

  RooPlot* mZHFrameMC = mZH.frame();

  dataSetZjetsSB.plotOn(mZHFrameMC, Binning(binsmZH));
  ext_model_ZHSB.plotOn(mZHFrameMC, VisualizeError(*mZHSB_result), FillColor(kYellow));
  dataSetZjetsSB.plotOn(mZHFrameMC, Binning(binsmZH), LineColor(kBlue+3), MarkerColor(kBlue+3));
  ext_model_ZHSB.plotOn(mZHFrameMC, LineStyle(7), LineColor(kBlue));

  dataSetZjetsSG.plotOn(mZHFrameMC, Binning(binsmZH));
  ext_model_ZHSG.plotOn(mZHFrameMC, VisualizeError(*mZHSG_result), FillColor(kYellow));
  dataSetZjetsSG.plotOn(mZHFrameMC, Binning(binsmZH), LineColor(kRed+3), MarkerColor(kRed+3));
  ext_model_ZHSG.plotOn(mZHFrameMC, LineStyle(7), LineColor(kRed));

  TLegend* leg0 = new TLegend(0.60,0.67,0.85,0.85);

  leg0->AddEntry(mZHFrameMC->findObject(mZHFrameMC->nameOf(2)), "MC side band", "lep");
  leg0->AddEntry(mZHFrameMC->findObject(mZHFrameMC->nameOf(6)), "MC signal region", "lep");
  leg0->AddEntry(mZHFrameMC->findObject(mZHFrameMC->nameOf(3)), "Fit curve of side band", "l");
  leg0->AddEntry(mZHFrameMC->findObject(mZHFrameMC->nameOf(7)), "Fit curve of signal region", "l");
  leg0->Draw();

  mZHFrameMC->addObject(leg0);
  mZHFrameMC->SetMinimum(0);
  mZHFrameMC->SetTitle("");

  /*-*-*-*-*-*-*/

  RooPlot* mZHFrame = mZH.frame();

  dataSetDataSB.plotOn(mZHFrame, Binning(binsmZH));
  ext_model_ZH .plotOn(mZHFrame, VisualizeError(*mZH_result), FillColor(kYellow));
  dataSetDataSB.plotOn(mZHFrame, Binning(binsmZH));
  ext_model_ZH .plotOn(mZHFrame, LineStyle(7), LineColor(kBlue));

  TLegend* leg1 = new TLegend(0.60,0.72,0.85,0.85);

  leg1->AddEntry(mZHFrame->findObject(mZHFrame->nameOf(2)), "Data side band", "lep");
  leg1->AddEntry(mZHFrame->findObject(mZHFrame->nameOf(3)), "Fit curve with errors", "l");
  leg1->Draw();
  
  mZHFrame->addObject(leg1);
  mZHFrame->SetMinimum(0);
  mZHFrame->SetTitle("");

  /*-*-*-*-*-*-*/

  RooPlot* predictFrame = mZH.frame();

  dataSetDataSG.plotOn(predictFrame, Binning(binsmZH));
  model_sigData.plotOn(predictFrame, Normalization(normFactor, RooAbsReal::NumEvent), LineStyle(7), LineColor(kRed));

  TLegend* leg2 = new TLegend(0.60,0.72,0.85,0.85);

  leg2->AddEntry(predictFrame->findObject(predictFrame->nameOf(0)), "Data signal region", "lep");
  leg2->AddEntry(predictFrame->findObject(predictFrame->nameOf(1)), "Predicted backgrounds", "l");
  leg2->Draw();
  
  predictFrame->addObject(leg2);
  predictFrame->SetMinimum(0);
  predictFrame->SetTitle("");

  /*-*-*-*-*-*-*/

  TLatex* lar = new TLatex();

  lar->SetTextSize(0.035);
  lar->SetLineWidth(5);

  TCanvas* cv = new TCanvas("cv","",0,0,1000,800);

  cv->cd();
  mZHFrameMC->Draw();
  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.65, 0.60, Form("%s  %s btag", channel.data(), catcut.data()));
  cv->Print(Form("rooFit_forData_%s_cat%s.pdf(", channel.data(), catcut.data()));

  cv->Clear();
  cv->cd();
  mZHFrame->Draw();
  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.65, 0.63, Form("%s  %s btag", channel.data(), catcut.data()));
  cv->Print(Form("rooFit_forData_%s_cat%s.pdf", channel.data(), catcut.data()));

  cv->Clear();
  cv->cd();
  predictFrame->Draw();
  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.65, 0.63, Form("%s  %s btag", channel.data(), catcut.data()));
  cv->Print(Form("rooFit_forData_%s_cat%s.pdf", channel.data(), catcut.data()));
  
  cv->Clear();
  cv->cd();
  mJetFrame->Draw();
  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.65, 0.63, Form("%s  %s btag", channel.data(), catcut.data()));
  cv->Print(Form("rooFit_forData_%s_cat%s.pdf)", channel.data(), catcut.data()));

}
