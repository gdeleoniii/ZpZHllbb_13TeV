R__LOAD_LIBRARY(PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(PDFs/PdfDiagonalizer_cc.so)
using namespace RooFit;

void rooFitData(string channel, string catcut, bool removeMinor=true){

  // Suppress all the INFO message

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  // Input files and sum all backgrounds

  TChain* treeData  = new TChain("tree");
  TChain* treeZjets = new TChain("tree");

  if( channel == "ele" ){

    treeData->Add(Form("%s/data/SingleElectron-Run2015D-v1_toyMC.root",  channel.data()));
    treeData->Add(Form("%s/data/SingleElectron-Run2015D-v4_toyMC.root", channel.data()));

  }

  else if( channel == "mu" ){

    treeData->Add(Form("%s/data/SingleMuon-Run2015D-v1_toyMC.root",  channel.data()));
    treeData->Add(Form("%s/data/SingleMuon-Run2015D-v4_toyMC.root", channel.data()));

  }

  else return;

  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-100to200_13TeV_toyMC.root", channel.data()));
  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-200to400_13TeV_toyMC.root", channel.data()));
  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-400to600_13TeV_toyMC.root", channel.data()));
  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-600toInf_13TeV_toyMC.root", channel.data()));

  // To remove minor background contribution in data set (weight is -1)

  if( removeMinor ){

    treeData->Add(Form("%s/VV/WW_TuneCUETP8M1_13TeV_toyMC.root", channel.data()));
    treeData->Add(Form("%s/VV/WZ_TuneCUETP8M1_13TeV_toyMC.root", channel.data()));
    treeData->Add(Form("%s/VV/ZZ_TuneCUETP8M1_13TeV_toyMC.root", channel.data()));
    treeData->Add(Form("%s/TT/TT_TuneCUETP8M1_13TeV_toyMC.root", channel.data()));
    treeData->Add(Form("%s/ZH/ZH_HToBB_ZToLL_M125_13TeV_toyMC.root", channel.data()));

  }

  // Define all the variables from the trees

  RooRealVar cat ("cat", "", 0, 2);
  RooRealVar mJet("prmass", "M_{jet}",  30.,  300., "GeV");
  RooRealVar mZH ("mllbb",   "M_{ZH}", 900., 3000., "GeV");
  RooRealVar evWeight("evweight", "", -1.e3, 1.e3);

  // Set the range in zh mass 

  mZH.setRange("fullRange", 900., 3000.);

  RooBinning binsmZH(21, 900, 3000);

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

  // Create a dataset from a tree -> to process an unbinned likelihood fitting

  RooDataSet dataSetData   ("dataSetData",    "dataSetData",    variables, Cut(catCut),           WeightVar(evWeight), Import(*treeData));
  RooDataSet dataSetDataSB ("dataSetDataSB",  "dataSetDataSB",  variables, Cut(catCut && sbCut),  WeightVar(evWeight), Import(*treeData));
  RooDataSet dataSetDataSG ("dataSetDataSG",  "dataSetDataSG",  variables, Cut(catCut && sigCut), WeightVar(evWeight), Import(*treeData));
  RooDataSet dataSetZjets  ("dataSetZjets",   "dataSetZjets",   variables, Cut(catCut),           WeightVar(evWeight), Import(*treeZjets));
  RooDataSet dataSetZjetsSB("dataSetZjetsSB", "dataSetZjetsSB", variables, Cut(catCut && sbCut),  WeightVar(evWeight), Import(*treeZjets));  
  RooDataSet dataSetZjetsSG("dataSetZjetsSG", "dataSetZjetsSG", variables, Cut(catCut && sigCut), WeightVar(evWeight), Import(*treeZjets));
  
  // Total events number

  float totalMcEv   = dataSetZjets.sumEntries();
  float totalSBMcEv = dataSetZjetsSB.sumEntries();
  float totalSGMcEv = dataSetZjetsSG.sumEntries();
  float totalDataEv = dataSetData.sumEntries();

  RooRealVar nMcEvents  ("nMcEvents",   "nMcEvents",   0., 9999999.);
  RooRealVar nSBMcEvents("nSBMcEvents", "nSBMcEvents", 0., 9999999.);
  RooRealVar nSGMcEvents("nSGMcEvents", "nSGMcEvents", 0., 9999999.);
  RooRealVar nDataEvents("nDataEvents", "nDataEvents", 0., 9999999.);

  nMcEvents.setVal(totalMcEv);
  nMcEvents.setConstant(true);
  
  nSBMcEvents.setVal(totalSBMcEv);
  nSBMcEvents.setConstant(true);
  
  nSGMcEvents.setVal(totalSGMcEv);
  nSGMcEvents.setConstant(true);
  
  nDataEvents.setVal(totalDataEv);
  nDataEvents.setConstant(true);

  // Side band jet mass in data

  RooRealVar lamda("lamda", "lamda", -0.02,  -5.,   5.);

  RooRealVar constantSB("constantSB", "constantSB", -0.02,  -1.,   0.);
  RooRealVar offsetSB  ("offsetSB",   "offsetSB",      30, -50., 200.);
  RooRealVar widthSB   ("widthSB",    "widthSB",      100,   0., 200.);

  offsetSB.setConstant(true);

  //RooErfExpPdf model_mJetSB("model_mJetSB", "model_mJetSB", mJet, constantSB, offsetSB, widthSB);
  RooExponential model_mJetSB("model_mJetSB", "model_mJetSB", mJet, lamda);

  RooExtendPdf ext_model_mJetSB("ext_model_mJetSB", "ext_model_mJetSB", model_mJetSB, nSBMcEvents);

  RooFitResult* mJetSB_result = ext_model_mJetSB.fitTo(dataSetDataSB, SumW2Error(true), Extended(true), Range("lowSB,highSB"), Strategy(2), Minimizer("Minuit2"), Save(1));

  // Normalize factor to normalize the background in signal region of data

  RooAbsReal* nSIGFit = ext_model_mJetSB.createIntegral(RooArgSet(mJet), NormSet(mJet), Range("signal"));
  RooAbsReal* nSBFit  = ext_model_mJetSB.createIntegral(RooArgSet(mJet), NormSet(mJet), Range("lowSB,highSB"));

  float normFactor = dataSetDataSB.sumEntries()*(nSIGFit->getVal()/nSBFit->getVal());
  
  // Plot the results on a frame

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

  // Alpha ratio part

  RooRealVar a("a", "a",  0., -1.,    1.);
  RooRealVar b("b", "b", 1000,  0., 4000.);
  
  RooGenericPdf model_ZHSB("model_ZHSB", "model_ZHSB", "TMath::Exp(@1*@0+@2/@0)", RooArgSet(mZH,a,b));
  RooGenericPdf model_ZHSG("model_ZHSG", "model_ZHSG", "TMath::Exp(@1*@0+@2/@0)", RooArgSet(mZH,a,b));
  RooGenericPdf model_ZH  ("model_ZH",   "model_ZH",   "TMath::Exp(@1*@0+@2/@0)", RooArgSet(mZH,a,b));

  RooExtendPdf ext_model_ZHSB("ext_model_ZHSB", "ext_model_ZHSB", model_ZHSB, nSBMcEvents);
  RooExtendPdf ext_model_ZHSG("ext_model_ZHSG", "ext_model_ZHSG", model_ZHSG, nSGMcEvents);
  RooExtendPdf ext_model_ZH  ("ext_model_ZH",   "ext_model_ZH",   model_ZH,   nDataEvents);

  // Fit ZH mass in side band  

  RooFitResult* mZHSB_result = ext_model_ZHSB.fitTo(dataSetZjetsSB, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));

  RooAbsReal* nZHSBFit = ext_model_ZHSB.createIntegral(RooArgSet(mZH), NormSet(mZH), Range("fullRange"));

  float p0 = a.getVal();
  float p1 = b.getVal();

  // Fit ZH mass in signal region

  RooFitResult* mZHSG_result = ext_model_ZHSG.fitTo(dataSetZjetsSG, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));

  RooAbsReal* nZHSGFit = ext_model_ZHSG.createIntegral(RooArgSet(mZH), NormSet(mZH), Range("fullRange"));

  float p2 = a.getVal();
  float p3 = b.getVal();

  // Fit ZH mass in side band region (data)

  RooFitResult* mZH_result = ext_model_ZH.fitTo(dataSetDataSB, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));

  // Draw the model of alpha ratio
  // Multiply the model of background in data side band with the model of alpha ratio to the a model of background in data signal region

  RooGenericPdf model_alpha  ("model_alpha", "model_alpha", Form("TMath::Exp(%f*@0+%f/@0)/TMath::Exp(%f*@0+%f/@0)", p2,p3,p0,p1), RooArgSet(mZH));
  RooProdPdf    model_sigData("model_sigData", "ext_model_ZH*model_alpha", RooArgList(ext_model_ZH,model_alpha));

  // Plot the results to a frame 

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

  TLatex* lar = new TLatex();

  lar->SetTextSize(0.035);
  lar->SetLineWidth(5);
  
  TCanvas* c = new TCanvas("c","",0,0,1000,800);

  c->cd();
  mZHFrameMC->Draw();
  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.65, 0.60, Form("%s  %s btag", channel.data(), catcut.data()));
  c->Print(Form("rooFit_forData_%s_cat%s.pdf(", channel.data(), catcut.data()));

  c->cd();
  mZHFrame->Draw();
  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.65, 0.63, Form("%s  %s btag", channel.data(), catcut.data()));
  c->Print(Form("rooFit_forData_%s_cat%s.pdf", channel.data(), catcut.data()));

  c->cd();
  predictFrame->Draw();
  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.65, 0.63, Form("%s  %s btag", channel.data(), catcut.data()));
  c->Print(Form("rooFit_forData_%s_cat%s.pdf", channel.data(), catcut.data()));
  
  c->cd();
  mJetFrame->Draw();
  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.65, 0.63, Form("%s  %s btag", channel.data(), catcut.data()));
  c->Print(Form("rooFit_forData_%s_cat%s.pdf)", channel.data(), catcut.data()));

}
