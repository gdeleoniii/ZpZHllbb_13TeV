R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/PdfDiagonalizer_cc.so)
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/readFitParam.h"
using namespace RooFit;

void rooFitData(string channel, string catcut, bool removeMinor=true){

  // Suppress all the INFO message

  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooMsgService::instance().setSilentMode(true);
  gROOT->ProcessLine("gErrorIgnoreLevel=kWarning;");

  // Input files and sum all backgrounds

  TChain* treeData  = new TChain("tree");
  TChain* treeZjets = new TChain("tree");

  if( channel == "ele" ){

    treeData->Add("data/SingleElectron-Run2015D-v1_eleMiniTree.root");
    treeData->Add("data/SingleElectron-Run2015D-v4_eleMiniTree.root");

  }

  else if( channel == "mu" ){

    treeData->Add("data/SingleMuon-Run2015D-v1_muMiniTree.root");
    treeData->Add("data/SingleMuon-Run2015D-v4_muMiniTree.root");

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
  RooRealVar mZH ("mllbb", "M_{ZH}", 750., 4300., "GeV");
  RooRealVar evWeight("evweight", "", 0., 1.e3);

  // Set the range in zh mass and in jet mass

  mZH .setRange("fullRange", 750., 4300.);
  mJet.setRange("allRange", 30., 300.);
  mJet.setRange("lowSB",    30.,  65.);
  mJet.setRange("highSB",  135., 300.);
  mJet.setRange("signal",  105., 135.);

  RooBinning binsmZH(71, 750, 4300);
  RooBinning binsmJet(54, 30, 300);

  RooArgSet variables(cat, mJet, mZH, evWeight);

  TCut catCut = Form("cat==%s", catcut.data());
  TCut sbCut  = "prmass>30 && !(prmass>65 && prmass<135) && prmass<300";
  TCut sigCut = "prmass>105 && prmass<135";

  // Create a dataset from a tree -> to process an unbinned likelihood fitting

  RooDataSet dataSetData   ("dataSetData",    "dataSetData",    variables, Cut(catCut),           WeightVar(evWeight), Import(*treeData));
  RooDataSet dataSetDataSB ("dataSetDataSB",  "dataSetDataSB",  variables, Cut(catCut && sbCut),  WeightVar(evWeight), Import(*treeData));
  RooDataSet dataSetDataSG ("dataSetDataSG",  "dataSetDataSG",  variables, Cut(catCut && sigCut), WeightVar(evWeight), Import(*treeData));
  RooDataSet dataSetZjetsSB("dataSetZjetsSB", "dataSetZjetsSB", variables, Cut(catCut && sbCut),  WeightVar(evWeight), Import(*treeZjets));  
  RooDataSet dataSetZjetsSG("dataSetZjetsSG", "dataSetZjetsSG", variables, Cut(catCut && sigCut), WeightVar(evWeight), Import(*treeZjets));
  
  // Total events number

  RooRealVar nSBMcEvents  ("nSBMcEvents",   "nSBMcEvents",   0., 1.e9);
  RooRealVar nSGMcEvents  ("nSGMcEvents",   "nSGMcEvents",   0., 1.e9);
  RooRealVar nSBDataEvents("nSBDataEvents", "nSBDataEvents", 0., 1.e9);

  nSBMcEvents.setVal(dataSetZjetsSB.sumEntries());
  nSGMcEvents.setVal(dataSetZjetsSG.sumEntries());
  nSBDataEvents.setVal(dataSetDataSB.sumEntries());

  // Alpha ratio part  
  // Set fit parameters

  RooRealVar sbVara("sbVara", "sbVara", 1.00, 1.e4);
  RooRealVar sbVarb("sbVarb", "sbVarb", 0.01, 0.09);
  RooRealVar sgVara("sgVara", "sgVara", 1.00, 1.e4);
  RooRealVar sgVarb("sgVarb", "sgVarb", 0.01, 0.09);
  RooRealVar daVara("daVara", "daVara", 1.00, 1.e4);
  RooRealVar daVarb("daVarb", "daVarb", 0.01, 0.09);

  // Fix parameter "a"

  sbVara.setVal(108.3);
  sgVara.setVal(108.3);
  daVara.setVal(164.5);  

  sbVara.setConstant(true);
  sgVara.setConstant(true);
  daVara.setConstant(true);  

  // Fit ZH mass in MC side band

  RooGenericPdf model_ZHSB("model_ZHSB", "model_ZHSB", "exp(-@0/(@1+@2*@0))", RooArgSet(mZH,sbVara,sbVarb));
  RooExtendPdf  ext_model_ZHSB("ext_model_ZHSB", "ext_model_ZHSB", model_ZHSB, nSBMcEvents);
  RooFitResult* mZHSB_result = ext_model_ZHSB.fitTo(dataSetZjetsSB, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));

  // Fit ZH mass in MC signal region

  RooGenericPdf model_ZHSG("model_ZHSG", "model_ZHSG", "exp(-@0/(@1+@2*@0))", RooArgSet(mZH,sgVara,sgVarb));
  RooExtendPdf  ext_model_ZHSG("ext_model_ZHSG", "ext_model_ZHSG", model_ZHSG, nSGMcEvents);
  RooFitResult* mZHSG_result = ext_model_ZHSG.fitTo(dataSetZjetsSG, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));

  // Fit ZH mass in data side band

  RooGenericPdf model_ZH("model_ZH", "model_ZH", "exp(-@0/(@1+@2*@0))", RooArgSet(mZH,daVara,daVarb));
  RooExtendPdf  ext_model_ZH("ext_model_ZH", "ext_model_ZH", model_ZH, nSBDataEvents);
  RooFitResult* mZH_result = ext_model_ZH.fitTo(dataSetDataSB, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));

  // Multiply the model of background in data side band with the model of alpha ratio to the a model of background in data signal region

  float constant = ext_model_ZHSB.createIntegral(mZH)->getVal()/ext_model_ZHSG.createIntegral(mZH)->getVal();

  RooGenericPdf alpha_display("alpha_display", "alpha_display", Form("%f*exp(-@0/(%f+%f*@0))/exp(-@0/(%f+%f*@0))", constant, sgVara.getVal(), sgVarb.getVal(), sbVara.getVal(), sbVarb.getVal()), mZH);
  RooFormulaVar model_alpha("model_alpha", Form("%f*exp(-mllbb/(%f+%f*mllbb))/exp(-mllbb/(%f+%f*mllbb))", constant, sgVara.getVal(), sgVarb.getVal(), sbVara.getVal(), sbVarb.getVal()), mZH);
  RooEffProd    model_predicted("model_predicted", "model_predicted", model_ZH, model_alpha);

  // Fit jet mass in data side band

  RooRealVar     lamda("lamda", "lamda", -0.015, -0.04, -0.01);
  RooExponential model_mJetSB("model_mJetSB", "model_mJetSB", mJet, lamda);
  RooExtendPdf   ext_model_mJetSB("ext_model_mJetSB", "ext_model_mJetSB", model_mJetSB, nSBDataEvents);
  RooFitResult*  mJetSB_result = ext_model_mJetSB.fitTo(dataSetDataSB, SumW2Error(true), Extended(true), Range("lowSB,highSB"), Strategy(2), Minimizer("Minuit2"), Save(1));

  // Normalize factor to normalize the background in signal region of data

  RooAbsReal* nSIGFit = model_mJetSB.createIntegral(RooArgSet(mJet), NormSet(mJet), Range("signal"));
  RooAbsReal* nSBFit  = model_mJetSB.createIntegral(RooArgSet(mJet), NormSet(mJet), Range("lowSB,highSB"));

  RooRealVar normFactor("normFactor", "normFactor", 0., 1.e9);
  normFactor.setVal(nSBDataEvents.getVal()*(nSIGFit->getVal()/nSBFit->getVal()));

  fprintf(stdout, "sbVara=%f\nsbVarb=%f\nsgVara=%f\nsgVarb=%f\ndaVara=%f\ndaVarb=%f\nlamda=%f\n", sbVara.getVal(), sbVarb.getVal(), sgVara.getVal(), sgVarb.getVal(), daVara.getVal(), daVarb.getVal(), lamda.getVal());

  // Plot the results on frame 

  RooPlot* mcSBmZhFrame        = mZH.frame();
  RooPlot* mcSGmZhFrame        = mZH.frame();
  RooPlot* alphaFrame          = mZH.frame();
  RooPlot* dataSBmZhFrame      = mZH.frame();
  RooPlot* dataSBmJetFrame     = mJet.frame();
  RooPlot* expectedFrame       = mZH.frame(); 
  RooPlot* mcSBmZhPullFrame    = mZH.frame();
  RooPlot* mcSGmZhPullFrame    = mZH.frame();
  RooPlot* dataSBmZhPullFrame  = mZH.frame();
  RooPlot* dataSBmJetPullFrame = mJet.frame();

  dataSetZjetsSB.plotOn(mcSBmZhFrame, Binning(binsmZH));
  ext_model_ZHSB.plotOn(mcSBmZhFrame, VisualizeError(*mZHSB_result,1,false), FillStyle(3002));
  dataSetZjetsSB.plotOn(mcSBmZhFrame, Binning(binsmZH));
  ext_model_ZHSB.plotOn(mcSBmZhFrame, LineColor(kBlue));

  dataSetZjetsSG.plotOn(mcSGmZhFrame, Binning(binsmZH));
  ext_model_ZHSG.plotOn(mcSGmZhFrame, VisualizeError(*mZHSG_result,1,false), FillStyle(3002));
  dataSetZjetsSG.plotOn(mcSGmZhFrame, Binning(binsmZH));
  ext_model_ZHSG.plotOn(mcSGmZhFrame, LineColor(kBlue));
  
  ext_model_ZHSB.plotOn(alphaFrame, Normalization(1, RooAbsReal::NumEvent), LineColor(kBlue));
  ext_model_ZHSG.plotOn(alphaFrame, Normalization(1, RooAbsReal::NumEvent), LineColor(kRed));
  alpha_display .plotOn(alphaFrame, Normalization(1, RooAbsReal::NumEvent), LineColor(kBlack));

  dataSetDataSB.plotOn(dataSBmZhFrame, Binning(binsmZH));
  ext_model_ZH .plotOn(dataSBmZhFrame, VisualizeError(*mZH_result,1,false), FillStyle(3002));
  dataSetDataSB.plotOn(dataSBmZhFrame, Binning(binsmZH));
  ext_model_ZH .plotOn(dataSBmZhFrame, LineColor(kBlue));

  dataSetDataSB   .plotOn(dataSBmJetFrame, Binning(binsmJet));
  ext_model_mJetSB.plotOn(dataSBmJetFrame, Range("allRange"), VisualizeError(*mJetSB_result,1,false), FillStyle(3002));
  dataSetDataSB   .plotOn(dataSBmJetFrame, Binning(binsmJet));
  ext_model_mJetSB.plotOn(dataSBmJetFrame, Range("allRange"));

  dataSetDataSG.plotOn(expectedFrame, Binning(binsmZH));
  model_predicted.plotOn(expectedFrame, Normalization(normFactor.getVal(), RooAbsReal::NumEvent), LineColor(kRed+1));
  // Using RooAbsReal::NumEvent in order to consider the bin width of data set. Equivalent to (normFactor*binWidth) if using RooAbsReal::Raw.

  // Output the results

  TLatex lar;
  lar.SetTextSize(0.03);
  lar.SetLineWidth(5);

  float up_height = 0.82;
  float dw_height = (1-up_height)*1.445;

  TCanvas c0("c0","",0,0,1000,800);
  
  c0.Divide(1,2);

  TPad* c0_up = (TPad*)c0.GetListOfPrimitives()->FindObject("c0_1");
  TPad* c0_dw = (TPad*)c0.GetListOfPrimitives()->FindObject("c0_2"); 

  c0_up->SetPad(0,1-up_height,1,1);
  c0_dw->SetPad(0,0,1,dw_height);
  c0_dw->SetBottomMargin(0.25);
  c0_up->cd()->SetLogy(1);

  mcSBmZhFrame->SetTitle("");
  mcSBmZhFrame->SetMinimum(1e-4);
  mcSBmZhFrame->SetMaximum(catcut=="1"?1000:100);
  mcSBmZhFrame->GetXaxis()->SetTitle("");
  mcSBmZhFrame->GetXaxis()->SetLabelOffset(999);
  mcSBmZhFrame->Draw();

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.65, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s b-tag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, "Z+jets in sidebands");

  c0_up->RedrawAxis();
  c0_dw->cd()->SetLogy(0);

  mcSBmZhPullFrame->addObject(mcSBmZhFrame->pullHist(), "P");
  mcSBmZhPullFrame->SetTitle("");
  mcSBmZhPullFrame->GetYaxis()->SetTitle("Pulls");
  mcSBmZhPullFrame->GetYaxis()->SetTitleOffset(0.25);
  mcSBmZhPullFrame->GetXaxis()->SetLabelSize(0.125);
  mcSBmZhPullFrame->GetXaxis()->SetTitleSize(0.125);
  mcSBmZhPullFrame->GetYaxis()->SetLabelSize(0.125);
  mcSBmZhPullFrame->GetYaxis()->SetTitleSize(0.125);
  mcSBmZhPullFrame->GetYaxis()->SetNdivisions(505);
  mcSBmZhPullFrame->SetMinimum(-4);
  mcSBmZhPullFrame->SetMaximum(4);
  mcSBmZhPullFrame->Draw();

  c0.Draw();
  c0.Print(Form("rooFit_forData_%s_cat%s.pdf(", channel.data(), catcut.data()));

  TCanvas c1("c1","",0,0,1000,800);

  c1.Divide(1,2);

  TPad* c1_up = (TPad*)c1.GetListOfPrimitives()->FindObject("c1_1");
  TPad* c1_dw = (TPad*)c1.GetListOfPrimitives()->FindObject("c1_2");

  c1_up->SetPad(0,1-up_height,1,1);
  c1_dw->SetPad(0,0,1,dw_height);
  c1_dw->SetBottomMargin(0.25);
  c1_up->cd()->SetLogy(1);

  mcSGmZhFrame->SetTitle("");
  mcSGmZhFrame->SetMinimum(1e-4);
  mcSGmZhFrame->SetMaximum(catcut=="1"?100:10);
  mcSGmZhFrame->GetXaxis()->SetTitle("");
  mcSGmZhFrame->GetXaxis()->SetLabelOffset(999);
  mcSGmZhFrame->Draw();

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.65, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s b-tag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, "Z+jets in signal region");

  c1_up->RedrawAxis();
  c1_dw->cd()->SetLogy(0);

  mcSGmZhPullFrame->addObject(mcSGmZhFrame->pullHist(), "P");
  mcSGmZhPullFrame->SetTitle("");
  mcSGmZhPullFrame->GetYaxis()->SetTitle("Pulls");
  mcSGmZhPullFrame->GetYaxis()->SetTitleOffset(0.25);
  mcSGmZhPullFrame->GetXaxis()->SetLabelSize(0.125);
  mcSGmZhPullFrame->GetXaxis()->SetTitleSize(0.125);
  mcSGmZhPullFrame->GetYaxis()->SetLabelSize(0.125);
  mcSGmZhPullFrame->GetYaxis()->SetTitleSize(0.125);
  mcSGmZhPullFrame->GetYaxis()->SetNdivisions(505);
  mcSGmZhPullFrame->SetMinimum(-4);
  mcSGmZhPullFrame->SetMaximum(4);
  mcSGmZhPullFrame->Draw();

  c1.Draw();
  c1.Print(Form("rooFit_forData_%s_cat%s.pdf", channel.data(), catcut.data()));  

  TCanvas c2("c2","",0,0,1000,800);
  
  c2.Divide(1,2);

  TPad* c2_up = (TPad*)c2.GetListOfPrimitives()->FindObject("c2_1");
  TPad* c2_dw = (TPad*)c2.GetListOfPrimitives()->FindObject("c2_2"); 

  c2_up->SetPad(0,1-up_height,1,1);
  c2_dw->SetPad(0,0,1,dw_height);
  c2_dw->SetBottomMargin(0.25);
  c2_up->cd()->SetLogy(1);

  dataSBmZhFrame->SetTitle("");
  dataSBmZhFrame->SetMinimum(1e-4);
  dataSBmZhFrame->SetMaximum(catcut=="1"?100:10);
  dataSBmZhFrame->GetXaxis()->SetTitle("");
  dataSBmZhFrame->GetXaxis()->SetLabelOffset(999);
  dataSBmZhFrame->Draw();

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar.DrawLatexNDC(0.65, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s btag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, "data in sidebands");

  c2_up->RedrawAxis();
  c2_dw->cd()->SetLogy(0);

  dataSBmZhPullFrame->addObject(dataSBmZhFrame->pullHist(), "P");
  dataSBmZhPullFrame->SetTitle("");
  dataSBmZhPullFrame->GetYaxis()->SetTitle("Pulls");
  dataSBmZhPullFrame->GetYaxis()->SetTitleOffset(0.25);
  dataSBmZhPullFrame->GetXaxis()->SetLabelSize(0.125);
  dataSBmZhPullFrame->GetXaxis()->SetTitleSize(0.125);
  dataSBmZhPullFrame->GetYaxis()->SetLabelSize(0.125);
  dataSBmZhPullFrame->GetYaxis()->SetTitleSize(0.125);
  dataSBmZhPullFrame->GetYaxis()->SetNdivisions(505);
  dataSBmZhPullFrame->SetMinimum(-4);
  dataSBmZhPullFrame->SetMaximum(4);
  dataSBmZhPullFrame->Draw();

  c2.Draw();
  c2.Print(Form("rooFit_forData_%s_cat%s.pdf", channel.data(), catcut.data()));

  TCanvas c3("c3","",0,0,1000,800);
  
  c3.Divide(1,2);

  TPad* c3_up = (TPad*)c3.GetListOfPrimitives()->FindObject("c3_1");
  TPad* c3_dw = (TPad*)c3.GetListOfPrimitives()->FindObject("c3_2"); 

  c3_up->SetPad(0,1-up_height,1,1);
  c3_dw->SetPad(0,0,1,dw_height);
  c3_dw->SetBottomMargin(0.25);
  c3_up->cd()->SetLogy(1);

  dataSBmJetFrame->SetTitle("");
  dataSBmJetFrame->SetMinimum(1e-2);
  dataSBmJetFrame->SetMaximum(100);
  dataSBmJetFrame->GetXaxis()->SetTitle("");
  dataSBmJetFrame->GetXaxis()->SetLabelOffset(999);
  dataSBmJetFrame->Draw();

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar.DrawLatexNDC(0.60, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s btag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, "data jet mass in sidebands");
  lar.DrawLatexNDC(0.15, 0.78, Form("Normalization factor: %.3f", normFactor.getVal()));
  
  c3_up->RedrawAxis();
  c3_dw->cd()->SetLogy(0);

  dataSBmJetPullFrame->addObject(dataSBmJetFrame->pullHist(), "P");
  dataSBmJetPullFrame->SetTitle("");
  dataSBmJetPullFrame->GetYaxis()->SetTitle("Pulls");
  dataSBmJetPullFrame->GetYaxis()->SetTitleOffset(0.25);
  dataSBmJetPullFrame->GetXaxis()->SetLabelSize(0.125);
  dataSBmJetPullFrame->GetXaxis()->SetTitleSize(0.125);
  dataSBmJetPullFrame->GetYaxis()->SetLabelSize(0.125);
  dataSBmJetPullFrame->GetYaxis()->SetTitleSize(0.125);
  dataSBmJetPullFrame->GetYaxis()->SetNdivisions(505);
  dataSBmJetPullFrame->SetMinimum(-4);
  dataSBmJetPullFrame->SetMaximum(4);
  dataSBmJetPullFrame->Draw();

  c3.Draw();
  c3.Print(Form("rooFit_forData_%s_cat%s.pdf", channel.data(), catcut.data()));

  TCanvas cv("cv","",0,0,1000,800);
  TLegend leg(0.60,0.70,0.85,0.80);

  cv.cd();
  alphaFrame->SetTitle("");
  alphaFrame->GetYaxis()->SetTitle("Normalization");
  alphaFrame->GetYaxis()->SetTitleOffset(1.3);
  alphaFrame->Draw();
  leg.AddEntry(alphaFrame->findObject(alphaFrame->nameOf(0)), "bkg. fit in sidebands", "l");
  leg.AddEntry(alphaFrame->findObject(alphaFrame->nameOf(1)), "bkg. fit in signal region", "l");
  leg.AddEntry(alphaFrame->findObject(alphaFrame->nameOf(2)), "#alpha function (y=exp(#frac{-x}{a+bx}))", "l");
  leg.SetBorderSize(0);
  leg.Draw();
  alphaFrame->addObject(&leg);
  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.60, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.62, 0.82, Form("%s, %s b-tag", channel.data(), catcut.data()));
  cv.Print(Form("rooFit_forData_%s_cat%s.pdf", channel.data(), catcut.data()));

  cv.Clear();
  cv.cd()->SetLogy(1);
  expectedFrame->SetTitle("");
  expectedFrame->SetMinimum(1e-4);
  expectedFrame->SetMaximum(catcut=="1"?100:10);
  expectedFrame->Draw();
  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar.DrawLatexNDC(0.60, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s btag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, "expected background in data signal region");
  cv.Print(Form("rooFit_forData_%s_cat%s.pdf)", channel.data(), catcut.data()));

}
