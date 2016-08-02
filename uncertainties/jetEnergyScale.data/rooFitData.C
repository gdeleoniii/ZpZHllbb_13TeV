R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/PdfDiagonalizer_cc.so)
using namespace RooFit;

void rooFitData(string channel, string catcut, bool removeMinor=true){

  // Suppress all the INFO message

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooMsgService::instance().setSilentMode(true);

  // Input files and sum all backgrounds

  TChain* treeData  = new TChain("tree");

  if( channel == "ele" ){

    treeData->Add(Form("%s/data/SingleElectron-Run2015D-v1_toyMC.root",  channel.data()));
    treeData->Add(Form("%s/data/SingleElectron-Run2015D-v4_toyMC.root", channel.data()));

  }

  else if( channel == "mu" ){

    treeData->Add(Form("%s/data/SingleMuon-Run2015D-v1_toyMC.root",  channel.data()));
    treeData->Add(Form("%s/data/SingleMuon-Run2015D-v4_toyMC.root", channel.data()));

  }

  else return;

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
  RooRealVar mJet("prmass", "M_{jet}", 30., 300., "GeV");
  RooRealVar mZH ("mllbb", "M_{ZH}", 800., 4000., "GeV");
  RooRealVar evWeight("evweight", "", 0., 1.e3);

  // Set the range in zh mass and in jet mass

  mJet.setRange("allRange", 30., 300.);
  mJet.setRange("lowSB",    30.,  65.);
  mJet.setRange("highSB",  135., 300.);
  mJet.setRange("signal",  105., 135.);

  RooBinning binsmJet(27, 30, 300);
  RooArgSet variables(cat, mJet, mZH, evWeight);

  TCut catCut = Form("cat==%s", catcut.c_str());
  TCut sbCut  = "prmass>30 && !(prmass>65 && prmass<135) && prmass<300";

  // Create a dataset from a tree -> to process an unbinned likelihood fitting
 
  RooDataSet dataSetDataSB ("dataSetDataSB",  "dataSetDataSB",  variables, Cut(catCut && sbCut),  WeightVar(evWeight), Import(*treeData));
 
  // Side band jet mass in data

  /*
    RooRealVar constantSB("constantSB", "constantSB", -0.02,  -1.,   0.);
    RooRealVar offsetSB  ("offsetSB",   "offsetSB",      30, -50., 200.);
    RooRealVar widthSB   ("widthSB",    "widthSB",      100,   0., 200.);
    offsetSB.setConstant(true);
    RooErfExpPdf model_mJetSB("model_mJetSB", "model_mJetSB", mJet, constantSB, offsetSB, widthSB);
  */

  RooRealVar lamda("lamda", "lamda", -0.02, -0.5, -0.0001);
  RooExponential model_mJetSB("model_mJetSB", "model_mJetSB", mJet, lamda);
  RooExtendPdf ext_model_mJetSB("ext_model_mJetSB", "ext_model_mJetSB", model_mJetSB, nSBMcEvents);
  RooFitResult* mJetSB_result = ext_model_mJetSB.fitTo(dataSetDataSB, SumW2Error(true), Extended(true), Range("lowSB,highSB"), Strategy(2), Minimizer("Minuit2"), Save(1));

  // Normalize factor to normalize the background in signal region of data

  RooAbsReal* nSIGFit = ext_model_mJetSB.createIntegral(RooArgSet(mJet), NormSet(mJet), Range("signal"));
  RooAbsReal* nSBFit  = ext_model_mJetSB.createIntegral(RooArgSet(mJet), NormSet(mJet), Range("lowSB,highSB"));

  float normFactor = dataSetDataSB.sumEntries()*(nSIGFit->getVal()/nSBFit->getVal());
 
  fprintf(stdout, "The normalization factor is %g\n", normFactor);

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

  TLatex* lar = new TLatex();

  lar->SetTextSize(0.035);
  lar->SetLineWidth(5);

  TCanvas* cv = new TCanvas("cv","",0,0,1000,800);

  cv->cd();
  mJetFrame->Draw();
  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.65, 0.63, Form("%s  %s btag", channel.data(), catcut.data()));
  cv->Print(Form("rooFit_forData_%s_cat%s.pdf", channel.data(), catcut.data()));

}
