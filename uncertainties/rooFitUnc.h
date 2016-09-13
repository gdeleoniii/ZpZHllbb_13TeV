R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/PdfDiagonalizer_cc.so)
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/readFitParam.h"
using namespace RooFit;

void rooFitUnc(string channel, string catcut, string region, TF1** f_alpha, TH1** h_shape, int i){

  // Suppress all the INFO message

  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooMsgService::instance().setSilentMode(true);
  gROOT->ProcessLine("gErrorIgnoreLevel=kWarning;");

  // Input files and sum all backgrounds

  TChain* treeData  = new TChain("tree");
  TChain* treeZjets = new TChain("tree");

  if( channel == "ele" ){

    treeData->Add("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/data/SingleElectron-Run2015D-v1_eleMiniTree.root");
    treeData->Add("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/data/SingleElectron-Run2015D-v4_eleMiniTree.root");

  }

  else if( channel == "mu" ){

    treeData->Add("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/data/SingleMuon-Run2015D-v1_muMiniTree.root");
    treeData->Add("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/data/SingleMuon-Run2015D-v4_muMiniTree.root");

  }

  else return;

  treeZjets->Add(Form("Zjets/DYJetsToLL_M-50_HT-100to200_13TeV_%s_%sMiniTree.root", region.data(), channel.data()));
  treeZjets->Add(Form("Zjets/DYJetsToLL_M-50_HT-200to400_13TeV_%s_%sMiniTree.root", region.data(), channel.data()));
  treeZjets->Add(Form("Zjets/DYJetsToLL_M-50_HT-400to600_13TeV_%s_%sMiniTree.root", region.data(), channel.data()));
  treeZjets->Add(Form("Zjets/DYJetsToLL_M-50_HT-600toInf_13TeV_%s_%sMiniTree.root", region.data(), channel.data()));

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

  param myVal(channel.data(), catcut.data());

  RooRealVar sbVara("sbVara", "sbVara", myVal.value("sbVara"), myVal.value("sbVaraMin"), myVal.value("sbVaraMax"));
  RooRealVar sbVarb("sbVarb", "sbVarb", myVal.value("sbVarb"), myVal.value("sbVarbMin"), myVal.value("sbVarbMax"));
  RooRealVar sgVara("sgVara", "sgVara", myVal.value("sgVara"), myVal.value("sgVaraMin"), myVal.value("sgVaraMax"));
  RooRealVar sgVarb("sgVarb", "sgVarb", myVal.value("sgVarb"), myVal.value("sgVarbMin"), myVal.value("sgVarbMax"));
  RooRealVar daVara("daVara", "daVara", myVal.value("daVara"), myVal.value("daVaraMin"), myVal.value("daVaraMax"));
  RooRealVar daVarb("daVarb", "daVarb", myVal.value("daVarb"), myVal.value("daVarbMin"), myVal.value("daVarbMax"));

  // ZH mass in MC side band

  RooGenericPdf model_ZHSB("model_ZHSB", "model_ZHSB", "exp(-@0/(@1+@2*@0))", RooArgSet(mZH,sbVara,sbVarb));
  RooExtendPdf  ext_model_ZHSB("ext_model_ZHSB", "ext_model_ZHSB", model_ZHSB, nSBMcEvents);

  // ZH mass in MC signal region

  RooGenericPdf model_ZHSG("model_ZHSG", "model_ZHSG", "exp(-@0/(@1+@2*@0))", RooArgSet(mZH,sgVara,sgVarb));
  RooExtendPdf  ext_model_ZHSG("ext_model_ZHSG", "ext_model_ZHSG", model_ZHSG, nSGMcEvents);

  // ZH mass in data side band 

  RooGenericPdf model_ZH("model_ZH", "model_ZH", "exp(-@0/(@1+@2*@0))", RooArgSet(mZH,daVara,daVarb));
  RooExtendPdf  ext_model_ZH("ext_model_ZH", "ext_model_ZH", model_ZH, nSBDataEvents);

  // Make category to fit signal and sideband

  RooCategory samples("samples", "samples");

  samples.defineType("mcSideband");
  samples.defineType("mcSignal");
  samples.defineType("dataSideband");

  RooDataSet dataSetCombine("dataSetCombine", "dataSetCombine", variables, Index(samples), Import("mcSideband", dataSetZjetsSB), Import("mcSignal", dataSetZjetsSG), Import("dataSideband", dataSetDataSB), WeightVar(evWeight));

  RooSimultaneous modelCombine("modelCombine", "modelCombine", samples);

  modelCombine.addPdf(ext_model_ZHSB, "mcSideband");
  modelCombine.addPdf(ext_model_ZHSG, "mcSignal");
  modelCombine.addPdf(ext_model_ZH,   "dataSideband");

  RooFitResult* combineResult = modelCombine.fitTo(dataSetCombine, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));

  // Multiply the model of background in data side band with the model of alpha ratio to the a model of background in data signal region

  float constant = ext_model_ZHSB.createIntegral(mZH)->getVal()/ext_model_ZHSG.createIntegral(mZH)->getVal();

  RooFormulaVar model_alpha("model_alpha", Form("%f*exp(-mllbb/(%f+%f*mllbb))/exp(-mllbb/(%f+%f*mllbb))", constant, sgVara.getVal(), sgVarb.getVal(), sbVara.getVal(), sbVarb.getVal()), mZH);
  RooEffProd    model_predicted("model_predicted", "model_predicted", model_ZH, model_alpha);

  *f_alpha = new TF1(Form("f_alpha%i",i), "[0]*TMath::Exp(-x/([1]+[2]*x))/TMath::Exp(-x/([3]+[4]*x))", 750, 4300);

  (*f_alpha)->SetParameters(constant, sgVara.getVal(), sgVarb.getVal(), sbVara.getVal(), sbVarb.getVal());

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

  *h_shape = model_predicted.createHistogram(Form("h_shape%i",i), mZH, Binning(binsmZH), Extended(true));

}
