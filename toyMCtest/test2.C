R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/PdfDiagonalizer_cc.so)
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/readFitParam.h"
using namespace RooFit;

void test2(string channel, string catcut, bool removeMinor=true){

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
  // RooFitResult* mZHSB_result = ext_model_ZHSB.fitTo(dataSetZjetsSB, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));

  // Fit ZH mass in MC signal region

  RooGenericPdf model_ZHSG("model_ZHSG", "model_ZHSG", "exp(-@0/(@1+@2*@0))", RooArgSet(mZH,sgVara,sgVarb));
  RooExtendPdf  ext_model_ZHSG("ext_model_ZHSG", "ext_model_ZHSG", model_ZHSG, nSGMcEvents);
  // RooFitResult* mZHSG_result = ext_model_ZHSG.fitTo(dataSetZjetsSG, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));

  // Fit ZH mass in data side band

  RooGenericPdf model_ZH("model_ZH", "model_ZH", "exp(-@0/(@1+@2*@0))", RooArgSet(mZH,daVara,daVarb));
  RooExtendPdf  ext_model_ZH("ext_model_ZH", "ext_model_ZH", model_ZH, nSBDataEvents);
  // RooFitResult* mZH_result = ext_model_ZH.fitTo(dataSetDataSB, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));





  RooCategory sample("sample", "sample");

  sample.defineType("sideband");
  sample.defineType("signal");

  RooDataSet comb("comb", "comb",  variables, Index(sample), Import("sideband", dataSetZjetsSB), Import("signal", dataSetZjetsSG), WeightVar(evWeight));

  RooSimultaneous simPdf("simPdf", "simPdf", sample);

  simPdf.addPdf(ext_model_ZHSB, "sideband");
  simPdf.addPdf(ext_model_ZHSG, "signal");

  RooFitResult* results = simPdf.fitTo(comb, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));


  cout << simPdf.getPropagatedError(*results) << endl;



  fprintf(stdout, "sbVara=%f\nsbVarb=%f\nsgVara=%f\nsgVarb=%f\ndaVara=%f\ndaVarb=%f\n", sbVara.getVal(), sbVarb.getVal(), sgVara.getVal(), sgVarb.getVal(), daVara.getVal(), daVarb.getVal());

  RooPlot* fframe = mZH.frame();
  RooPlot* xframe = mZH.frame();

  comb.plotOn(fframe, Cut("sample==sample::signal"), Binning(binsmZH));
  simPdf.plotOn(fframe, Slice(sample, "signal"), ProjWData(sample, comb), VisualizeError(*results,1,false), FillStyle(3002));
  //comb.plotOn(fframe, Cut("sample==sample::signal"), Binning(binsmZH));
  //simPdf.plotOn(fframe, Slice(sample, "signal"), ProjWData(sample, comb));
  //comb.plotOn(xframe, Cut("sample==sample::sideband"), Binning(binsmZH));
  dataSetZjetsSG.plotOn(xframe, Binning(binsmZH));


  TCanvas c0("c0","",0,0,1000,800);
  c0.cd()->SetLogy(1);

  fframe->SetMinimum(1e-4);
  fframe->SetMaximum(1e3);
  fframe->Draw();

  c0.Print("test.pdf(");

  c0.Clear();
  c0.cd()->SetLogy(1);

  xframe->SetMinimum(1e-4);
  xframe->SetMaximum(1e3);
  xframe->Draw();

  c0.Print("test.pdf)");




  // Multiply the model of background in data side band with the model of alpha ratio to the a model of background in data signal region
  /*
  float constant = ext_model_ZHSB.createIntegral(mZH)->getVal()/ext_model_ZHSG.createIntegral(mZH)->getVal();

  RooGenericPdf alpha_display("alpha_display", "alpha_display", Form("%f*exp(-@0/(%f+%f*@0))/exp(-@0/(%f+%f*@0))", constant, sgVara.getVal(), sgVarb.getVal(), sbVara.getVal(), sbVarb.getVal()), mZH);
  RooFormulaVar model_alpha("model_alpha", Form("%f*exp(-mllbb/(%f+%f*mllbb))/exp(-mllbb/(%f+%f*mllbb))", constant, sgVara.getVal(), sgVarb.getVal(), sbVara.getVal(), sbVarb.getVal()), mZH);
  RooEffProd    model_predicted("model_predicted", "model_predicted", model_ZH, model_alpha);
  */



}
