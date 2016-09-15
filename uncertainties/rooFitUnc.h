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

  TChain* tree_Data = new TChain("tree");
  TChain* tree_Dom  = new TChain("tree");
  TChain* tree_Sub  = new TChain("tree");

  // Data

  if( channel == "ele" ){

    tree_Data->Add("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/data/SingleElectron-Run2015D-v1_eleMiniTree.root");
    tree_Data->Add("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/data/SingleElectron-Run2015D-v4_eleMiniTree.root");

  }

  else if( channel == "mu" ){

    tree_Data->Add("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/data/SingleMuon-Run2015D-v1_muMiniTree.root");
    tree_Data->Add("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/data/SingleMuon-Run2015D-v4_muMiniTree.root");

  }

  // Dominant background

  tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-100to200_13TeV_%sMiniTree.root", channel.data()));
  tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-200to400_13TeV_%sMiniTree.root", channel.data()));
  tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-400to600_13TeV_%sMiniTree.root", channel.data()));
  tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-600toInf_13TeV_%sMiniTree.root", channel.data()));

  // Subdominant background

  tree_Sub->Add(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/minor/WW_TuneCUETP8M1_13TeV_%sMiniTree.root",     channel.data()));
  tree_Sub->Add(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/minor/WZ_TuneCUETP8M1_13TeV_%sMiniTree.root",     channel.data()));
  tree_Sub->Add(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/minor/ZZ_TuneCUETP8M1_13TeV_%sMiniTree.root",     channel.data()));
  tree_Sub->Add(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/minor/TT_TuneCUETP8M1_13TeV_%sMiniTree.root",     channel.data()));
  tree_Sub->Add(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/minor/ZH_HToBB_ZToLL_M125_13TeV_%sMiniTree.root", channel.data()));

  // Define all the variables from the trees

  RooRealVar cat     ("cat", "", 0, 2);
  RooRealVar mJet    ("prmass", "M_{jet}", 30., 300., "GeV");
  RooRealVar mZH     ("mllbb", "M_{ZH}", 750., 4300., "GeV");
  RooRealVar evWeight("evweight", "", 0., 1.e3);

  // Set the range in zh mass and in jet mass

  mZH .setRange("all",    750., 4300.);
  mJet.setRange("all",    30.,  300.);
  mJet.setRange("lowSB",  30.,  65.);
  mJet.setRange("highSB", 135.,  300.);
  mJet.setRange("signal", 105.,  135.);

  RooBinning bin_mZH (71, 750, 4300);
  RooBinning bin_mJet(54, 30, 300);

  RooArgSet variables(cat, mJet, mZH, evWeight);

  TCut catCut = Form("cat==%s", catcut.data());
  TCut sbCut  = "prmass>30 && !(prmass>65 && prmass<135) && prmass<300";
  TCut sigCut = "prmass>105 && prmass<135";

  // Create a dataset from a tree -> to process an unbinned likelihood fitting

  RooDataSet set_Data  ("set_Data",   "set_Data",   variables, Cut(catCut),           WeightVar(evWeight), Import(*tree_Data));
  RooDataSet set_sbData("set_sbData", "set_sbData", variables, Cut(catCut && sbCut),  WeightVar(evWeight), Import(*tree_Data));
  RooDataSet set_sgData("set_sgData", "set_sgData", variables, Cut(catCut && sigCut), WeightVar(evWeight), Import(*tree_Data));
  RooDataSet set_sbDom ("set_sbDom",  "set_sbDom",  variables, Cut(catCut && sbCut),  WeightVar(evWeight), Import(*tree_Dom));  
  RooDataSet set_sgDom ("set_sgDom",  "set_sgDom",  variables, Cut(catCut && sigCut), WeightVar(evWeight), Import(*tree_Dom));
  RooDataSet set_sbSub ("set_sbSub",  "set_sbSub",  variables, Cut(catCut && sbCut),  WeightVar(evWeight), Import(*tree_Sub));
  RooDataSet set_sgSub ("set_sgSub",  "set_sgSub",  variables, Cut(catCut && sigCut), WeightVar(evWeight), Import(*tree_Sub));
  
  // Total events number

  RooRealVar nEv_sbDom ("nEv_sbDom",  "nEv_sbDom",  0., 1.e9);
  RooRealVar nEv_sgDom ("nEv_sgDom",  "nEv_sgDom",  0., 1.e9);
  RooRealVar nEv_sbSub ("nEv_sbSub",  "nEv_sbSub",  0., 1.e9);
  RooRealVar nEv_sgSub ("nEv_sgSub",  "nEv_sgSub",  0., 1.e9);
  RooRealVar nEv_sbData("nEv_sbData", "nEv_sbData", 0., 1.e9);

  nEv_sbDom .setVal(set_sbDom.sumEntries());
  nEv_sgDom .setVal(set_sbDom.sumEntries());
  nEv_sbSub .setVal(set_sbSub .sumEntries());
  nEv_sgSub .setVal(set_sgSub .sumEntries());
  nEv_sbData.setVal(set_sbData.sumEntries());

  // Set fit parameters

  param myVal(channel.data(), catcut.data());

  RooRealVar a_domSb("a_domSb", "a_domSb", myVal.value("a_domSb"), myVal.value("a_domSbMin"), myVal.value("a_domSbMax"));
  RooRealVar b_domSb("b_domSb", "b_domSb", myVal.value("b_domSb"), myVal.value("b_domSbMin"), myVal.value("b_domSbMax"));
  RooRealVar a_domSg("a_domSg", "a_domSg", myVal.value("a_domSg"), myVal.value("a_domSgMin"), myVal.value("a_domSgMax"));
  RooRealVar b_domSg("b_domSg", "b_domSg", myVal.value("b_domSg"), myVal.value("b_domSgMin"), myVal.value("b_domSgMax"));
  RooRealVar a_subSb("a_subSb", "a_subSb", myVal.value("a_subSb"), myVal.value("a_subSbMin"), myVal.value("a_subSbMax"));
  RooRealVar b_subSb("b_subSb", "b_subSb", myVal.value("b_subSb"), myVal.value("b_subSbMin"), myVal.value("b_subSbMax"));
  RooRealVar a_subSg("a_subSg", "a_subSg", myVal.value("a_subSg"), myVal.value("a_subSgMin"), myVal.value("a_subSgMax"));
  RooRealVar b_subSg("b_subSg", "b_subSg", myVal.value("b_subSg"), myVal.value("b_subSgMin"), myVal.value("b_subSgMax"));
  RooRealVar a_datSb("a_datSb", "a_datSb", myVal.value("a_datSb"), myVal.value("a_datSbMin"), myVal.value("a_datSbMax"));
  RooRealVar b_datSb("b_datSb", "b_datSb", myVal.value("b_datSb"), myVal.value("b_datSbMin"), myVal.value("b_datSbMax"));

  // Create pdf for ZH mass

  RooGenericPdf pdf_sbDomZh ("pdf_sbDomZh",  "pdf_sbDomZh",  "exp(-@0/(@1+@2*@0))", RooArgSet(mZH,a_domSb,b_domSb));
  RooGenericPdf pdf_sgDomZh ("pdf_sgDomZh",  "pdf_sgDomZh",  "exp(-@0/(@1+@2*@0))", RooArgSet(mZH,a_domSg,b_domSg));
  RooGenericPdf pdf_sbSubZh ("pdf_sbSubZh",  "pdf_sbSubZh",  "exp(-@0/(@1+@2*@0))", RooArgSet(mZH,a_subSb,b_subSb));
  RooGenericPdf pdf_sgSubZh ("pdf_sgSubZh",  "pdf_sgSubZh",  "exp(-@0/(@1+@2*@0))", RooArgSet(mZH,a_subSg,b_subSg));
  RooGenericPdf pdf_sbDataZh("pdf_sbDataZh", "pdf_sbDataZh", "exp(-@0/(@1+@2*@0))", RooArgSet(mZH,a_datSb,b_datSb));

  // Extended pdf from RooGenericPdf

  RooExtendPdf ext_sbDomZh ("ext_sbDomZh",  "ext_sbDomZh",  pdf_sbDomZh,  nEv_sbDom);
  RooExtendPdf ext_sgDomZh ("ext_sgDomZh",  "ext_sgDomZh",  pdf_sgDomZh,  nEv_sgDom);
  RooExtendPdf ext_sbSubZh ("ext_sbSubZh",  "ext_sbSubZh",  pdf_sbSubZh,  nEv_sbSub);
  RooExtendPdf ext_sgSubZh ("ext_sgSubZh",  "ext_sgSubZh",  pdf_sgSubZh,  nEv_sgSub);
  RooExtendPdf ext_sbDataZh("ext_sbDataZh", "ext_sbDataZh", pdf_sbDataZh, nEv_sbData);

  // Make category to fit dominant background signal/sideband and data sideband

  RooCategory cat_domData("cat_domData", "cat_domData");

  cat_domData.defineType("dom_SB");
  cat_domData.defineType("dom_SG");
  cat_domData.defineType("data_SB");

  RooDataSet cmb_domData("cmb_domData", "cmb_domData", variables, Index(cat_domData), Import("dom_SB", set_sbDom), Import("dom_SG", set_sgDom), Import("data_SB", set_sbData), WeightVar(evWeight));

  RooSimultaneous pdf_domData("pdf_domData", "pdf_domData", cat_domData);

  pdf_domData.addPdf(ext_sbDomZh, "dom_SB");
  pdf_domData.addPdf(ext_sgDomZh, "dom_SG");
  pdf_domData.addPdf(ext_sbDataZh,"data_SB");

  RooFitResult* res_domData = pdf_domData.fitTo(cmb_domData, SumW2Error(true), Extended(true), Range("all"), Strategy(2), Minimizer("Minuit2"), Save(1));

  // Make category to fit subdominant background signal/sideband

  RooCategory cat_sub("cat_sub", "cat_sub");

  cat_sub.defineType("sub_SB");
  cat_sub.defineType("sub_SG");

  RooDataSet cmb_sub("cmb_sub", "cmb_sub", variables, Index(cat_sub), Import("sub_SB", set_sbSub), Import("sub_SG", set_sgSub), WeightVar(evWeight)); 

  RooSimultaneous pdf_sub("pdf_sub", "pdf_sub", cat_sub);

  pdf_sub.addPdf(ext_sbSubZh, "sub_SB");
  pdf_sub.addPdf(ext_sgSubZh, "sub_SG");

  RooFitResult* res_sub = pdf_sub.fitTo(cmb_sub, SumW2Error(true), Extended(true), Range("all"), Strategy(2), Minimizer("Minuit2"), Save(1));

  // Multiply the model of background in data side band with the model of alpha ratio to the a model of background in data signal region

  float constant = ext_sbDomZh.createIntegral(mZH)->getVal()/ext_sgDomZh.createIntegral(mZH)->getVal();

  RooGenericPdf alpha_display("alpha_display", "alpha_display", Form("%f*exp(-@0/(@1+@2*@0))/exp(-@0/(@3+@4*@0))", constant), RooArgSet(mZH,a_domSg,b_domSg,a_domSb,b_domSb)); 
  RooFormulaVar pdf_alpha("pdf_alpha", Form("%f*exp(-mllbb/(%f+%f*mllbb))/exp(-mllbb/(%f+%f*mllbb))", constant, a_domSg.getVal(), b_domSg.getVal(), a_domSb.getVal(), b_domSb.getVal()), mZH);
  RooEffProd    pdf_predicted("pdf_predicted", "pdf_predicted", pdf_sbDataZh, pdf_alpha);

  *f_alpha = new TF1(Form("f_alpha%i",i), "[0]*exp(-x/([1]+[2]*x))/exp(-x/([3]+[4]*x))", 750, 4300);

  (*f_alpha)->SetParameters(constant, a_domSg.getVal(), b_domSg.getVal(), a_domSb.getVal(), b_domSb.getVal());

  // Fit jet mass in data side band

  RooRealVar     lamda("lamda", "lamda", myVal.value("lamda"), myVal.value("lamdaMin"), myVal.value("lamdaMax"));
  RooExponential pdf_sbDataJet("pdf_sbDataJet", "pdf_sbDataJet", mJet, lamda);
  RooExtendPdf   ext_sbDataJet("ext_sbDataJet", "ext_sbDataJet", pdf_sbDataJet, nEv_sbData);
  RooFitResult*  res_sbDataJet = ext_sbDataJet.fitTo(set_sbData, SumW2Error(true), Extended(true), Range("lowSB,highSB"), Strategy(2), Minimizer("Minuit2"), Save(1));

  // Normalize factor to normalize the background in signal region of data

  RooAbsReal* nSIGFit = ext_sbDataJet.createIntegral(RooArgSet(mJet), NormSet(mJet), Range("signal"));
  RooAbsReal* nSBFit  = ext_sbDataJet.createIntegral(RooArgSet(mJet), NormSet(mJet), Range("lowSB,highSB"));

  RooRealVar normFactor("normFactor", "normFactor", 0., 1.e9);
  normFactor.setVal(nEv_sbData.getVal()*(nSIGFit->getVal()/nSBFit->getVal()));

  RooFormulaVar normFormula("normFormula", "normFormula", "@0*@1/@2", RooArgList(nEv_sbData, *nSIGFit, *nSBFit));

  fprintf(stdout, "a_domSb = %.3f +- %.3f\n", a_domSb.getVal(), a_domSb.getError());
  fprintf(stdout, "b_domSb = %.3f +- %.3f\n", b_domSb.getVal(), b_domSb.getError());
  fprintf(stdout, "a_domSg = %.3f +- %.3f\n", a_domSg.getVal(), a_domSg.getError());
  fprintf(stdout, "b_domSg = %.3f +- %.3f\n", b_domSg.getVal(), b_domSg.getError());
  fprintf(stdout, "a_subSb = %.3f +- %.3f\n", a_subSb.getVal(), a_subSb.getError());
  fprintf(stdout, "b_subSb = %.3f +- %.3f\n", b_subSb.getVal(), b_subSb.getError());
  fprintf(stdout, "a_subSg = %.3f +- %.3f\n", a_subSg.getVal(), a_subSg.getError());
  fprintf(stdout, "b_subSg = %.3f +- %.3f\n", b_subSg.getVal(), b_subSg.getError());
  fprintf(stdout, "a_datSb = %.3f +- %.3f\n", a_datSb.getVal(), a_datSb.getError());
  fprintf(stdout, "b_datSb = %.3f +- %.3f\n", b_datSb.getVal(), b_datSb.getError());
  fprintf(stdout, "lamda   = %.3f +- %.3f\n", lamda  .getVal(), lamda  .getError());

  *h_shape = pdf_predicted.createHistogram(Form("h_shape%i",i), mZH, Binning(bin_mZH), Extended(true));

}
