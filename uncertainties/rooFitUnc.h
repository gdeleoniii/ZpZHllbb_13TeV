R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/PdfDiagonalizer_cc.so)
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/readFitParam.h"
using namespace RooFit;

void rooFitUnc(string channel, string catcut, string region, TF1** f_alpha, TF1** f_predict, TH1** h_shape, int i){

  // Suppress all the INFO message

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooMsgService::instance().setSilentMode(true);
  gROOT->ProcessLine("gErrorIgnoreLevel=kWarning;");

  // Input files and sum all data/backgrounds

  TChain* tree_Data = new TChain("tree");
  TChain* tree_Dom  = new TChain("tree");
  TChain* tree_Sub1 = new TChain("tree");
  TChain* tree_Sub2 = new TChain("tree");

  // Data

  if( channel == "ele" ){

    tree_Data->Add("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/data/SingleElectron-Run2015D-v1_eleMiniTree.root");
    tree_Data->Add("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/data/SingleElectron-Run2015D-v4_eleMiniTree.root");

  }

  else if( channel == "mu" ){

    tree_Data->Add("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/data/SingleMuon-Run2015D-v1_muMiniTree.root");
    tree_Data->Add("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/data/SingleMuon-Run2015D-v4_muMiniTree.root");

  }

  // Dominant and subdominant background

  tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-100to200_13TeV_%s_%sMiniTree.root", region.data(), channel.data()));
  tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-200to400_13TeV_%s_%sMiniTree.root", region.data(), channel.data()));
  tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-400to600_13TeV_%s_%sMiniTree.root", region.data(), channel.data()));
  tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-600toInf_13TeV_%s_%sMiniTree.root", region.data(), channel.data()));

  tree_Sub1->Add(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/minor/WW_TuneCUETP8M1_13TeV_%sMiniTree.root", channel.data()));
  tree_Sub1->Add(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/minor/WZ_TuneCUETP8M1_13TeV_%sMiniTree.root", channel.data()));
  tree_Sub1->Add(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/minor/ZZ_TuneCUETP8M1_13TeV_%sMiniTree.root", channel.data()));
  tree_Sub1->Add(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/minor/TT_TuneCUETP8M1_13TeV_%sMiniTree.root", channel.data()));
  tree_Sub2->Add(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/minor/ZH_HToBB_ZToLL_M125_13TeV_%sMiniTree.root", channel.data()));

  // Define all the variables from the trees

  RooRealVar cat      ("cat", "", 0, 2);
  RooRealVar mJet     ("prmass", "M_{jet}", 30., 300., "GeV");
  RooRealVar mZH      ("mllbb", "M_{ZH}", 750., 4300., "GeV");
  RooRealVar evWeight ("evweight", "", 0., 1.e3);

  // Set the range in zh mass and in jet mass

  mZH .setRange("All",  750., 4300.);
  mJet.setRange("All",   30.,  300.);
  mJet.setRange("SB_l",  30.,   65.);
  mJet.setRange("SB_h", 135.,  300.);
  mJet.setRange("SG",   105.,  135.);

  TCut cut_bTag = Form("cat==%s", catcut.data());
  TCut cut_sb   = "prmass>30 && !(prmass>65 && prmass<135) && prmass<300";
  TCut cut_sg   = "prmass>105 && prmass<135";

  // Create a dataset from a tree to process an unbinned likelihood fitting

  RooDataSet set_sbDom ("set_sbDom",  "set_sbDom",  RooArgSet(cat, mJet, mZH, evWeight), Cut(cut_bTag && cut_sb), Import(*tree_Dom),  WeightVar(evWeight)); 
  RooDataSet set_sgDom ("set_sgDom",  "set_sgDom",  RooArgSet(cat, mJet, mZH, evWeight), Cut(cut_bTag && cut_sg), Import(*tree_Dom),  WeightVar(evWeight));
  RooDataSet set_sbSub1("set_sbSub1", "set_sbSub1", RooArgSet(cat, mJet, mZH, evWeight), Cut(cut_bTag && cut_sb), Import(*tree_Sub1), WeightVar(evWeight));
  RooDataSet set_sgSub1("set_sgSub1", "set_sgSub1", RooArgSet(cat, mJet, mZH, evWeight), Cut(cut_bTag && cut_sg), Import(*tree_Sub1), WeightVar(evWeight));
  RooDataSet set_sbSub2("set_sbSub2", "set_sbSub2", RooArgSet(cat, mJet, mZH, evWeight), Cut(cut_bTag && cut_sb), Import(*tree_Sub2), WeightVar(evWeight));
  RooDataSet set_sgSub2("set_sgSub2", "set_sgSub2", RooArgSet(cat, mJet, mZH, evWeight), Cut(cut_bTag && cut_sg), Import(*tree_Sub2), WeightVar(evWeight)); 
  RooDataSet set_sbData("set_sbData", "set_sbData", RooArgSet(cat, mJet, mZH, evWeight), Cut(cut_bTag && cut_sb), Import(*tree_Data), WeightVar(evWeight));
  RooDataSet set_sgData("set_sgData", "set_sgData", RooArgSet(cat, mJet, mZH, evWeight), Cut(cut_bTag && cut_sg), Import(*tree_Data), WeightVar(evWeight));

  // Total events number

  RooRealVar nEv_sbDom ("nEv_sbDom",  "nEv_sbDom",  0, 1e9);
  RooRealVar nEv_sgDom ("nEv_sgDom",  "nEv_sgDom",  0, 1e9);
  RooRealVar nEv_sbSub1("nEv_sbSub1", "nEv_sbSub1", 0, 1e9);
  RooRealVar nEv_sgSub1("nEv_sgSub1", "nEv_sgSub1", 0, 1e9);
  RooRealVar nEv_sbSub2("nEv_sbSub2", "nEv_sbSub2", 0, 1e9);
  RooRealVar nEv_sgSub2("nEv_sgSub2", "nEv_sgSub2", 0, 1e9);
  RooRealVar nEv_sbData("nEv_sbData", "nEv_sbData", 0, 1e9);

  nEv_sbDom .setVal(set_sbDom .sumEntries());
  nEv_sgDom .setVal(set_sbDom .sumEntries());
  nEv_sbSub1.setVal(set_sbSub1.sumEntries());
  nEv_sgSub1.setVal(set_sgSub1.sumEntries());
  nEv_sbSub2.setVal(set_sbSub2.sumEntries());
  nEv_sgSub2.setVal(set_sgSub2.sumEntries());
  nEv_sbData.setVal(set_sbData.sumEntries());

  // Set fit parameters for ZH mass

  param myVal(channel.data(), catcut.data());

  RooRealVar a_domSb ("a_domSb",  "a_domSb",  myVal.value("a_domSb"),  myVal.value("a_domSbMin"),  myVal.value("a_domSbMax"));
  RooRealVar b_domSb ("b_domSb",  "b_domSb",  myVal.value("b_domSb"),  myVal.value("b_domSbMin"),  myVal.value("b_domSbMax"));
  RooRealVar a_domSg ("a_domSg",  "a_domSg",  myVal.value("a_domSg"),  myVal.value("a_domSgMin"),  myVal.value("a_domSgMax"));
  RooRealVar b_domSg ("b_domSg",  "b_domSg",  myVal.value("b_domSg"),  myVal.value("b_domSgMin"),  myVal.value("b_domSgMax"));
  RooRealVar a_dataSb("a_dataSb", "a_dataSb", myVal.value("a_dataSb"), myVal.value("a_dataSbMin"), myVal.value("a_dataSbMax"));
  RooRealVar b_dataSb("b_dataSb", "b_dataSb", myVal.value("b_dataSb"), myVal.value("b_dataSbMin"), myVal.value("b_dataSbMax"));
  RooRealVar a_sub1Sb("a_sub1Sb", "a_sub1Sb", myVal.value("a_sub1Sb"), myVal.value("a_sub1SbMin"), myVal.value("a_sub1SbMax"));
  RooRealVar a_sub1Sg("a_sub1Sg", "a_sub1Sg", myVal.value("a_sub1Sg"), myVal.value("a_sub1SgMin"), myVal.value("a_sub1SgMax"));
  RooRealVar a_sub2Sb("a_sub2Sb", "a_sub2Sb", myVal.value("a_sub2Sb"), myVal.value("a_sub2SbMin"), myVal.value("a_sub2SbMax"));
  RooRealVar a_sub2Sg("a_sub2Sg", "a_sub2Sg", myVal.value("a_sub2Sg"), myVal.value("a_sub2SgMin"), myVal.value("a_sub2SgMax"));

  // Create pdf for ZH mass

  RooGenericPdf pdf_sbDomZh ("pdf_sbDomZh",  "pdf_sbDomZh",  "exp(-@0/(@1+@2*@0))", RooArgSet(mZH,a_domSb,b_domSb));
  RooGenericPdf pdf_sgDomZh ("pdf_sgDomZh",  "pdf_sgDomZh",  "exp(-@0/(@1+@2*@0))", RooArgSet(mZH,a_domSg,b_domSg));
  RooGenericPdf pdf_sbDataZh("pdf_sbDataZh", "pdf_sbDataZh", "exp(-@0/(@1+@2*@0))", RooArgSet(mZH,a_dataSb,b_dataSb));
  RooGenericPdf pdf_sbSub1Zh("pdf_sbSub1Zh", "pdf_sbSub1Zh", "exp(-@0/@1)",         RooArgSet(mZH,a_sub1Sb));
  RooGenericPdf pdf_sgSub1Zh("pdf_sgSub1Zh", "pdf_sgSub1Zh", "exp(-@0/@1)",         RooArgSet(mZH,a_sub1Sg));
  RooGenericPdf pdf_sbSub2Zh("pdf_sbSub2Zh", "pdf_sbSub2Zh", "exp(-@0/@1)",         RooArgSet(mZH,a_sub2Sb));
  RooGenericPdf pdf_sgSub2Zh("pdf_sgSub2Zh", "pdf_sgSub2Zh", "exp(-@0/@1)",         RooArgSet(mZH,a_sub2Sg));

  // Extended pdf from RooGenericPdf

  RooExtendPdf ext_sbDomZh ("ext_sbDomZh",  "ext_sbDomZh",  pdf_sbDomZh,  nEv_sbDom);
  RooExtendPdf ext_sgDomZh ("ext_sgDomZh",  "ext_sgDomZh",  pdf_sgDomZh,  nEv_sgDom);
  RooExtendPdf ext_sbSub1Zh("ext_sbSub1Zh", "ext_sbSub1Zh", pdf_sbSub1Zh, nEv_sbSub1);
  RooExtendPdf ext_sgSub1Zh("ext_sgSub1Zh", "ext_sgSub1Zh", pdf_sgSub1Zh, nEv_sgSub1);
  RooExtendPdf ext_sbSub2Zh("ext_sbSub2Zh", "ext_sbSub2Zh", pdf_sbSub2Zh, nEv_sbSub2);
  RooExtendPdf ext_sgSub2Zh("ext_sgSub2Zh", "ext_sgSub2Zh", pdf_sgSub2Zh, nEv_sgSub2);
  RooExtendPdf ext_sbDataZh("ext_sbDataZh", "ext_sbDataZh", pdf_sbDataZh, nEv_sbData);

  // Make category to fit dominant/subdominant background signal/sideband and data sideband

  RooCategory cat_combine("cat_combine", "cat_combine");

  cat_combine.defineType("dom_SB");
  cat_combine.defineType("dom_SG");
  cat_combine.defineType("sub1_SB");
  cat_combine.defineType("sub1_SG");
  cat_combine.defineType("sub2_SB");
  cat_combine.defineType("sub2_SG");
  cat_combine.defineType("data_SB");

  map<string, RooDataSet*> mapSet;

  mapSet.insert(pair<string, RooDataSet*>("dom_SB",  &set_sbDom));
  mapSet.insert(pair<string, RooDataSet*>("dom_SG",  &set_sgDom));
  mapSet.insert(pair<string, RooDataSet*>("sub1_SB", &set_sbSub1));
  mapSet.insert(pair<string, RooDataSet*>("sub1_SG", &set_sgSub1));
  mapSet.insert(pair<string, RooDataSet*>("sub2_SB", &set_sbSub2));
  mapSet.insert(pair<string, RooDataSet*>("sub2_SG", &set_sgSub2));
  mapSet.insert(pair<string, RooDataSet*>("data_SB", &set_sbData));

  RooDataSet cmb_combine("cmb_combine", "cmb_combine", RooArgSet(cat, mJet, mZH, evWeight), Index(cat_combine), Import(mapSet), WeightVar(evWeight));

  RooSimultaneous pdf_combine("pdf_combine", "pdf_combine", cat_combine);

  pdf_combine.addPdf(ext_sbDomZh,  "dom_SB");
  pdf_combine.addPdf(ext_sgDomZh,  "dom_SG");
  pdf_combine.addPdf(ext_sbSub1Zh, "sub1_SB");
  pdf_combine.addPdf(ext_sgSub1Zh, "sub1_SG");
  pdf_combine.addPdf(ext_sbSub2Zh, "sub2_SB");
  pdf_combine.addPdf(ext_sgSub2Zh, "sub2_SG");
  pdf_combine.addPdf(ext_sbDataZh, "data_SB");

  pdf_combine.fitTo(cmb_combine, SumW2Error(true), Extended(true), Range("All"), Strategy(2), Minimizer("Minuit2"), Save(1));

  // Multiply the model of background in data side band with the model of alpha ratio to the a model of background in data signal region
  // predicted background = (sbDataZh - sbSub1Zh - sbSub2Zh) * alpha + sgSub1Zh + sgSub2Zh

  float alpConst = ext_sbDomZh.createIntegral(mZH)->getVal()/ext_sgDomZh.createIntegral(mZH)->getVal();

  *f_alpha   = new TF1(Form("f_alpha%i",i), "[0]*exp(-x/([1]+[2]*x))/exp(-x/([3]+[4]*x))", 750, 4300);
  *f_predict = new TF1(Form("f_predict%i",i), "([0]*exp(-x/([1]+[2]*x))-[3]*exp(-x/[4])-[5]*exp(-x/[6]))*[7]*exp(-x/([8]+[9]*x))/exp(-x/([10]+[11]*x))+[12]*exp(-x/[13])+[14]*exp(-x/[15])", 750, 4300);

  double param_alpha[5]    = {alpConst, a_domSg.getVal(), b_domSg.getVal(), a_domSb.getVal(), b_domSb.getVal()};
  double param_predict[16] = {nEv_sbData.getVal(), a_dataSb.getVal(), b_dataSb.getVal(), nEv_sbSub1.getVal(), a_sub1Sb.getVal(), nEv_sbSub2.getVal(), a_sub2Sb.getVal(), alpConst, a_domSg.getVal(), b_domSg.getVal(), a_domSb.getVal(), b_domSb.getVal(), nEv_sgSub1.getVal(), a_sub1Sg.getVal(), nEv_sgSub2.getVal(), a_sub2Sg.getVal()};

  (*f_alpha)  ->SetParameters(param_alpha);
  (*f_predict)->SetParameters(param_predict);

  // jet mass in data side band

  RooRealVar    j_data("j_data", "j_data", myVal.value("j_data"), myVal.value("j_dataMin"), myVal.value("j_dataMax"));
  RooGenericPdf pdf_dataJet("pdf_dataJet", "pdf_dataJet", "exp(-@0/@1)", RooArgSet(mJet, j_data));
  RooExtendPdf  ext_dataJet("ext_dataJet", "ext_dataJet", pdf_dataJet, nEv_sbData);

  ext_dataJet.fitTo(set_sbData, SumW2Error(true), Extended(true), Range("SB_l,SB_h"), Strategy(2), Minimizer("Minuit2"), Save(1));

  // Normalize factor to normalize the background in signal region of data

  RooAbsReal* nFit_sg = ext_dataJet.createIntegral(RooArgSet(mJet), NormSet(mJet), Range("SG"));
  RooAbsReal* nFit_sb = ext_dataJet.createIntegral(RooArgSet(mJet), NormSet(mJet), Range("SB_l,SB_h"));

  // Since the statistic of 2015 data is low, the jet mass distribution in 2 btag is consider as a flat distribution

  float normFactorVal = (catcut=="1") ? nEv_sbData.getVal()*(nFit_sg->getVal()/nFit_sb->getVal()) : 6;

  *h_shape = (*f_predict)->CreateHistogram();
  (*h_shape)->Scale(normFactorVal/(*h_shape)->Integral());
  (*h_shape)->SetBins(71, 750, 4300);

  fprintf(stdout, "nw = %i\t", i);
  fprintf(stdout, "a_domSb = %.3f +- %.3f\t", a_domSb.getVal(), a_domSb.getError());
  fprintf(stdout, "b_domSb = %.3f +- %.3f\t", b_domSb.getVal(), b_domSb.getError());
  fprintf(stdout, "a_domSg = %.3f +- %.3f\t", a_domSg.getVal(), a_domSg.getError());
  fprintf(stdout, "b_domSg = %.3f +- %.3f\n", b_domSg.getVal(), b_domSg.getError());

  fprintf(stdout, "real norm = %.3f\t", normFactorVal);
  fprintf(stdout, "hist norm = %.3f\n", (*h_shape)->Integral());

}