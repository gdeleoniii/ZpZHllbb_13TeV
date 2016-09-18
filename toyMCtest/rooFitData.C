R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/PdfDiagonalizer_cc.so)
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/readFitParam.h"
using namespace RooFit;

void rooFitData(string channel, string catcut){

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

    tree_Data->Add("data/SingleElectron-Run2015D-v1_eleMiniTree.root");
    tree_Data->Add("data/SingleElectron-Run2015D-v4_eleMiniTree.root");

  }

  else if( channel == "mu" ){

    tree_Data->Add("data/SingleMuon-Run2015D-v1_muMiniTree.root");
    tree_Data->Add("data/SingleMuon-Run2015D-v4_muMiniTree.root");

  }

  // Dominant and subdominant background

  tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-100to200_13TeV_%sMiniTree.root", channel.data()));
  tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-200to400_13TeV_%sMiniTree.root", channel.data()));
  tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-400to600_13TeV_%sMiniTree.root", channel.data()));
  tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-600toInf_13TeV_%sMiniTree.root", channel.data()));

  tree_Sub1->Add(Form("minor/WW_TuneCUETP8M1_13TeV_%sMiniTree.root", channel.data()));
  tree_Sub1->Add(Form("minor/WZ_TuneCUETP8M1_13TeV_%sMiniTree.root", channel.data()));
  tree_Sub1->Add(Form("minor/ZZ_TuneCUETP8M1_13TeV_%sMiniTree.root", channel.data()));
  tree_Sub1->Add(Form("minor/TT_TuneCUETP8M1_13TeV_%sMiniTree.root", channel.data()));
  tree_Sub2->Add(Form("minor/ZH_HToBB_ZToLL_M125_13TeV_%sMiniTree.root", channel.data()));

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

  RooBinning bin_mZH (71, 750, 4300);
  RooBinning bin_mJet(54, 30, 300);

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

  RooDataSet set_combine("set_combine", "set_combine", RooArgSet(cat, mJet, mZH, evWeight), Index(cat_combine), Import(mapSet), WeightVar(evWeight));

  RooSimultaneous pdf_combine("pdf_combine", "pdf_combine", cat_combine);

  pdf_combine.addPdf(ext_sbDomZh,  "dom_SB");
  pdf_combine.addPdf(ext_sgDomZh,  "dom_SG");
  pdf_combine.addPdf(ext_sbSub1Zh, "sub1_SB");
  pdf_combine.addPdf(ext_sgSub1Zh, "sub1_SG");
  pdf_combine.addPdf(ext_sbSub2Zh, "sub2_SB");
  pdf_combine.addPdf(ext_sgSub2Zh, "sub2_SG");
  pdf_combine.addPdf(ext_sbDataZh, "data_SB");

  RooFitResult* res_combine = pdf_combine.fitTo(set_combine, SumW2Error(true), Extended(true), Range("All"), Strategy(2), Minimizer("Minuit2"), Save(1));

  // Multiply the model of background in data side band with the model of alpha ratio to the a model of background in data signal region
  // predicted background = (sbDataZh - sbSub1Zh - sbSub2Zh) * alpha + sgSub1Zh + sgSub2Zh

  RooArgSet arg_combine(mZH);

  arg_combine.add(a_dataSb);
  arg_combine.add(b_dataSb);
  arg_combine.add(a_sub1Sb);
  arg_combine.add(a_sub2Sb);
  arg_combine.add(a_domSg);
  arg_combine.add(b_domSg);
  arg_combine.add(a_domSb);
  arg_combine.add(b_domSb);
  arg_combine.add(a_sub1Sg);
  arg_combine.add(a_sub2Sg);

  float alpConst = ext_sbDomZh.createIntegral(mZH)->getVal()/ext_sgDomZh.createIntegral(mZH)->getVal();

  RooGenericPdf pdf_alpha("pdf_alpha", "pdf_alpha", Form("%f*exp(-@0/(@1+@2*@0))/exp(-@0/(@3+@4*@0))", alpConst), RooArgSet(mZH,a_domSg,b_domSg,a_domSb,b_domSb));
  RooGenericPdf pdf_predict("pdf_predict", "pdf_predict", Form("(%f*exp(-@0/(@1+@2*@0))-%f*exp(-@0/@3)-%f*exp(-@0/@4))*%f*exp(-@0/(@5+@6*@0))/exp(-@0/(@7+@8*@0))+%f*exp(-@0/@9)+%f*exp(-@0/@10)", nEv_sbData.getVal(), nEv_sbSub1.getVal(), nEv_sbSub2.getVal(), alpConst, nEv_sgSub1.getVal(), nEv_sgSub2.getVal()), arg_combine);

  // jet mass in data side band

  RooRealVar    j_data("j_data", "j_data", myVal.value("j_data"), myVal.value("j_dataMin"), myVal.value("j_dataMax"));
  RooGenericPdf pdf_dataJet("pdf_dataJet", "pdf_dataJet", "exp(-@0/@1)", RooArgSet(mJet, j_data));
  RooExtendPdf  ext_dataJet("ext_dataJet", "ext_dataJet", pdf_dataJet, nEv_sbData);

  RooFitResult* res_dataJet = ext_dataJet.fitTo(set_sbData, SumW2Error(true), Extended(true), Range("SB_l,SB_h"), Strategy(2), Minimizer("Minuit2"), Save(1));

  // Normalize factor to normalize the background in signal region of data

  RooAbsReal* nFit_sg = ext_dataJet.createIntegral(RooArgSet(mJet), NormSet(mJet), Range("SG"));
  RooAbsReal* nFit_sb = ext_dataJet.createIntegral(RooArgSet(mJet), NormSet(mJet), Range("SB_l,SB_h"));

  // Since the statistic of 2015 data is low, the jet mass distribution in 2 btag is consider as a flat distribution

  RooFormulaVar normFormula("normFormula", "normFormula", "@0*@1/@2", RooArgList(nEv_sbData, *nFit_sg, *nFit_sb));

  float normFactorVal = nEv_sbData.getVal()*(nFit_sg->getVal()/nFit_sb->getVal());
  float normFactorUnc = (catcut=="1") ? normFormula.getPropagatedError(*res_dataJet) : fabs(6-normFactorVal);

  RooRealVar normFactor("normFactor", "normFactor", 0, 1e9);
  normFactor.setVal( (catcut=="1") ? normFactorVal : 6 );

  fprintf(stdout, "nEv_sbDom  = %.3f +- %.3f\n", nEv_sbDom .getVal(), nEv_sbDom .getError());
  fprintf(stdout, "nEv_sgDom  = %.3f +- %.3f\n", nEv_sgDom .getVal(), nEv_sgDom .getError());
  fprintf(stdout, "nEv_sbSub1 = %.3f +- %.3f\n", nEv_sbSub1.getVal(), nEv_sbSub1.getError());
  fprintf(stdout, "nEv_sgSub1 = %.3f +- %.3f\n", nEv_sgSub1.getVal(), nEv_sgSub1.getError());
  fprintf(stdout, "nEv_sbSub2 = %.3f +- %.3f\n", nEv_sbSub2.getVal(), nEv_sbSub2.getError());
  fprintf(stdout, "nEv_sgSub2 = %.3f +- %.3f\n", nEv_sgSub2.getVal(), nEv_sgSub2.getError());
  fprintf(stdout, "nEv_sbData = %.3f +- %.3f\n", nEv_sbData.getVal(), nEv_sbData.getError());
  fprintf(stdout, "a_domSb    = %.3f +- %.3f\n", a_domSb   .getVal(), a_domSb   .getError());
  fprintf(stdout, "b_domSb    = %.3f +- %.3f\n", b_domSb   .getVal(), b_domSb   .getError());
  fprintf(stdout, "a_domSg    = %.3f +- %.3f\n", a_domSg   .getVal(), a_domSg   .getError());
  fprintf(stdout, "b_domSg    = %.3f +- %.3f\n", b_domSg   .getVal(), b_domSg   .getError());
  fprintf(stdout, "a_dataSb   = %.3f +- %.3f\n", a_dataSb  .getVal(), a_dataSb  .getError());
  fprintf(stdout, "b_dataSb   = %.3f +- %.3f\n", b_dataSb  .getVal(), b_dataSb  .getError());
  fprintf(stdout, "a_sub1Sb   = %.3f +- %.3f\n", a_sub1Sb  .getVal(), a_sub1Sb  .getError());
  fprintf(stdout, "a_sub1Sg   = %.3f +- %.3f\n", a_sub1Sg  .getVal(), a_sub1Sg  .getError());
  fprintf(stdout, "a_sub2Sb   = %.3f +- %.3f\n", a_sub2Sb  .getVal(), a_sub2Sb  .getError());
  fprintf(stdout, "a_sub2Sg   = %.3f +- %.3f\n", a_sub2Sg  .getVal(), a_sub2Sg  .getError());
  fprintf(stdout, "j_data     = %.3f +- %.3f\n", j_data    .getVal(), j_data    .getError());

  // Plot the results on frame 

  RooPlot* frm_sbDomZh        = mZH.frame();
  RooPlot* frm_sgDomZh        = mZH.frame();
  RooPlot* frm_sbSub1Zh       = mZH.frame();
  RooPlot* frm_sgSub1Zh       = mZH.frame();
  RooPlot* frm_sbSub2Zh       = mZH.frame();
  RooPlot* frm_sgSub2Zh       = mZH.frame();
  RooPlot* frm_sbDataZh       = mZH.frame();
  RooPlot* frm_alpha          = mZH.frame();
  RooPlot* frm_predict        = mZH.frame();
  RooPlot* frm_dataJet        = mJet.frame();
  RooPlot* frm_sbDomZh_pull   = mZH.frame();
  RooPlot* frm_sgDomZh_pull   = mZH.frame();
  RooPlot* frm_sbSub1Zh_pull  = mZH.frame();
  RooPlot* frm_sgSub1Zh_pull  = mZH.frame();
  RooPlot* frm_sbSub2Zh_pull  = mZH.frame();
  RooPlot* frm_sgSub2Zh_pull  = mZH.frame();
  RooPlot* frm_sbDataZh_pull  = mZH.frame();
  RooPlot* frm_dataJet_pull   = mJet.frame();

  set_combine.plotOn(frm_sbDomZh, Cut("cat_combine==cat_combine::dom_SB"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_combine.plotOn(frm_sbDomZh, Slice(cat_combine,"dom_SB"), ProjWData(cat_combine,set_combine), VisualizeError(*res_combine,1,false), FillStyle(3002));
  set_combine.plotOn(frm_sbDomZh, Cut("cat_combine==cat_combine::dom_SB"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_combine.plotOn(frm_sbDomZh, Slice(cat_combine,"dom_SB"), ProjWData(cat_combine,set_combine), LineColor(kBlue));

  set_combine.plotOn(frm_sgDomZh, Cut("cat_combine==cat_combine::dom_SG"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_combine.plotOn(frm_sgDomZh, Slice(cat_combine,"dom_SG"), ProjWData(cat_combine,set_combine), VisualizeError(*res_combine,1,false), FillStyle(3002));
  set_combine.plotOn(frm_sgDomZh, Cut("cat_combine==cat_combine::dom_SG"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_combine.plotOn(frm_sgDomZh, Slice(cat_combine,"dom_SG"), ProjWData(cat_combine,set_combine), LineColor(kBlue));

  set_combine.plotOn(frm_sbSub1Zh, Cut("cat_combine==cat_combine::sub1_SB"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_combine.plotOn(frm_sbSub1Zh, Slice(cat_combine,"sub1_SB"), ProjWData(cat_combine,set_combine), VisualizeError(*res_combine,1,false), FillStyle(3002));
  set_combine.plotOn(frm_sbSub1Zh, Cut("cat_combine==cat_combine::sub1_SB"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_combine.plotOn(frm_sbSub1Zh, Slice(cat_combine,"sub1_SB"), ProjWData(cat_combine,set_combine), LineColor(kBlue));

  set_combine.plotOn(frm_sgSub1Zh, Cut("cat_combine==cat_combine::sub1_SG"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_combine.plotOn(frm_sgSub1Zh, Slice(cat_combine,"sub1_SG"), ProjWData(cat_combine,set_combine), VisualizeError(*res_combine,1,false), FillStyle(3002));
  set_combine.plotOn(frm_sgSub1Zh, Cut("cat_combine==cat_combine::sub1_SG"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_combine.plotOn(frm_sgSub1Zh, Slice(cat_combine,"sub1_SG"), ProjWData(cat_combine,set_combine), LineColor(kBlue));

  set_combine.plotOn(frm_sbSub2Zh, Cut("cat_combine==cat_combine::sub2_SB"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_combine.plotOn(frm_sbSub2Zh, Slice(cat_combine,"sub2_SB"), ProjWData(cat_combine,set_combine), VisualizeError(*res_combine,1,false), FillStyle(3002));
  set_combine.plotOn(frm_sbSub2Zh, Cut("cat_combine==cat_combine::sub2_SB"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_combine.plotOn(frm_sbSub2Zh, Slice(cat_combine,"sub2_SB"), ProjWData(cat_combine,set_combine), LineColor(kBlue));

  set_combine.plotOn(frm_sgSub2Zh, Cut("cat_combine==cat_combine::sub2_SG"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_combine.plotOn(frm_sgSub2Zh, Slice(cat_combine,"sub2_SG"), ProjWData(cat_combine,set_combine), VisualizeError(*res_combine,1,false), FillStyle(3002));
  set_combine.plotOn(frm_sgSub2Zh, Cut("cat_combine==cat_combine::sub2_SG"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_combine.plotOn(frm_sgSub2Zh, Slice(cat_combine,"sub2_SG"), ProjWData(cat_combine,set_combine), LineColor(kBlue));

  set_combine.plotOn(frm_sbDataZh, Cut("cat_combine==cat_combine::data_SB"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_combine.plotOn(frm_sbDataZh, Slice(cat_combine,"data_SB"), ProjWData(cat_combine,set_combine), VisualizeError(*res_combine,1,false), FillStyle(3002));
  set_combine.plotOn(frm_sbDataZh, Cut("cat_combine==cat_combine::data_SB"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_combine.plotOn(frm_sbDataZh, Slice(cat_combine,"data_SB"), ProjWData(cat_combine,set_combine), LineColor(kBlue));

  pdf_alpha  .plotOn(frm_alpha, VisualizeError(*res_combine,1,false), FillStyle(3002), FillColor(kBlack));
  pdf_alpha  .plotOn(frm_alpha, LineColor(kBlack));
  ext_sbDomZh.plotOn(frm_alpha, Normalization(1, RooAbsReal::NumEvent), LineColor(kBlue));
  ext_sgDomZh.plotOn(frm_alpha, Normalization(1, RooAbsReal::NumEvent), LineColor(kRed));

  set_sgData .plotOn(frm_predict, DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_predict.plotOn(frm_predict, VisualizeError(*res_combine,1,false), Normalization(normFactor.getVal(), RooAbsReal::NumEvent), FillStyle(3002));
  set_sgData .plotOn(frm_predict, DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_predict.plotOn(frm_predict, Normalization(normFactor.getVal(), RooAbsReal::NumEvent), LineColor(kBlue));
  // Using RooAbsReal::NumEvent in order to consider the bin width of data set. Equivalent to (normFactor*binWidth) if using RooAbsReal::Raw.

  set_sbData .plotOn(frm_dataJet, DataError(RooAbsData::SumW2), Binning(bin_mJet));
  if( catcut == "1" ){
    ext_dataJet.plotOn(frm_dataJet, Range("All"), VisualizeError(*res_dataJet,1,false), FillStyle(3002));
    set_sbData .plotOn(frm_dataJet, DataError(RooAbsData::SumW2), Binning(bin_mJet));
    ext_dataJet.plotOn(frm_dataJet, Range("All"));
  }

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

  frm_sbDomZh->SetTitle("");
  frm_sbDomZh->SetMinimum(1e-4);
  frm_sbDomZh->SetMaximum(catcut=="1"?1000:100);
  frm_sbDomZh->GetXaxis()->SetTitle("");
  frm_sbDomZh->GetXaxis()->SetLabelOffset(999);
  frm_sbDomZh->Draw();

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.65, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s b-tag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, "Z+jets in sidebands");

  c0_up->RedrawAxis();
  c0_dw->cd()->SetLogy(0);

  frm_sbDomZh_pull->addObject(frm_sbDomZh->pullHist(), "P");
  frm_sbDomZh_pull->SetTitle("");
  frm_sbDomZh_pull->GetYaxis()->SetTitle("Pulls");
  frm_sbDomZh_pull->GetYaxis()->SetTitleOffset(0.25);
  frm_sbDomZh_pull->GetXaxis()->SetLabelSize(0.125);
  frm_sbDomZh_pull->GetXaxis()->SetTitleSize(0.125);
  frm_sbDomZh_pull->GetYaxis()->SetLabelSize(0.125);
  frm_sbDomZh_pull->GetYaxis()->SetTitleSize(0.125);
  frm_sbDomZh_pull->GetYaxis()->SetNdivisions(505);
  frm_sbDomZh_pull->SetMinimum(-4);
  frm_sbDomZh_pull->SetMaximum(4);
  frm_sbDomZh_pull->Draw();

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

  frm_sgDomZh->SetTitle("");
  frm_sgDomZh->SetMinimum(1e-4);
  frm_sgDomZh->SetMaximum(catcut=="1"?100:10);
  frm_sgDomZh->GetXaxis()->SetTitle("");
  frm_sgDomZh->GetXaxis()->SetLabelOffset(999);
  frm_sgDomZh->Draw();

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.65, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s b-tag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, "Z+jets in signal region");

  c1_up->RedrawAxis();
  c1_dw->cd()->SetLogy(0);

  frm_sgDomZh_pull->addObject(frm_sgDomZh->pullHist(), "P");
  frm_sgDomZh_pull->SetTitle("");
  frm_sgDomZh_pull->GetYaxis()->SetTitle("Pulls");
  frm_sgDomZh_pull->GetYaxis()->SetTitleOffset(0.25);
  frm_sgDomZh_pull->GetXaxis()->SetLabelSize(0.125);
  frm_sgDomZh_pull->GetXaxis()->SetTitleSize(0.125);
  frm_sgDomZh_pull->GetYaxis()->SetLabelSize(0.125);
  frm_sgDomZh_pull->GetYaxis()->SetTitleSize(0.125);
  frm_sgDomZh_pull->GetYaxis()->SetNdivisions(505);
  frm_sgDomZh_pull->SetMinimum(-4);
  frm_sgDomZh_pull->SetMaximum(4);
  frm_sgDomZh_pull->Draw();

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

  frm_sbDataZh->SetTitle("");
  frm_sbDataZh->SetMinimum(1e-1);
  frm_sbDataZh->SetMaximum(catcut=="1"?100:10);
  frm_sbDataZh->GetXaxis()->SetTitle("");
  frm_sbDataZh->GetXaxis()->SetLabelOffset(999);
  frm_sbDataZh->Draw();

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar.DrawLatexNDC(0.65, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s btag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, "data in sidebands");

  c2_up->RedrawAxis();
  c2_dw->cd()->SetLogy(0);

  frm_sbDataZh_pull->addObject(frm_sbDataZh->pullHist(), "P");
  frm_sbDataZh_pull->SetTitle("");
  frm_sbDataZh_pull->GetYaxis()->SetTitle("Pulls");
  frm_sbDataZh_pull->GetYaxis()->SetTitleOffset(0.25);
  frm_sbDataZh_pull->GetXaxis()->SetLabelSize(0.125);
  frm_sbDataZh_pull->GetXaxis()->SetTitleSize(0.125);
  frm_sbDataZh_pull->GetYaxis()->SetLabelSize(0.125);
  frm_sbDataZh_pull->GetYaxis()->SetTitleSize(0.125);
  frm_sbDataZh_pull->GetYaxis()->SetNdivisions(505);
  frm_sbDataZh_pull->SetMinimum(-4);
  frm_sbDataZh_pull->SetMaximum(4);
  frm_sbDataZh_pull->Draw();

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

  frm_sbSub1Zh->SetTitle("");
  frm_sbSub1Zh->SetMinimum(1e-4);
  frm_sbSub1Zh->SetMaximum(10);
  frm_sbSub1Zh->GetXaxis()->SetTitle("");
  frm_sbSub1Zh->GetXaxis()->SetLabelOffset(999);
  frm_sbSub1Zh->Draw();

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.60, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s btag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, "Subdominant background (VV+t#bar{t}) in sidebands");
  
  c3_up->RedrawAxis();
  c3_dw->cd()->SetLogy(0);

  frm_sbSub1Zh_pull->addObject(frm_sbSub1Zh->pullHist(), "P");
  frm_sbSub1Zh_pull->SetTitle("");
  frm_sbSub1Zh_pull->GetYaxis()->SetTitle("Pulls");
  frm_sbSub1Zh_pull->GetYaxis()->SetTitleOffset(0.25);
  frm_sbSub1Zh_pull->GetXaxis()->SetLabelSize(0.125);
  frm_sbSub1Zh_pull->GetXaxis()->SetTitleSize(0.125);
  frm_sbSub1Zh_pull->GetYaxis()->SetLabelSize(0.125);
  frm_sbSub1Zh_pull->GetYaxis()->SetTitleSize(0.125);
  frm_sbSub1Zh_pull->GetYaxis()->SetNdivisions(505);
  frm_sbSub1Zh_pull->SetMinimum(-4);
  frm_sbSub1Zh_pull->SetMaximum(4);
  frm_sbSub1Zh_pull->Draw();

  c3.Draw();
  c3.Print(Form("rooFit_forData_%s_cat%s.pdf", channel.data(), catcut.data()));

  TCanvas c4("c4","",0,0,1000,800);

  c4.Divide(1,2);

  TPad* c4_up = (TPad*)c4.GetListOfPrimitives()->FindObject("c4_1");
  TPad* c4_dw = (TPad*)c4.GetListOfPrimitives()->FindObject("c4_2");

  c4_up->SetPad(0,1-up_height,1,1);
  c4_dw->SetPad(0,0,1,dw_height);
  c4_dw->SetBottomMargin(0.25);
  c4_up->cd()->SetLogy(1);

  frm_sgSub1Zh->SetTitle("");
  frm_sgSub1Zh->SetMinimum(1e-4);
  frm_sgSub1Zh->SetMaximum(10);
  frm_sgSub1Zh->GetXaxis()->SetTitle("");
  frm_sgSub1Zh->GetXaxis()->SetLabelOffset(999);
  frm_sgSub1Zh->Draw();

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.65, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s b-tag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, "Subdominant (VV+t#bar{t}) in signal region");

  c4_up->RedrawAxis();
  c4_dw->cd()->SetLogy(0);

  frm_sgSub1Zh_pull->addObject(frm_sgSub1Zh->pullHist(), "P");
  frm_sgSub1Zh_pull->SetTitle("");
  frm_sgSub1Zh_pull->GetYaxis()->SetTitle("Pulls");
  frm_sgSub1Zh_pull->GetYaxis()->SetTitleOffset(0.25);
  frm_sgSub1Zh_pull->GetXaxis()->SetLabelSize(0.125);
  frm_sgSub1Zh_pull->GetXaxis()->SetTitleSize(0.125);
  frm_sgSub1Zh_pull->GetYaxis()->SetLabelSize(0.125);
  frm_sgSub1Zh_pull->GetYaxis()->SetTitleSize(0.125);
  frm_sgSub1Zh_pull->GetYaxis()->SetNdivisions(505);
  frm_sgSub1Zh_pull->SetMinimum(-4);
  frm_sgSub1Zh_pull->SetMaximum(4);
  frm_sgSub1Zh_pull->Draw();

  c4.Draw();
  c4.Print(Form("rooFit_forData_%s_cat%s.pdf", channel.data(), catcut.data()));  

  TCanvas c5("c5","",0,0,1000,800);

  c5.Divide(1,2);

  TPad* c5_up = (TPad*)c5.GetListOfPrimitives()->FindObject("c5_1");
  TPad* c5_dw = (TPad*)c5.GetListOfPrimitives()->FindObject("c5_2");

  c5_up->SetPad(0,1-up_height,1,1);
  c5_dw->SetPad(0,0,1,dw_height);
  c5_dw->SetBottomMargin(0.25);
  c5_up->cd()->SetLogy(1);

  frm_sbSub2Zh->SetTitle("");
  frm_sbSub2Zh->SetMinimum(1e-4);
  frm_sbSub2Zh->SetMaximum(1);
  frm_sbSub2Zh->GetXaxis()->SetTitle("");
  frm_sbSub2Zh->GetXaxis()->SetLabelOffset(999);
  frm_sbSub2Zh->Draw();

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.60, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s btag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, "Subdominant background (ZH) in sidebands");

  c5_up->RedrawAxis();
  c5_dw->cd()->SetLogy(0);

  frm_sbSub2Zh_pull->addObject(frm_sbSub2Zh->pullHist(), "P");
  frm_sbSub2Zh_pull->SetTitle("");
  frm_sbSub2Zh_pull->GetYaxis()->SetTitle("Pulls");
  frm_sbSub2Zh_pull->GetYaxis()->SetTitleOffset(0.25);
  frm_sbSub2Zh_pull->GetXaxis()->SetLabelSize(0.125);
  frm_sbSub2Zh_pull->GetXaxis()->SetTitleSize(0.125);
  frm_sbSub2Zh_pull->GetYaxis()->SetLabelSize(0.125);
  frm_sbSub2Zh_pull->GetYaxis()->SetTitleSize(0.125);
  frm_sbSub2Zh_pull->GetYaxis()->SetNdivisions(505);
  frm_sbSub2Zh_pull->SetMinimum(-4);
  frm_sbSub2Zh_pull->SetMaximum(4);
  frm_sbSub2Zh_pull->Draw();

  c5.Draw();
  c5.Print(Form("rooFit_forData_%s_cat%s.pdf", channel.data(), catcut.data()));

  TCanvas c6("c6","",0,0,1000,800);

  c6.Divide(1,2);

  TPad* c6_up = (TPad*)c6.GetListOfPrimitives()->FindObject("c6_1");
  TPad* c6_dw = (TPad*)c6.GetListOfPrimitives()->FindObject("c6_2");

  c6_up->SetPad(0,1-up_height,1,1);
  c6_dw->SetPad(0,0,1,dw_height);
  c6_dw->SetBottomMargin(0.25);
  c6_up->cd()->SetLogy(1);

  frm_sgSub2Zh->SetTitle("");
  frm_sgSub2Zh->SetMinimum(1e-4);
  frm_sgSub2Zh->SetMaximum(1);
  frm_sgSub2Zh->GetXaxis()->SetTitle("");
  frm_sgSub2Zh->GetXaxis()->SetLabelOffset(999);
  frm_sgSub2Zh->Draw();

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.65, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s b-tag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, "Subdominant in (ZH) signal region");

  c6_up->RedrawAxis();
  c6_dw->cd()->SetLogy(0);

  frm_sgSub2Zh_pull->addObject(frm_sgSub2Zh->pullHist(), "P");
  frm_sgSub2Zh_pull->SetTitle("");
  frm_sgSub2Zh_pull->GetYaxis()->SetTitle("Pulls");
  frm_sgSub2Zh_pull->GetYaxis()->SetTitleOffset(0.25);
  frm_sgSub2Zh_pull->GetXaxis()->SetLabelSize(0.125);
  frm_sgSub2Zh_pull->GetXaxis()->SetTitleSize(0.125);
  frm_sgSub2Zh_pull->GetYaxis()->SetLabelSize(0.125);
  frm_sgSub2Zh_pull->GetYaxis()->SetTitleSize(0.125);
  frm_sgSub2Zh_pull->GetYaxis()->SetNdivisions(505);
  frm_sgSub2Zh_pull->SetMinimum(-4);
  frm_sgSub2Zh_pull->SetMaximum(4);
  frm_sgSub2Zh_pull->Draw();

  c6.Draw();
  c6.Print(Form("rooFit_forData_%s_cat%s.pdf", channel.data(), catcut.data()));

  TCanvas c7("c7","",0,0,1000,800);
  
  c7.Divide(1,2);

  TPad* c7_up = (TPad*)c7.GetListOfPrimitives()->FindObject("c7_1");
  TPad* c7_dw = (TPad*)c7.GetListOfPrimitives()->FindObject("c7_2"); 

  c7_up->SetPad(0,1-up_height,1,1);
  c7_dw->SetPad(0,0,1,dw_height);
  c7_dw->SetBottomMargin(0.25);
  c7_up->cd()->SetLogy(1);

  frm_dataJet->SetTitle("");
  frm_dataJet->SetMinimum(1e-1);
  frm_dataJet->SetMaximum(catcut=="1"?1e3:1e2);
  frm_dataJet->GetXaxis()->SetTitle("");
  frm_dataJet->GetXaxis()->SetLabelOffset(999);
  frm_dataJet->Draw();

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar.DrawLatexNDC(0.60, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s btag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, "data jet mass in sidebands");
  lar.DrawLatexNDC(0.15, 0.78, Form("Normalization factor: %.3f#pm%.3f", normFactor.getVal(), normFactorUnc));
  
  c7_up->RedrawAxis();
  c7_dw->cd()->SetLogy(0);

  frm_dataJet_pull->addObject(frm_dataJet->pullHist(), "P");
  frm_dataJet_pull->SetTitle("");
  frm_dataJet_pull->GetYaxis()->SetTitle("Pulls");
  frm_dataJet_pull->GetYaxis()->SetTitleOffset(0.25);
  frm_dataJet_pull->GetXaxis()->SetLabelSize(0.125);
  frm_dataJet_pull->GetXaxis()->SetTitleSize(0.125);
  frm_dataJet_pull->GetYaxis()->SetLabelSize(0.125);
  frm_dataJet_pull->GetYaxis()->SetTitleSize(0.125);
  frm_dataJet_pull->GetYaxis()->SetNdivisions(505);
  frm_dataJet_pull->SetMinimum(-4);
  frm_dataJet_pull->SetMaximum(4);
  frm_dataJet_pull->Draw();

  c7.Draw();
  c7.Print(Form("rooFit_forData_%s_cat%s.pdf", channel.data(), catcut.data()));

  TCanvas c8("c8","",0,0,1000,800);
  TLegend leg(0.60,0.70,0.85,0.80);

  c8.cd();
  frm_alpha->SetTitle("");
  frm_alpha->GetYaxis()->SetTitle("Normalization");
  frm_alpha->GetYaxis()->SetTitleOffset(1.3);
  frm_alpha->Draw();
  leg.AddEntry(frm_alpha->findObject(frm_alpha->nameOf(2)), "bkg. fit in sidebands", "l");
  leg.AddEntry(frm_alpha->findObject(frm_alpha->nameOf(3)), "bkg. fit in signal region", "l");
  leg.AddEntry(frm_alpha->findObject(frm_alpha->nameOf(1)), "#alpha function (y=exp(#frac{-x}{a+bx}))", "l");
  leg.SetBorderSize(0);
  leg.Draw();
  frm_alpha->addObject(&leg);
  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.60, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.62, 0.82, Form("%s, %s b-tag", channel.data(), catcut.data()));
  c8.Print(Form("rooFit_forData_%s_cat%s.pdf", channel.data(), catcut.data()));

  c8.Clear();
  c8.cd()->SetLogy(1);
  frm_predict->SetTitle("");
  frm_predict->SetMinimum(1e-2);
  frm_predict->SetMaximum(catcut=="1"?100:10);
  frm_predict->Draw();
  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar.DrawLatexNDC(0.60, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s btag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, "expected background in data signal region");
  c8.Print(Form("rooFit_forData_%s_cat%s.pdf)", channel.data(), catcut.data()));

}
