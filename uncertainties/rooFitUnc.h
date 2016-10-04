R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/PdfDiagonalizer_cc.so)
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/readFitParam.h"
using namespace RooFit;

void rooFitUnc(string channel, string catcut, string region, TF1** f_alpha, TH1** h_shape, int i, string tag=""){

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

  string str1 = "";
  string str2 = "";

  if( tag=="jes" ){
    str1 = ""; 
    str2 = region+"_";
  }
  else if( tag=="pdf" ){
    str1 = ""; 
    str2 = "";
  }
  else{
    str1 = "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/"; 
    str2 = "";
  }

  if( channel == "ele" ){
    
    tree_Data->Add(Form("%sdata/SingleElectron-Run2015D-v1_%seleMiniTree.root", str1.data(), str2.data()));
    tree_Data->Add(Form("%sdata/SingleElectron-Run2015D-v4_%seleMiniTree.root", str1.data(), str2.data()));
    
  }
  
  else if( channel == "mu" ){
    
    tree_Data->Add(Form("%sdata/SingleMuon-Run2015D-v1_%smuMiniTree.root", str1.data(), str2.data()));
    tree_Data->Add(Form("%sdata/SingleMuon-Run2015D-v4_%smuMiniTree.root", str1.data(), str2.data()));
    
  }
  
  // Dominant and subdominant background

  if( tag!="pdf" ){

    tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-100to200_13TeV_%s_%sMiniTree.root", region.data(), channel.data()));
    tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-200to400_13TeV_%s_%sMiniTree.root", region.data(), channel.data()));
    tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-400to600_13TeV_%s_%sMiniTree.root", region.data(), channel.data()));
    tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-600toInf_13TeV_%s_%sMiniTree.root", region.data(), channel.data()));
    
  }

  else{

    tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-100to200_13TeV_%sMiniTree.root", channel.data()));
    tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-200to400_13TeV_%sMiniTree.root", channel.data()));
    tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-400to600_13TeV_%sMiniTree.root", channel.data()));
    tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-600toInf_13TeV_%sMiniTree.root", channel.data()));

  }

  string str3 = (tag=="pdf") ? "" : "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/";
  
  tree_Sub1->Add(Form("%sminor/WW_TuneCUETP8M1_13TeV_%sMiniTree.root",     str3.data(), channel.data()));
  tree_Sub1->Add(Form("%sminor/WZ_TuneCUETP8M1_13TeV_%sMiniTree.root",     str3.data(), channel.data()));
  tree_Sub1->Add(Form("%sminor/ZZ_TuneCUETP8M1_13TeV_%sMiniTree.root",     str3.data(), channel.data()));
  tree_Sub1->Add(Form("%sminor/TT_TuneCUETP8M1_13TeV_%sMiniTree.root",     str3.data(), channel.data()));
  tree_Sub2->Add(Form("%sminor/ZH_HToBB_ZToLL_M125_13TeV_%sMiniTree.root", str3.data(), channel.data()));
  
  // Define all the variables from the trees

  RooRealVar cat("cat", "", 0, 2);
  RooRealVar mJet("prmass", "M_{jet}", 30., 300., "GeV");
  RooRealVar mZH("mllbb", "M_{ZH}", 750., 4300., "GeV");
  RooRealVar evWeight(((tag=="pdf") ?  Form("evweight%02i",i) : "evweight"), "", 0., 1.e3);

  // Set the range in zh mass and in jet mass

  mZH .setRange("All",  750., 4300.);
  mJet.setRange("All",   30.,  300.);
  mJet.setRange("SB_l",  30.,   65.);
  mJet.setRange("SB_h", 135.,  300.);
  mJet.setRange("SG",   105.,  135.);

  RooBinning bin_mZH(71, 750, 4300);

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

  float nHist_sbDom  = set_sbDom .sumEntries();
  float nHist_sgDom  = set_sgDom .sumEntries();
  float nHist_sbSub1 = set_sbSub1.sumEntries();
  float nHist_sgSub1 = set_sgSub1.sumEntries();
  float nHist_sbSub2 = set_sbSub2.sumEntries();
  float nHist_sgSub2 = set_sgSub2.sumEntries();
  float nHist_sbData = set_sbData.sumEntries();

  RooRealVar nEv_sbDom ("nEv_sbDom",  "nEv_sbDom",  nHist_sbDom,  nHist_sbDom*0.5,  nHist_sbDom*1.5);
  RooRealVar nEv_sgDom ("nEv_sgDom",  "nEv_sgDom",  nHist_sgDom,  nHist_sgDom*0.5,  nHist_sgDom*1.5);
  RooRealVar nEv_sbSub1("nEv_sbSub1", "nEv_sbSub1", nHist_sbSub1, nHist_sbSub1*0.5, nHist_sbSub1*1.5);
  RooRealVar nEv_sgSub1("nEv_sgSub1", "nEv_sgSub1", nHist_sgSub1, nHist_sgSub1*0.5, nHist_sgSub1*1.5);
  RooRealVar nEv_sbSub2("nEv_sbSub2", "nEv_sbSub2", nHist_sbSub2, nHist_sbSub2*0.5, nHist_sbSub2*1.5);
  RooRealVar nEv_sgSub2("nEv_sgSub2", "nEv_sgSub2", nHist_sgSub2, nHist_sgSub2*0.5, nHist_sgSub2*1.5);
  RooRealVar nEv_sbData("nEv_sbData", "nEv_sbData", nHist_sbData, nHist_sbData*0.5, nHist_sbData*1.5);
  RooRealVar nEv_forJet("nEv_forJet", "nEv_forJet", nHist_sbData, nHist_sbData*0.5, nHist_sbData*1.5);

  // Set fit parameters for ZH mass

  param myVal(channel.data(), catcut.data());

  RooRealVar a_domSb ("a_domSb",  "a_domSb",  myVal.value("a_domSb"),  myVal.value("a_domSbMin"),  myVal.value("a_domSbMax"));
  RooRealVar b_domSb ("b_domSb",  "b_domSb",  myVal.value("b_domSb"),  myVal.value("b_domSbMin"),  myVal.value("b_domSbMax"));
  RooRealVar a_domSg ("a_domSg",  "a_domSg",  myVal.value("a_domSg"),  myVal.value("a_domSgMin"),  myVal.value("a_domSgMax"));
  RooRealVar b_domSg ("b_domSg",  "b_domSg",  myVal.value("b_domSg"),  myVal.value("b_domSgMin"),  myVal.value("b_domSgMax"));
  RooRealVar a_dataSb("a_dataSb", "a_dataSb", myVal.value("a_dataSb"), myVal.value("a_dataSbMin"), myVal.value("a_dataSbMax"));
  RooRealVar b_dataSb("b_dataSb", "b_dataSb", myVal.value("b_dataSb"), myVal.value("b_dataSbMin"), myVal.value("b_dataSbMax"));
  RooRealVar a_sub1Sb("a_sub1Sb", "a_sub1Sb", myVal.value("a_sub1Sb"), myVal.value("a_sub1SbMin"), myVal.value("a_sub1SbMax"));
  RooRealVar b_sub1Sb("b_sub1Sb", "b_sub1Sb", myVal.value("b_sub1Sb"), myVal.value("b_sub1SbMin"), myVal.value("b_sub1SbMax"));
  RooRealVar a_sub1Sg("a_sub1Sg", "a_sub1Sg", myVal.value("a_sub1Sg"), myVal.value("a_sub1SgMin"), myVal.value("a_sub1SgMax"));
  RooRealVar b_sub1Sg("b_sub1Sg", "b_sub1Sg", myVal.value("b_sub1Sg"), myVal.value("b_sub1SgMin"), myVal.value("b_sub1SgMax"));
  RooRealVar a_sub2Sb("a_sub2Sb", "a_sub2Sb", myVal.value("a_sub2Sb"), myVal.value("a_sub2SbMin"), myVal.value("a_sub2SbMax"));
  RooRealVar b_sub2Sb("b_sub2Sb", "b_sub2Sb", myVal.value("b_sub2Sb"), myVal.value("b_sub2SbMin"), myVal.value("b_sub2SbMax"));
  RooRealVar a_sub2Sg("a_sub2Sg", "a_sub2Sg", myVal.value("a_sub2Sg"), myVal.value("a_sub2SgMin"), myVal.value("a_sub2SgMax"));
  RooRealVar b_sub2Sg("b_sub2Sg", "b_sub2Sg", myVal.value("b_sub2Sg"), myVal.value("b_sub2SgMin"), myVal.value("b_sub2SgMax"));

  // Create pdf for ZH mass

  RooGenericPdf pdf_sbDomZh ("pdf_sbDomZh",  "pdf_sbDomZh",  "exp(-@0/(@1+@2*@0))", RooArgSet(mZH, a_domSb , b_domSb));
  RooGenericPdf pdf_sgDomZh ("pdf_sgDomZh",  "pdf_sgDomZh",  "exp(-@0/(@1+@2*@0))", RooArgSet(mZH, a_domSg , b_domSg));
  RooGenericPdf pdf_sbSub1Zh("pdf_sbSub1Zh", "pdf_sbSub1Zh", "exp(-@0/(@1+@2*@0))", RooArgSet(mZH, a_sub1Sb, b_sub1Sb));
  RooGenericPdf pdf_sgSub1Zh("pdf_sgSub1Zh", "pdf_sgSub1Zh", "exp(-@0/(@1+@2*@0))", RooArgSet(mZH, a_sub1Sg, b_sub1Sg));
  RooGenericPdf pdf_sbSub2Zh("pdf_sbSub2Zh", "pdf_sbSub2Zh", "exp(-@0/(@1+@2*@0))", RooArgSet(mZH, a_sub2Sb, b_sub2Sb));
  RooGenericPdf pdf_sgSub2Zh("pdf_sgSub2Zh", "pdf_sgSub2Zh", "exp(-@0/(@1+@2*@0))", RooArgSet(mZH, a_sub2Sg, b_sub2Sg));
  RooGenericPdf pdf_sbDataZh("pdf_sbDataZh", "pdf_sbDataZh", "exp(-@0/(@1+@2*@0))", RooArgSet(mZH, a_dataSb, b_dataSb));

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

  RooDataSet set_combine("set_combine", "set_combine", RooArgSet(cat, mJet, mZH, evWeight), Index(cat_combine), WeightVar(evWeight), Import(mapSet));

  RooSimultaneous pdf_combine("pdf_combine", "pdf_combine", cat_combine);

  pdf_combine.addPdf(ext_sbDomZh,  "dom_SB");
  pdf_combine.addPdf(ext_sgDomZh,  "dom_SG");
  pdf_combine.addPdf(ext_sbSub1Zh, "sub1_SB");
  pdf_combine.addPdf(ext_sgSub1Zh, "sub1_SG");
  pdf_combine.addPdf(ext_sbSub2Zh, "sub2_SB");
  pdf_combine.addPdf(ext_sgSub2Zh, "sub2_SG");
  pdf_combine.addPdf(ext_sbDataZh, "data_SB");

  RooLinkedList cmdList;

  cmdList.Add(Save(true).Clone());
  cmdList.Add(Minos(true).Clone());  
  cmdList.Add(Offset(true).Clone());
  cmdList.Add(Extended(true).Clone());
  cmdList.Add(SumW2Error(false).Clone());
  cmdList.Add(NumCPU(8).Clone());
  cmdList.Add(Strategy(2).Clone());
  cmdList.Add(Range("All").Clone());
  cmdList.Add(Minimizer("Minuit2","migrad").Clone());

  pdf_combine.fitTo(set_combine, cmdList);

  // Multiply the model of background in data side band with the model of alpha ratio to the a model of background in data signal region
  // predicted background = (sbDataZh - sbSub1Zh - sbSub2Zh) * alpha + sgSub1Zh + sgSub2Zh

  double param_alpha[4] = {a_domSg.getVal(), 
			   b_domSg.getVal(), 
			   a_domSb.getVal(), 
			   b_domSb.getVal()};
  
  double param_predict[19] = {nEv_sbData.getVal(), 
			      a_dataSb  .getVal(),
			      b_dataSb  .getVal(),
			      nEv_sbSub1.getVal(), 
			      a_sub1Sb  .getVal(),
			      b_sub1Sb  .getVal(), 
			      nEv_sbSub2.getVal(),
			      a_sub2Sb  .getVal(),
			      b_sub2Sb  .getVal(), 
			      a_domSg   .getVal(), 
			      b_domSg   .getVal(),
			      a_domSb   .getVal(), 
			      b_domSb   .getVal(),
			      nEv_sgSub1.getVal(),
			      a_sub1Sg  .getVal(),
			      b_sub1Sg  .getVal(), 
			      nEv_sgSub2.getVal(), 
			      a_sub2Sg  .getVal(),
			      b_sub1Sg  .getVal()};

  // Normalization correction

  float corr_sbDom  = 1/ext_sbDomZh .createIntegral(mZH)->getVal();
  float corr_sgDom  = 1/ext_sgDomZh .createIntegral(mZH)->getVal();
  float corr_sbSub1 = 1/ext_sbSub1Zh.createIntegral(mZH)->getVal();
  float corr_sgSub1 = 1/ext_sgSub1Zh.createIntegral(mZH)->getVal();
  float corr_sbSub2 = 1/ext_sbSub2Zh.createIntegral(mZH)->getVal();
  float corr_sgSub2 = 1/ext_sgSub2Zh.createIntegral(mZH)->getVal();
  float corr_sbData = 1/ext_sbDataZh.createIntegral(mZH)->getVal();

  string eqn_alpha   = Form("(%f*exp(-x/([0]+[1]*x)))/(%f*exp(-x/([2]+[3]*x)))", corr_sgDom, corr_sbDom);
  string eqn_predict = Form("%f*(%f*[0]*exp(-x/([1]+[2]*x))-%f*[3]*exp(-x/([4]+[5]*x))-%f*[6]*exp(-x/([7]+[8]*x))*(%f*exp(-x/([9]+[10]*x)))/(%f*exp(-x/([11]+[12]*x)))+%f*[13]*exp(-x/([14]+[15]*x))+%f*[16]*exp(-x/([17]+[18]*x)))", bin_mZH.binWidth(1), corr_sbData, corr_sbSub1, corr_sbSub2, corr_sgDom, corr_sbDom, corr_sgSub1, corr_sgSub2);


  *f_alpha = new TF1(Form("f_alpha%i",i), eqn_alpha.data(), bin_mZH.lowBound(), bin_mZH.highBound());
  (*f_alpha)->SetParameters(param_alpha);

  TF1 f_predict("f_predict", eqn_predict.data(), bin_mZH.lowBound(), bin_mZH.highBound());
  f_predict.SetParameters(param_predict);

  // jet mass in data side band

  RooRealVar    j_data("j_data", "j_data", myVal.value("j_data"), myVal.value("j_dataMin"), myVal.value("j_dataMax"));
  RooGenericPdf pdf_dataJet("pdf_dataJet", "pdf_dataJet", "exp(-@0/@1)", RooArgSet(mJet, j_data));
  RooExtendPdf  ext_dataJet("ext_dataJet", "ext_dataJet", pdf_dataJet, nEv_forJet);

  RooLinkedList cmdjetList;

  cmdjetList.Add(Save(true).Clone());
  cmdjetList.Add(Minos(true).Clone());
  cmdjetList.Add(Offset(true).Clone());
  cmdjetList.Add(Extended(true).Clone());
  cmdjetList.Add(SumW2Error(false).Clone());
  cmdjetList.Add(NumCPU(8).Clone());
  cmdjetList.Add(Strategy(2).Clone());
  cmdjetList.Add(Range("SB_l,SB_h").Clone());
  cmdjetList.Add(Minimizer("Minuit2","migrad").Clone());

  ext_dataJet.fitTo(set_sbData, cmdjetList);

  // Normalize factor to normalize the background in signal region of data

  RooAbsReal* nFit_sg = ext_dataJet.createIntegral(mJet, Range("SG"));
  RooAbsReal* nFit_sb = ext_dataJet.createIntegral(mJet, Range("SB_l,SB_h"));

  // Since the statistic of 2015 data is low, the jet mass distribution in 2 btag is consider as a flat distribution

  float normFactorVal = (catcut=="1") ? nEv_sbData.getVal()*(nFit_sg->getVal()/nFit_sb->getVal()) : 6;

  // Convert TF1 to TH1
  
  *h_shape = new TH1D(Form("h_shape%i",i), "", bin_mZH.numBins(), bin_mZH.lowBound(), bin_mZH.highBound());

  float a = bin_mZH.lowBound();
  float b = a + bin_mZH.binWidth(1);
  
  for( int n = 1; n <= bin_mZH.numBins(); ++n ){

    (*h_shape)->SetBinContent(n, f_predict.Integral(a, b)/bin_mZH.binWidth(1));

    a += bin_mZH.binWidth(1);
    b += bin_mZH.binWidth(1);
    
  }

  (*h_shape)->Scale(normFactorVal/(*h_shape)->Integral());

  fprintf(stdout, "****nw=%i****\n", i);
  fprintf(stdout, "a_domSb=%.3f+-%.3f\n", a_domSb.getVal(), a_domSb.getError());
  fprintf(stdout, "b_domSb=%.3f+-%.3f\n", b_domSb.getVal(), b_domSb.getError());
  fprintf(stdout, "a_domSg=%.3f+-%.3f\n", a_domSg.getVal(), a_domSg.getError());
  fprintf(stdout, "b_domSg=%.3f+-%.3f\n", b_domSg.getVal(), b_domSg.getError());

  if( tag == "jes" ){
    fprintf(stdout, "nEv_sbData=%.3f+-%.3f\n", nEv_sbData.getVal(), nEv_sbData.getError());
    fprintf(stdout, "a_dataSb=%.3f+-%.3f\n", a_dataSb.getVal(), a_dataSb.getError());
    fprintf(stdout, "b_dataSb=%.3f+-%.3f\n", b_dataSb.getVal(), b_dataSb.getError());
    fprintf(stdout, "j_data=%.3f+-%.3f\n", j_data.getVal(), j_data.getError());
    fprintf(stdout, "normFactor=%.3f\n", normFactorVal);
  }

}
