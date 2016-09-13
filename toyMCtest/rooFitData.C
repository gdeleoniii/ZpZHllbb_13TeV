R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/PdfDiagonalizer_cc.so)
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/readFitParam.h"
using namespace RooFit;

void rooFitData(string channel, string catcut){

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

    tree_Data->Add("data/SingleElectron-Run2015D-v1_eleMiniTree.root");
    tree_Data->Add("data/SingleElectron-Run2015D-v4_eleMiniTree.root");

  }

  else if( channel == "mu" ){

    tree_Data->Add("data/SingleMuon-Run2015D-v1_muMiniTree.root");
    tree_Data->Add("data/SingleMuon-Run2015D-v4_muMiniTree.root");

  }

  // Dominant background

  tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-100to200_13TeV_%sMiniTree.root", channel.data()));
  tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-200to400_13TeV_%sMiniTree.root", channel.data()));
  tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-400to600_13TeV_%sMiniTree.root", channel.data()));
  tree_Dom->Add(Form("Zjets/DYJetsToLL_M-50_HT-600toInf_13TeV_%sMiniTree.root", channel.data()));

  // Subdominant background

  tree_Sub->Add(Form("minor/WW_TuneCUETP8M1_13TeV_%sMiniTree.root",     channel.data()));
  tree_Sub->Add(Form("minor/WZ_TuneCUETP8M1_13TeV_%sMiniTree.root",     channel.data()));
  tree_Sub->Add(Form("minor/ZZ_TuneCUETP8M1_13TeV_%sMiniTree.root",     channel.data()));
  tree_Sub->Add(Form("minor/TT_TuneCUETP8M1_13TeV_%sMiniTree.root",     channel.data()));
  tree_Sub->Add(Form("minor/ZH_HToBB_ZToLL_M125_13TeV_%sMiniTree.root", channel.data()));

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

  // Plot the results on frame 

  RooPlot* frm_sbDomZh       = mZH.frame();
  RooPlot* frm_sgDomZh       = mZH.frame();
  RooPlot* frm_alpha         = mZH.frame();
  RooPlot* frm_sbDataZh      = mZH.frame();
  RooPlot* frm_sbDataJet     = mJet.frame();
  RooPlot* frm_expected      = mZH.frame();
  RooPlot* frm_sbSubZh       = mZH.frame();
  RooPlot* frm_sgSubZh       = mZH.frame(); 
  RooPlot* frm_sbDomZh_pull  = mZH.frame();
  RooPlot* frm_sgDomZh_pull  = mZH.frame();
  RooPlot* frm_sbSubZh_pull  = mZH.frame();
  RooPlot* frm_sgSubZh_pull  = mZH.frame();
  RooPlot* frm_sbDataZh_pull = mZH.frame();
  RooPlot* frm_sbDataJet_pull= mJet.frame();

  cmb_domData.plotOn(frm_sbDomZh, Cut("cat_domData==cat_domData::dom_SB"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_domData.plotOn(frm_sbDomZh, Slice(cat_domData,"dom_SB"), ProjWData(cat_domData,cmb_domData), VisualizeError(*res_domData,1,false), FillStyle(3002));
  cmb_domData.plotOn(frm_sbDomZh, Cut("cat_domData==cat_domData::dom_SB"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_domData.plotOn(frm_sbDomZh, Slice(cat_domData,"dom_SB"), ProjWData(cat_domData,cmb_domData), LineColor(kBlue));

  cmb_domData.plotOn(frm_sgDomZh, Cut("cat_domData==cat_domData::dom_SG"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_domData.plotOn(frm_sgDomZh, Slice(cat_domData,"dom_SG"), ProjWData(cat_domData,cmb_domData), VisualizeError(*res_domData,1,false), FillStyle(3002));
  cmb_domData.plotOn(frm_sgDomZh, Cut("cat_domData==cat_domData::dom_SG"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_domData.plotOn(frm_sgDomZh, Slice(cat_domData,"dom_SG"), ProjWData(cat_domData,cmb_domData), LineColor(kBlue));

  alpha_display.plotOn(frm_alpha, VisualizeError(*res_domData,1,false), FillStyle(3002), FillColor(kBlack));
  alpha_display.plotOn(frm_alpha, LineColor(kBlack));
  ext_sbDomZh  .plotOn(frm_alpha, Normalization(1, RooAbsReal::NumEvent), LineColor(kBlue));
  ext_sgDomZh  .plotOn(frm_alpha, Normalization(1, RooAbsReal::NumEvent), LineColor(kRed));
  
  cmb_domData.plotOn(frm_sbDataZh, Cut("cat_domData==cat_domData::data_SB"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_domData.plotOn(frm_sbDataZh, Slice(cat_domData,"data_SB"), ProjWData(cat_domData,cmb_domData), VisualizeError(*res_domData,1,false), FillStyle(3002));
  cmb_domData.plotOn(frm_sbDataZh, Cut("cat_domData==cat_domData::data_SB"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_domData.plotOn(frm_sbDataZh, Slice(cat_domData,"data_SB"), ProjWData(cat_domData,cmb_domData), LineColor(kBlue));
  
  set_sbData   .plotOn(frm_sbDataJet, DataError(RooAbsData::SumW2), Binning(bin_mJet));
  ext_sbDataJet.plotOn(frm_sbDataJet, Range("all"), VisualizeError(*res_sbDataJet,1,false), FillStyle(3002));
  set_sbData   .plotOn(frm_sbDataJet, DataError(RooAbsData::SumW2), Binning(bin_mJet));
  ext_sbDataJet.plotOn(frm_sbDataJet, Range("all"));
  
  set_sgData   .plotOn(frm_expected, DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_predicted.plotOn(frm_expected, VisualizeError(*res_domData,1,false), Normalization(normFactor.getVal(), RooAbsReal::NumEvent), FillStyle(3002));
  set_sgData   .plotOn(frm_expected, DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_predicted.plotOn(frm_expected, Normalization(normFactor.getVal(), RooAbsReal::NumEvent), LineColor(kBlue));
  // Using RooAbsReal::NumEvent in order to consider the bin width of data set. Equivalent to (normFactor*binWidth) if using RooAbsReal::Raw.
  
  cmb_sub.plotOn(frm_sbSubZh, Cut("cat_sub==cat_sub::sub_SB"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_sub.plotOn(frm_sbSubZh, Slice(cat_sub,"sub_SB"), ProjWData(cat_sub,cmb_sub), VisualizeError(*res_sub,1,false), FillStyle(3002));
  cmb_sub.plotOn(frm_sbSubZh, Cut("cat_sub==cat_sub::sub_SB"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_sub.plotOn(frm_sbSubZh, Slice(cat_sub,"sub_SB"), ProjWData(cat_sub,cmb_sub), LineColor(kBlue));
 
  cmb_sub.plotOn(frm_sgSubZh, Cut("cat_sub==cat_sub::sub_SG"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_sub.plotOn(frm_sgSubZh, Slice(cat_sub,"sub_SG"), ProjWData(cat_sub,cmb_sub), VisualizeError(*res_sub,1,false), FillStyle(3002));
  cmb_sub.plotOn(frm_sgSubZh, Cut("cat_sub==cat_sub::sub_SG"), DataError(RooAbsData::SumW2), Binning(bin_mZH));
  pdf_sub.plotOn(frm_sgSubZh, Slice(cat_sub,"sub_SG"), ProjWData(cat_sub,cmb_sub), LineColor(kBlue));

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

  frm_sbDataJet->SetTitle("");
  frm_sbDataJet->SetMinimum(1e-1);
  frm_sbDataJet->SetMaximum(100);
  frm_sbDataJet->GetXaxis()->SetTitle("");
  frm_sbDataJet->GetXaxis()->SetLabelOffset(999);
  frm_sbDataJet->Draw();

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar.DrawLatexNDC(0.60, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s btag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, "data jet mass in sidebands");
  lar.DrawLatexNDC(0.15, 0.78, Form("Normalization factor: %.3f#pm%.3f", normFactor.getVal(), normFormula.getPropagatedError(*res_sbDataJet)));
  
  c3_up->RedrawAxis();
  c3_dw->cd()->SetLogy(0);

  frm_sbDataJet_pull->addObject(frm_sbDataJet->pullHist(), "P");
  frm_sbDataJet_pull->SetTitle("");
  frm_sbDataJet_pull->GetYaxis()->SetTitle("Pulls");
  frm_sbDataJet_pull->GetYaxis()->SetTitleOffset(0.25);
  frm_sbDataJet_pull->GetXaxis()->SetLabelSize(0.125);
  frm_sbDataJet_pull->GetXaxis()->SetTitleSize(0.125);
  frm_sbDataJet_pull->GetYaxis()->SetLabelSize(0.125);
  frm_sbDataJet_pull->GetYaxis()->SetTitleSize(0.125);
  frm_sbDataJet_pull->GetYaxis()->SetNdivisions(505);
  frm_sbDataJet_pull->SetMinimum(-4);
  frm_sbDataJet_pull->SetMaximum(4);
  frm_sbDataJet_pull->Draw();

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

  frm_sbSubZh->SetTitle("");
  frm_sbSubZh->SetMinimum(1e-4);
  frm_sbSubZh->SetMaximum(10);
  frm_sbSubZh->GetXaxis()->SetTitle("");
  frm_sbSubZh->GetXaxis()->SetLabelOffset(999);
  frm_sbSubZh->Draw();

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.60, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s btag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, "Subdominant background in sidebands");
  
  c4_up->RedrawAxis();
  c4_dw->cd()->SetLogy(0);

  frm_sbSubZh_pull->addObject(frm_sbSubZh->pullHist(), "P");
  frm_sbSubZh_pull->SetTitle("");
  frm_sbSubZh_pull->GetYaxis()->SetTitle("Pulls");
  frm_sbSubZh_pull->GetYaxis()->SetTitleOffset(0.25);
  frm_sbSubZh_pull->GetXaxis()->SetLabelSize(0.125);
  frm_sbSubZh_pull->GetXaxis()->SetTitleSize(0.125);
  frm_sbSubZh_pull->GetYaxis()->SetLabelSize(0.125);
  frm_sbSubZh_pull->GetYaxis()->SetTitleSize(0.125);
  frm_sbSubZh_pull->GetYaxis()->SetNdivisions(505);
  frm_sbSubZh_pull->SetMinimum(-4);
  frm_sbSubZh_pull->SetMaximum(4);
  frm_sbSubZh_pull->Draw();

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

  frm_sgSubZh->SetTitle("");
  frm_sgSubZh->SetMinimum(1e-4);
  frm_sgSubZh->SetMaximum(10);
  frm_sgSubZh->GetXaxis()->SetTitle("");
  frm_sgSubZh->GetXaxis()->SetLabelOffset(999);
  frm_sgSubZh->Draw();

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.65, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s b-tag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, "Subdominant in signal region");

  c5_up->RedrawAxis();
  c5_dw->cd()->SetLogy(0);

  frm_sgSubZh_pull->addObject(frm_sgSubZh->pullHist(), "P");
  frm_sgSubZh_pull->SetTitle("");
  frm_sgSubZh_pull->GetYaxis()->SetTitle("Pulls");
  frm_sgSubZh_pull->GetYaxis()->SetTitleOffset(0.25);
  frm_sgSubZh_pull->GetXaxis()->SetLabelSize(0.125);
  frm_sgSubZh_pull->GetXaxis()->SetTitleSize(0.125);
  frm_sgSubZh_pull->GetYaxis()->SetLabelSize(0.125);
  frm_sgSubZh_pull->GetYaxis()->SetTitleSize(0.125);
  frm_sgSubZh_pull->GetYaxis()->SetNdivisions(505);
  frm_sgSubZh_pull->SetMinimum(-4);
  frm_sgSubZh_pull->SetMaximum(4);
  frm_sgSubZh_pull->Draw();

  c5.Draw();
  c5.Print(Form("rooFit_forData_%s_cat%s.pdf", channel.data(), catcut.data()));  

  TCanvas cv("cv","",0,0,1000,800);
  TLegend leg(0.60,0.70,0.85,0.80);

  cv.cd();
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
  cv.Print(Form("rooFit_forData_%s_cat%s.pdf", channel.data(), catcut.data()));

  cv.Clear();
  cv.cd()->SetLogy(1);
  frm_expected->SetTitle("");
  frm_expected->SetMinimum(1e-1);
  frm_expected->SetMaximum(catcut=="1"?100:10);
  frm_expected->Draw();
  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar.DrawLatexNDC(0.60, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s btag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, "expected background in data signal region");
  cv.Print(Form("rooFit_forData_%s_cat%s.pdf)", channel.data(), catcut.data()));

}
