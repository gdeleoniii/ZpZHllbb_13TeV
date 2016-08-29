R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/PdfDiagonalizer_cc.so)
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/uncertainties/readFitParam.h"
using namespace RooFit;

void rooFitAlpha(string channel, string catcut){

  // Suppress all the INFO message

  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooMsgService::instance().setSilentMode(true);

  // Define all the variables from the trees

  RooRealVar cat ("cat", "", 0, 2);
  RooRealVar mJet("prmass", "", 30., 300.);
  RooRealVar evWeight("evweight", "", 0., 1.e3);
  RooRealVar mZH("mllbb", "M_{ZH}", 750., 4300., "GeV");
 
  mZH.setRange("fullRange", 750., 4300.);

  RooBinning binsmZH(71, 750, 4300);
 
  RooArgSet variables(cat, mJet, mZH, evWeight);
 
  TCut catCut = Form("cat==%s", catcut.c_str());
  TCut sbCut  = "prmass>30 && !(prmass>65 && prmass<135) && prmass<300";
  TCut sigCut = "prmass>105 && prmass<135";

  string region[3] = {"central","up","down"};
  float alpha[13][3];
  float Mzh[13] = {750,800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4300};

  // Histograms that stored in root file (for shape analysis)

  TH1* h_shape[3];

  for(int nw = 2; nw >= 0; --nw){

    // Input files and sum all backgrounds

    TChain* treeData  = new TChain("tree");
    TChain* treeZjets = new TChain("tree");
    
    if( channel == "ele" ){

      treeData->Add(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/data/SingleElectron-Run2015D-v1_%sMiniTree.root", channel.data()));
      treeData->Add(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/data/SingleElectron-Run2015D-v4_%sMiniTree.root", channel.data()));

    }

    else if( channel == "mu" ){

      treeData->Add(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/data/SingleMuon-Run2015D-v1_%sMiniTree.root", channel.data()));
      treeData->Add(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/data/SingleMuon-Run2015D-v4_%sMiniTree.root", channel.data()));

    }

    // To remove minor background contribution in data set (weight is -1)

    treeData->Add(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/minor/WW_TuneCUETP8M1_13TeV_%sMiniTree.root",     channel.data()));
    treeData->Add(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/minor/WZ_TuneCUETP8M1_13TeV_%sMiniTree.root",     channel.data()));
    treeData->Add(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/minor/ZZ_TuneCUETP8M1_13TeV_%sMiniTree.root",     channel.data()));
    treeData->Add(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/minor/TT_TuneCUETP8M1_13TeV_%sMiniTree.root",     channel.data()));
    treeData->Add(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/toyMCtest/minor/ZH_HToBB_ZToLL_M125_13TeV_%sMiniTree.root", channel.data()));

    // Z+jets background

    treeZjets->Add(Form("Zjets/DYJetsToLL_M-50_HT-100to200_13TeV_%s_%sMiniTree.root", region[nw].data(), channel.data()));
    treeZjets->Add(Form("Zjets/DYJetsToLL_M-50_HT-200to400_13TeV_%s_%sMiniTree.root", region[nw].data(), channel.data()));
    treeZjets->Add(Form("Zjets/DYJetsToLL_M-50_HT-400to600_13TeV_%s_%sMiniTree.root", region[nw].data(), channel.data()));
    treeZjets->Add(Form("Zjets/DYJetsToLL_M-50_HT-600toInf_13TeV_%s_%sMiniTree.root", region[nw].data(), channel.data()));
    
    // Create a dataset from a tree -> to process an unbinned likelihood fitting

    RooDataSet dataSetData   ("dataSetData",    "dataSetData",    variables, Cut(catCut),           WeightVar(evWeight), Import(*treeData));
    RooDataSet dataSetDataSB ("dataSetDataSB",  "dataSetDataSB",  variables, Cut(catCut && sbCut),  WeightVar(evWeight), Import(*treeData));
    RooDataSet dataSetDataSG ("dataSetDataSG",  "dataSetDataSG",  variables, Cut(catCut && sigCut), WeightVar(evWeight), Import(*treeData));
    RooDataSet dataSetZjetsSB("dataSetZjetsSB", "dataSetZjetsSB", variables, Cut(catCut && sbCut),  WeightVar(evWeight), Import(*treeZjets));  
    RooDataSet dataSetZjetsSG("dataSetZjetsSG", "dataSetZjetsSG", variables, Cut(catCut && sigCut), WeightVar(evWeight), Import(*treeZjets));

    // Total events number

    RooRealVar nSBMcEvents("nSBMcEvents", "nSBMcEvents", 0., 1.e10);
    RooRealVar nSGMcEvents("nSGMcEvents", "nSGMcEvents", 0., 1.e10);
    RooRealVar nDataEvents  ("nDataEvents",   "nDataEvents",   0., 1.e10);
    RooRealVar nSBDataEvents("nSBDataEvents", "nSBDataEvents", 0., 1.e10);

    nSBMcEvents.setVal(dataSetZjetsSB.sumEntries());
    nSBMcEvents.setConstant(true);
  
    nSGMcEvents.setVal(dataSetZjetsSG.sumEntries());
    nSGMcEvents.setConstant(true);

    nDataEvents.setVal(dataSetData.sumEntries());
    nDataEvents.setConstant(true);

    nSBDataEvents.setVal(dataSetDataSB.sumEntries());
    nSBDataEvents.setConstant(true);

    // Side band jet mass in data

    RooRealVar lamda("lamda", "lamda", -0.015, -0.04, -0.01);

    RooExponential model_mJetSB("model_mJetSB", "model_mJetSB", mJet, lamda);
    RooExtendPdf ext_model_mJetSB("ext_model_mJetSB", "ext_model_mJetSB", model_mJetSB, nSBDataEvents);

    RooFitResult* mJetSB_result = ext_model_mJetSB.fitTo(dataSetDataSB, SumW2Error(true), Extended(true), Range("lowSB,highSB"), Strategy(2), Minimizer("Minuit2"), Save(1));

    fprintf(stdout, "lamda value is %g\n", lamda.getVal());

    // Normalize factor to normalize the background in signal region of data

    RooAbsReal* nSIGFit = ext_model_mJetSB.createIntegral(RooArgSet(mJet), NormSet(mJet), Range("signal"));
    RooAbsReal* nSBFit  = ext_model_mJetSB.createIntegral(RooArgSet(mJet), NormSet(mJet), Range("lowSB,highSB"));
 
    RooRealVar normFactor("normFactor", "normFactor", 0., 1.e9);

    normFactor.setVal(nSBDataEvents.getVal()*(nSIGFit->getVal()/nSBFit->getVal()));
    normFactor.setConstant(true);
  
    // Alpha ratio part

    // set fit parameters

    RooRealVar sbVara("sbVara", "sbVara", param(channel.data(),catcut.data(),"sbVaraMin"), param(channel.data(),catcut.data(),"sbVaraMax"));
    RooRealVar sbVarb("sbVarb", "sbVarb", param(channel.data(),catcut.data(),"sbVarbMin"), param(channel.data(),catcut.data(),"sbVarbMax"));
    RooRealVar sgVara("sgVara", "sgVara", param(channel.data(),catcut.data(),"sgVaraMin"), param(channel.data(),catcut.data(),"sgVaraMax"));
    RooRealVar sgVarb("sgVarb", "sgVarb", param(channel.data(),catcut.data(),"sgVarbMin"), param(channel.data(),catcut.data(),"sgVarbMax"));
    RooRealVar daVara("daVara", "daVara", param(channel.data(),catcut.data(),"daVaraMin"), param(channel.data(),catcut.data(),"daVaraMax"));
    RooRealVar daVarb("daVarb", "daVarb", param(channel.data(),catcut.data(),"daVarbMin"), param(channel.data(),catcut.data(),"daVarbMax"));

    // Fit ZH mass in MC side band

    RooGenericPdf model_ZHSB("model_ZHSB", "model_ZHSB", "TMath::Exp(-@0/(@1+@2*@0))", RooArgSet(mZH,sbVara,sbVarb));
    RooExtendPdf ext_model_ZHSB("ext_model_ZHSB", "ext_model_ZHSB", model_ZHSB, nSBMcEvents);
    ext_model_ZHSB.fitTo(dataSetZjetsSB, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));

    // Fit ZH mass in MC signal region

    RooGenericPdf model_ZHSG("model_ZHSG", "model_ZHSG", "TMath::Exp(-@0/(@1+@2*@0))", RooArgSet(mZH,sgVara,sgVarb));
    RooExtendPdf ext_model_ZHSG("ext_model_ZHSG", "ext_model_ZHSG", model_ZHSG, nSGMcEvents);
    ext_model_ZHSG.fitTo(dataSetZjetsSG, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));

    // Fit ZH mass in data side band

    RooGenericPdf model_ZH("model_ZH", "model_ZH", "TMath::Exp(-@0/(@1+@2*@0))", RooArgSet(mZH,daVara,daVarb));
    RooExtendPdf  ext_model_ZH("ext_model_ZH", "ext_model_ZH", model_ZH, nDataEvents);
    ext_model_ZH.fitTo(dataSetDataSB, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));

    // Set the model of alpha ratio

    float normConst = ((TF1*)ext_model_ZHSB.asTF(mZH, RooArgList(sbVara, sbVarb)))->Integral(750,4300)/((TF1*)ext_model_ZHSG.asTF(mZH, RooArgList(sgVara, sgVarb)))->Integral(750,4300);

    TF1* f_alpha = new TF1("f_alpha", "[0]*TMath::Exp(-x/([1]+[2]*x))/TMath::Exp(-x/([3]+[4]*x))", 750, 4300);

    f_alpha->SetParameters(normConst, sgVara.getVal(), sgVarb.getVal(), sbVara.getVal(), sbVarb.getVal());

    for( int im = 0; im < 13; ++im ){
      alpha[im][nw] = f_alpha->Eval(Mzh[im]);
    }

    fprintf(stdout, "sbVara=%f\tsbVarb=%f\tsgVara=%f\tsgVarb=%f\n", sbVara.getVal(), sbVarb.getVal(), sgVara.getVal(), sgVarb.getVal());

    // predicted background model --> histogram

    RooGenericPdf model_alpha("model_alpha", "model_alpha", Form("%f*TMath::Exp(-@0/(%f+%f*@0))/TMath::Exp(-@0/(%f+%f*@0))", normConst, sgVara.getVal(), sgVarb.getVal(), sbVara.getVal(), sbVarb.getVal()), RooArgSet(mZH));
    RooProdPdf    model_sigData("model_sigData", "ext_model_ZH*model_alpha", RooArgList(ext_model_ZH,model_alpha));
    RooExtendPdf  ext_model_sigData("ext_model_sigData", "ext_model_sigData", model_sigData, normFactor);

    h_shape[nw] = ext_model_sigData.createHistogram(Form("h_shape%i",nw), mZH, Binning(binsmZH), Extended(true));

    delete treeZjets;
    delete f_alpha;

  } // end of weight loop

  // Store histograms in root file (for shape analysis)

  TFile f_shape("histo_mZH_bTagUnc.root", "recreate");

  h_shape[0]->Write("h_mZH_bTag_central");
  h_shape[1]->Write("h_mZH_bTag_up");
  h_shape[2]->Write("h_mZH_bTag_down");
   
  // Calculate uncertainty of each mass bin

  float Alpha[13], Unc[13], relativeUnc[13];

  for( int im = 0; im < 13; ++im ){

    Alpha[im] = alpha[im][0];
    Unc[im] = (fabs(alpha[im][1]-alpha[im][0])>fabs(alpha[im][2]-alpha[im][0])) ? fabs(alpha[im][1]-alpha[im][0]) : fabs(alpha[im][2]-alpha[im][0]);
    relativeUnc[im] = Unc[im]/Alpha[im];

    fprintf(stdout, "massPoint=%i\trelativeUnc=%f\n", (int)Mzh[im], relativeUnc[im]);

  } // end of mass points
  
  TGraphErrors *g_alpha = new TGraphErrors(13, Mzh, Alpha, 0, Unc);

  g_alpha->SetTitle("");
  g_alpha->GetXaxis()->SetLimits(750,4300);
  g_alpha->GetXaxis()->SetTitle("m_{ZH}(GeV)");
  g_alpha->GetYaxis()->SetTitle("#alpha Ratio");  
  g_alpha->GetYaxis()->SetTitleOffset(1.3);
  g_alpha->SetMinimum(0.05);
  g_alpha->SetMaximum(50);
  g_alpha->SetLineWidth(2);
  g_alpha->SetLineColor(kBlue);
  g_alpha->SetMarkerStyle(8);
  g_alpha->SetMarkerColor(kBlue);
  g_alpha->SetFillStyle(3002);
  
  TLatex lar;

  lar.SetTextSize(0.03);
  lar.SetLineWidth(5);

  TCanvas cv("cv","",0,0,1000,900);

  cv.cd()->SetLogy();

  g_alpha->Draw("Xac");
  g_alpha->Draw("3same");

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.60, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s b-tag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, "b-tagging scale factor");

  cv.Draw();
  cv.Print(Form("alpha_bTagScale_%s_cat%s.pdf", channel.data(), catcut.data()));

}
