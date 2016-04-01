R__LOAD_LIBRARY(PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(PDFs/PdfDiagonalizer_cc.so)
using namespace RooFit;
using namespace std;

void forData(string channel, string catcut, bool removeMinor=true){

  // Suppress all the INFO message

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  // Input files and sum all backgrounds

  TChain* treeData  = new TChain("tree");
  TChain* treeZjets = new TChain("tree");

  if( channel == "ele" ){

    treeData->Add(Form("%s/data/SingleElectron-Run2015D-05Oct2015-v1_toyMCnew.root",  channel.data()));
    treeData->Add(Form("%s/data/SingleElectron-Run2015D-PromptReco-V4_toyMCnew.root", channel.data()));

  }

  else if( channel == "mu" ){

    treeData->Add(Form("%s/data/SingleMuon-Run2015D-05Oct2015-v1_toyMCnew.root",  channel.data()));
    treeData->Add(Form("%s/data/SingleMuon-Run2015D-PromptReco-V4_toyMCnew.root", channel.data()));

  }

  else return;

  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-100to200_13TeV_toyMCnew.root", channel.data()));
  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-200to400_13TeV_toyMCnew.root", channel.data()));
  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-400to600_13TeV_toyMCnew.root", channel.data()));
  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-600toInf_13TeV_toyMCnew.root", channel.data()));

  // To remove minor background contribution in data set (weight is -1)

  if( removeMinor ){

    treeData->Add(Form("%s/VV/WW_TuneCUETP8M1_13TeV_toyMCnew.root", channel.data()));
    treeData->Add(Form("%s/VV/WZ_TuneCUETP8M1_13TeV_toyMCnew.root", channel.data()));
    treeData->Add(Form("%s/VV/ZZ_TuneCUETP8M1_13TeV_toyMCnew.root", channel.data()));
    treeData->Add(Form("%s/TT/TT_TuneCUETP8M1_13TeV_toyMCnew.root", channel.data()));

  }

  // Define all the variables from the trees

  RooRealVar cat ("cat", "", 0, 2);
  RooRealVar mJet("prmass", "M_{jet}",  30.,  300., "GeV");
  RooRealVar mZH ("mllbb",   "M_{ZH}", 750., 4000., "GeV");
  RooRealVar evWeight("evweight", "", 0, 1.e3);

  // Set the range in jet mass

  mJet.setRange("allRange", 30., 300.);
  mJet.setRange("lowSB",    30.,  65.);
  mJet.setRange("highSB",  135., 300.);
  mJet.setRange("signal",  105., 135.);

  RooBinning binsmJet(54, 30, 300);

  RooArgSet variables(cat, mJet, mZH, evWeight);

  TCut catCut = Form("cat==%s", catcut.c_str());
  TCut sbCut  = "prmass>30 && !(prmass>65 && prmass<135) && prmass<300";
  TCut sigCut = "prmass>105 && prmass<135";

  // Create a dataset from a tree -> to process an unbinned likelihood fitting

  RooDataSet dataSetData("dataSetData",
			 "dataSetData",
			 variables,
			 Cut(catCut),
			 WeightVar(evWeight),
			 Import(*treeData));

  for(int i=0; i< dataSetData.sumEntries();i++){
    dataSetData.get(i);
    cout << dataSetData.weight()<<endl;;

  }

  RooDataSet dataSetDataSB("dataSetDataSB", 
			   "dataSetDataSB",
			   variables,
			   Cut(catCut && sbCut),
			   WeightVar(evWeight),
			   Import(*treeData));

  RooDataSet dataSetZjetsSB("dataSetZjetsSB",
			    "dataSetZjetsSB",
			    variables,
			    Cut(catCut && sbCut),
			    WeightVar(evWeight),
			    Import(*treeZjets));
  
  RooDataSet dataSetZjetsSG("dataSetZjetsSG",
			    "dataSetZjetsSG",
			    variables,
			    Cut(catCut && sigCut),
			    WeightVar(evWeight),
			    Import(*treeZjets));
  
  // Total events number

  float totalMcEv   = dataSetZjetsSB.sumEntries() + dataSetZjetsSG.sumEntries();
  float totalDataEv = dataSetData.sumEntries();

  RooRealVar nMcEvents("nMcEvents", "nMcEvents", 0., 99999.);
  RooRealVar nDataEvents("nDataEvents", "nDataEvents", 0., 99999.);

  nMcEvents.setVal(totalMcEv);
  nMcEvents.setConstant(true);

  nDataEvents.setVal(totalDataEv);
  nDataEvents.setConstant(true);

  // Define the variables for pdf

  RooRealVar constant("constant", "constant", -0.02,  -1.,   0.);
  RooRealVar offset  ("offset",   "offset",     30., -50., 200.);
  RooRealVar width   ("width",    "width",     100.,   0., 200.);

  // Set the parameter fixed when fitting

  offset.setConstant(true);

  // Define the pdf and fit data
  
  RooErfExpPdf model_mJet("model_mJet", "model_mJet", mJet, constant, offset, width);

  RooExtendPdf ext_model_mJet("ext_model_mJet", "ext_model_mJet", model_mJet, nMcEvents);

  RooFitResult* mJet_result = ext_model_mJet.fitTo(dataSetZjetsSB,
						   SumW2Error(true),
						   Extended(true),
						   Range("lowSB,highSB"),
						   Strategy(2),
						   Minimizer("Minuit2"),
						   Save(1));

  RooAbsReal* nSIGFit = ext_model_mJet.createIntegral(RooArgSet(mJet), NormSet(mJet), Range("signal"));

  float normFactor = nSIGFit->getVal() * totalMcEv;
  
  // Plot the results on a frame

  RooPlot* mJetFrame = mJet.frame();
  
  dataSetZjetsSB.plotOn(mJetFrame, Binning(binsmJet));  
  ext_model_mJet.plotOn(mJetFrame, Range("allRange"), VisualizeError(*mJet_result), FillColor(kYellow));
  dataSetZjetsSB.plotOn(mJetFrame, Binning(binsmJet));  
  ext_model_mJet.plotOn(mJetFrame, Range("allRange"));


  /*******************************************/
  /*               Alpha ratio               */
  /*******************************************/

  mZH.setRange("fullRange", 750., 4000.);

  RooBinning binsmZH(65, 750, 4000);

  RooRealVar a("a", "a",  0., -1.,    1.);
  RooRealVar b("b", "b", 1000,  0., 4000.);
  
  RooGenericPdf model_ZHSB("model_ZHSB", "model_ZHSB", "TMath::Exp(@1*@0+@2/@0)", RooArgSet(mZH,a,b));
  RooGenericPdf model_ZHSG("model_ZHSG", "model_ZHSG", "TMath::Exp(@1*@0+@2/@0)", RooArgSet(mZH,a,b));
  RooGenericPdf model_ZH  ("model_ZH",   "model_ZH",   "TMath::Exp(@1*@0+@2/@0)", RooArgSet(mZH,a,b));

  RooExtendPdf ext_model_ZHSB("ext_model_ZHSB", "ext_model_ZHSB", model_ZHSB, nMcEvents);
  RooExtendPdf ext_model_ZHSG("ext_model_ZHSG", "ext_model_ZHSG", model_ZHSG, nMcEvents);
  RooExtendPdf ext_model_ZH  ("ext_model_ZH",   "ext_model_ZH",   model_ZH,   nDataEvents);

  // Fit ZH mass in side band  

  RooFitResult* mZHSB_result = ext_model_ZHSB.fitTo(dataSetZjetsSB, 
						    SumW2Error(true), 
						    Extended(true),
						    Range("fullRange"),
						    Strategy(2),
						    Minimizer("Minuit2"), 
						    Save(1));

  float p0 = a.getVal();
  float p1 = b.getVal();

  // Fit ZH mass in signal region

  RooFitResult* mZHSG_result = ext_model_ZHSG.fitTo(dataSetZjetsSG,
						    SumW2Error(true),
						    Extended(true),
						    Range("fullRange"),
						    Strategy(2),
						    Minimizer("Minuit2"),
						    Save(1));

  float p2 = a.getVal();
  float p3 = b.getVal();

  // Fit ZH mass in side band region (data)

  RooFitResult* mZH_result = ext_model_ZH.fitTo(dataSetDataSB,
						SumW2Error(true),
						Extended(true),
						Range("fullRange"),
						Strategy(2),
						Minimizer("Minuit2"),
						Save(1));

  // Plot the results to a frame

  RooPlot* mZHFrameMC = mZH.frame();
  RooPlot* mZHFrame   = mZH.frame();

  dataSetZjetsSB.plotOn(mZHFrameMC, Binning(binsmZH));
  ext_model_ZHSB.plotOn(mZHFrameMC, VisualizeError(*mZHSB_result), FillColor(kYellow));
  dataSetZjetsSB.plotOn(mZHFrameMC, Binning(binsmZH));
  ext_model_ZHSB.plotOn(mZHFrameMC, LineStyle(7), LineColor(kBlue));
  
  dataSetZjetsSG.plotOn(mZHFrameMC, Binning(binsmZH));
  ext_model_ZHSG.plotOn(mZHFrameMC, VisualizeError(*mZHSG_result), FillColor(kYellow));
  dataSetZjetsSG.plotOn(mZHFrameMC, Binning(binsmZH));
  ext_model_ZHSG.plotOn(mZHFrameMC, LineStyle(7), LineColor(kRed));

  dataSetDataSB.plotOn(mZHFrame, Binning(binsmZH));
  ext_model_ZH .plotOn(mZHFrame, VisualizeError(*mZH_result), FillColor(kYellow));
  dataSetDataSB.plotOn(mZHFrame, Binning(binsmZH));
  ext_model_ZH .plotOn(mZHFrame, LineStyle(7), LineColor(kBlue));
  
  // Draw the model of alpha ratio

  RooGenericPdf model_alpha("model_alpha", "model_alpha", Form("TMath::Exp(%f*@0+%f/@0)/TMath::Exp(%f*@0+%f/@0)", p2,p3,p0,p1), RooArgSet(mZH));
  
  RooPlot* alphaFrame = mZH.frame();
  model_alpha.plotOn(alphaFrame);

  // Multiply the model of background in data side band with the model of alpha ratio to the a model of background in data signal region 

  RooProdPdf model_sigData("model_sigData", "ext_model_ZH*model_alpha", RooArgList(ext_model_ZH,model_alpha));
  model_sigData.plotOn(mZHFrame, Normalization(normFactor, RooAbsReal::NumEvent), LineStyle(7), LineColor(kRed));

  /*******************************************/
  /*                 OUTPUT                  */
  /*******************************************/

  TCanvas* c = new TCanvas("c","",0,0,1000,800);

  c->cd();
  mZHFrameMC->Draw();
  c->Print(Form("rooFit_forData_%s_cat%s.pdf(", channel.data(), catcut.data()));

  c->cd();
  mZHFrame->Draw();
  c->Print(Form("rooFit_forData_%s_cat%s.pdf", channel.data(), catcut.data()));

  c->cd();
  mJetFrame->Draw();
  c->Print(Form("rooFit_forData_%s_cat%s.pdf)", channel.data(), catcut.data()));

}
