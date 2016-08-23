R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/PdfDiagonalizer_cc.so)
using namespace RooFit;

void rooFitNorm(string channel, string catcut, string study, bool anotherModel=false, bool removeMinor=true){

  // Suppress all the INFO message

  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooMsgService::instance().setSilentMode(true);

  // Define all the variables from the trees

  RooRealVar cat ("cat", "", 0, 2);
  RooRealVar mJet("prmass", "M_{jet}", 30., 300., "GeV");
  RooRealVar evWeight("evweight", "", 0., 1.e3);

  // Set the range in zh mass and in jet mass

  mJet.setRange("allRange", 30., 300.);
  mJet.setRange("lowSB",    30.,  65.);
  mJet.setRange("highSB",  135., 300.);
  mJet.setRange("signal",  105., 135.);

  RooBinning binsmJet(54, 30, 300);

  RooArgSet variables(cat, mJet, evWeight);

  TCut catCut = Form("cat==%s", catcut.c_str());
  TCut sbCut  = "prmass>30 && !(prmass>65 && prmass<135) && prmass<300";

  float normFactor[3] = {0};
  string region[3] = {"central","up","down"};

  RooPlot* mJetFrame = mJet.frame();

  for(int nw = 2; nw >= 0; --nw){

    if( (anotherModel || !removeMinor) && nw != 0 ) continue; // if study no remove minor backgroung or using another model

    // Input files and sum all backgrounds

    TChain* treeData  = new TChain("tree");

    if( channel == "ele" ){

      treeData->Add(Form("data/SingleElectron-Run2015D-v1_%s_%sMiniTree.root", region[nw].data(), channel.data()));
      treeData->Add(Form("data/SingleElectron-Run2015D-v4_%s_%sMiniTree.root", region[nw].data(), channel.data()));

    }

    else if( channel == "mu" ){

      treeData->Add(Form("data/SingleMuon-Run2015D-v1_%s_%sMiniTree.root", region[nw].data(), channel.data()));
      treeData->Add(Form("data/SingleMuon-Run2015D-v4_%s_%sMiniTree.root", region[nw].data(), channel.data()));

    }

    else return;

    // To remove minor background contribution in data set (weight is -1)

    if( removeMinor ){

      treeData->Add(Form("minor/WW_TuneCUETP8M1_13TeV_%s_%sMiniTree.root",     region[nw].data(), channel.data()));
      treeData->Add(Form("minor/WZ_TuneCUETP8M1_13TeV_%s_%sMiniTree.root",     region[nw].data(), channel.data()));
      treeData->Add(Form("minor/ZZ_TuneCUETP8M1_13TeV_%s_%sMiniTree.root",     region[nw].data(), channel.data()));
      treeData->Add(Form("minor/TT_TuneCUETP8M1_13TeV_%s_%sMiniTree.root",     region[nw].data(), channel.data()));
      treeData->Add(Form("minor/ZH_HToBB_ZToLL_M125_13TeV_%s_%sMiniTree.root", region[nw].data(), channel.data()));

    }
 
    RooDataSet dataSetDataSB("dataSetDataSB", "dataSetDataSB", variables, Cut(catCut && sbCut), WeightVar(evWeight), Import(*treeData));
    RooRealVar nSBDataEvents("nSBDataEvents", "nSBDataEvents", 0., 1.e10);

    nSBDataEvents.setVal(dataSetDataSB.sumEntries());
    nSBDataEvents.setConstant(true);

    if( !anotherModel ){ // using exponential function
    
      RooRealVar lamda("lamda", "lamda", -0.015, -0.04, -0.01);

      RooExponential model_mJetSB("model_mJetSB", "model_mJetSB", mJet, lamda);
      RooExtendPdf ext_model_mJetSB("ext_model_mJetSB", "ext_model_mJetSB", model_mJetSB, nSBDataEvents);
      
      ext_model_mJetSB.fitTo(dataSetDataSB, SumW2Error(true), Extended(true), Range("lowSB,highSB"), Strategy(2), Minimizer("Minuit2"), Save(1));

      // Normalize factor to normalize the background in signal region of data

      RooAbsReal* nSIGFit = ext_model_mJetSB.createIntegral(RooArgSet(mJet), NormSet(mJet), Range("signal"));
      RooAbsReal* nSBFit  = ext_model_mJetSB.createIntegral(RooArgSet(mJet), NormSet(mJet), Range("lowSB,highSB"));

      normFactor[nw] = dataSetDataSB.sumEntries()*(nSIGFit->getVal()/nSBFit->getVal());

      dataSetDataSB   .plotOn(mJetFrame, Binning(binsmJet), DataError(RooAbsData::SumW2), MarkerColor((nw==0)?kBlack:kRed), LineColor((nw==0)?kBlack:kRed));
      ext_model_mJetSB.plotOn(mJetFrame, Range("allRange"), LineColor((nw==0)?kGreen+1:kRed));

    }
    
    else{ // using error-exponential function
      
      RooRealVar constantSB("constantSB", "constantSB", -0.02,  -1.,   0.);
      RooRealVar widthSB   ("widthSB",    "widthSB",      100,   0., 200.);
      RooRealVar offsetSB  ("offsetSB",   "offsetSB",      30, -50., 200.);
      
      offsetSB.setConstant(true);
      
      RooErfExpPdf model_mJetSB("model_mJetSB", "model_mJetSB", mJet, constantSB, offsetSB, widthSB);
      RooExtendPdf ext_model_mJetSB("ext_model_mJetSB", "ext_model_mJetSB", model_mJetSB, nSBDataEvents);
      
      ext_model_mJetSB.fitTo(dataSetDataSB, SumW2Error(true), Extended(true), Range("lowSB,highSB"), Strategy(2), Minimizer("Minuit2"), Save(1));
      
      // Normalize factor to normalize the background in signal region of data
      
      RooAbsReal* nSIGFit = ext_model_mJetSB.createIntegral(RooArgSet(mJet), NormSet(mJet), Range("signal"));
      RooAbsReal* nSBFit  = ext_model_mJetSB.createIntegral(RooArgSet(mJet), NormSet(mJet), Range("lowSB,highSB"));
      
      normFactor[nw] = dataSetDataSB.sumEntries()*(nSIGFit->getVal()/nSBFit->getVal());
      
      dataSetDataSB   .plotOn(mJetFrame, Binning(binsmJet), DataError(RooAbsData::SumW2), MarkerColor((nw==0)?kBlack:kRed), LineColor((nw==0)?kBlack:kRed));
      ext_model_mJetSB.plotOn(mJetFrame, Range("allRange"), LineColor((nw==0)?kGreen+1:kRed));
      
    }

  } // end of weight loop 

  // Output the results

  TLatex lar;

  lar.SetTextSize(0.03);
  lar.SetLineWidth(5);

  TCanvas cv("cv","",0,0,1000,800);

  cv.cd()->SetLogy();
  mJetFrame->SetMinimum(1e-2);
  mJetFrame->SetMaximum(50);
  mJetFrame->SetTitle("");
  mJetFrame->Draw();
  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar.DrawLatexNDC(0.60, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s btag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, Form("Data side band: %s", study.data()));
  if( !anotherModel && removeMinor )
    lar.DrawLatexNDC(0.15, 0.78, Form("Norm factor: %.3f#pm%.3f", normFactor[0], TMath::Max(fabs(normFactor[2]-normFactor[0]),fabs(normFactor[1]-normFactor[0]))));
  else
    lar.DrawLatexNDC(0.15, 0.78, Form("Norm factor: %.3f", normFactor[0]));
  cv.Print(Form("jetEnScaleOnData_%s_%s_cat%s.pdf", study.data(), channel.data(), catcut.data()));

}
