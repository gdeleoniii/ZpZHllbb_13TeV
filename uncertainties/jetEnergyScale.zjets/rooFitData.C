R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/PdfDiagonalizer_cc.so)
using namespace RooFit;

void rooFitData(string channel, string catcut){

  // Suppress all the INFO message

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooMsgService::instance().setSilentMode(true);

  string region[3] = {"central","up","down"};
  RooRealVar mZH("mllbb", "M_{ZH}", 900., 3000., "GeV");

  RooPlot* alphaFrame = mZH.frame();
  RooPlot* mZHsbFrame = mZH.frame();
  RooPlot* mZHsgFrame = mZH.frame();

  RooBinning binsmZH(21, 900, 3000);

  for(int nw = 2; nw >= 0; --nw){

    //fprintf(stdout, "Using weight %i\n", nw);

    // Input files and sum all backgrounds

    TChain* treeZjets = new TChain("tree");
    
    treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-100to200_13TeV_%s_toyMC.root", channel.data(), region[nw].data()));
    treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-200to400_13TeV_%s_toyMC.root", channel.data(), region[nw].data()));
    treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-400to600_13TeV_%s_toyMC.root", channel.data(), region[nw].data()));
    treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-600toInf_13TeV_%s_toyMC.root", channel.data(), region[nw].data()));

    // Define all the variables from the trees 

    RooRealVar cat ("cat", "", 0, 2);
    RooRealVar mJet("prmass", "", 30., 300.);
    RooRealVar evWeight("evweight", "", -1.e10, 1.e10);

    // Set the range in zh mass 

    mZH.setRange("fullRange", 900., 3000.);

    RooArgSet variables(cat, mJet, mZH, evWeight);

    TCut catCut = Form("cat==%s", catcut.c_str());
    TCut sbCut  = "prmass>30 && !(prmass>65 && prmass<135) && prmass<300";
    TCut sigCut = "prmass>105 && prmass<135";

    // Create a dataset from a tree -> to process an unbinned likelihood fitting

    RooDataSet dataSetZjetsSB("dataSetZjetsSB", "dataSetZjetsSB", variables, Cut(catCut && sbCut),  WeightVar(evWeight), Import(*treeZjets));  
    RooDataSet dataSetZjetsSG("dataSetZjetsSG", "dataSetZjetsSG", variables, Cut(catCut && sigCut), WeightVar(evWeight), Import(*treeZjets));
  
    // Total events number

    float totalSBMcEv = dataSetZjetsSB.sumEntries();
    float totalSGMcEv = dataSetZjetsSG.sumEntries();

    RooRealVar nSBMcEvents("nSBMcEvents", "nSBMcEvents", 0., 9999999.);
    RooRealVar nSGMcEvents("nSGMcEvents", "nSGMcEvents", 0., 9999999.);

    nSBMcEvents.setVal(totalSBMcEv);
    nSBMcEvents.setConstant(true);
  
    nSGMcEvents.setVal(totalSGMcEv);
    nSGMcEvents.setConstant(true);
  
    // Alpha ratio part

    // Fit ZH mass in side band 

    float bmin, bmax;

    if( channel == "ele" ){
      bmin = (catcut=="1") ? 1300. : 0.;
      bmax = (catcut=="1") ? 1800. : 2.;
    }

    else if( channel == "mu" ){
      bmin = 600.; 
      bmax = 1200.;
    }
    
    RooRealVar a("a", "a", -0.002, -0.005, 0.);
    RooRealVar b("b", "b", (bmin+bmax)*0.5, bmin, bmax);

    RooGenericPdf model_ZHSB("model_ZHSB", "model_ZHSB", "TMath::Exp(@1*@0+@2/@0)", RooArgSet(mZH,a,b));
    RooExtendPdf ext_model_ZHSB("ext_model_ZHSB", "ext_model_ZHSB", model_ZHSB, nSBMcEvents);

    RooFitResult* mZHSB_result = ext_model_ZHSB.fitTo(dataSetZjetsSB, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));
    RooAbsReal* nZHSBFit = ext_model_ZHSB.createIntegral(RooArgSet(mZH), NormSet(mZH), Range("fullRange"));

    float p0 = a.getVal();
    float p1 = b.getVal();

    // Fit ZH mass in signal region

    float dmin, dmax;
    
    if( channel == "ele" ){
      dmin = (catcut=="1") ? 800.  : 2100.;
      dmax = (catcut=="1") ? 1400. : 2800.;
    }
    
    else if( channel == "mu" ){
      dmin = (catcut=="1") ? 0. : 1500.;
      dmax = (catcut=="1") ? 1. : 2500.;
    }
    
    RooRealVar c("c", "c", -0.002, -0.005, 0.);
    RooRealVar d("d", "d", (dmin+dmax)*0.5, dmin, dmax);

    RooGenericPdf model_ZHSG("model_ZHSG", "model_ZHSG", "TMath::Exp(@1*@0+@2/@0)", RooArgSet(mZH,c,d));
    RooExtendPdf ext_model_ZHSG("ext_model_ZHSG", "ext_model_ZHSG", model_ZHSG, nSGMcEvents);

    RooFitResult* mZHSG_result = ext_model_ZHSG.fitTo(dataSetZjetsSG, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));
    RooAbsReal* nZHSGFit = ext_model_ZHSG.createIntegral(RooArgSet(mZH), NormSet(mZH), Range("fullRange"));

    float p2 = c.getVal();
    float p3 = d.getVal();

    // Draw the model of alpha ratio

    RooGenericPdf model_alpha("model_alpha", "model_alpha", Form("TMath::Exp(%f*@0+%f/@0)/TMath::Exp(%f*@0+%f/@0)", p2,p3,p0,p1), RooArgSet(mZH));

    fprintf(stdout, "p0=%f\tp1=%f\tp2=%f\tp3=%f\n", p0,p1,p2,p3);

    // Plot the results to a frame 

    model_alpha.plotOn(alphaFrame, LineColor((nw==0)?kBlue:kYellow));

    dataSetZjetsSB.plotOn(mZHsbFrame, Binning(binsmZH), MarkerColor((nw==0)?kBlue:kCyan), LineColor((nw==0)?kBlue:kCyan));
    model_ZHSB.plotOn(mZHsbFrame, Range("fullRange"), LineColor((nw==0)?kBlue:kCyan));
    
    dataSetZjetsSG.plotOn(mZHsgFrame, Binning(binsmZH), MarkerColor((nw==0)?kBlue:kCyan), LineColor((nw==0)?kBlue:kCyan));
    model_ZHSG.plotOn(mZHsgFrame, Range("fullRange"), LineColor((nw==0)?kBlue:kCyan));

  }

  TLegend* leg = new TLegend(0.15,0.15,0.30,0.25);

  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(alphaFrame->findObject(alphaFrame->nameOf(0)), "central values", "l");
  leg->Draw();

  alphaFrame->addObject(leg);
  alphaFrame->SetTitle("");
  alphaFrame->SetMaximum(0.03);
  alphaFrame->GetYaxis()->SetTitle("#alpha Ratio");
  alphaFrame->GetYaxis()->SetTitleOffset(1.3);

  mZHsbFrame->SetTitle("");
  mZHsbFrame->SetMinimum(0);

  mZHsgFrame->SetTitle("");
  mZHsgFrame->SetMinimum(0);

  TLatex* lar = new TLatex();

  lar->SetTextSize(0.035);
  lar->SetLineWidth(5);
  
  TCanvas* c = new TCanvas("c","",0,0,1000,800);

  c->cd();
  alphaFrame->Draw();
  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.72, 0.80, Form("%s  %s btag", channel.data(), catcut.data()));
  lar->DrawLatexNDC(0.72, 0.75, "JES");
  c->Print(Form("alpha_jetEnScale_%s_cat%s.pdf(", channel.data(), catcut.data()));

  c->cd();
  mZHsbFrame->Draw();
  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.72, 0.80, Form("%s  %s btag", channel.data(), catcut.data()));
  lar->DrawLatexNDC(0.72, 0.75, "side band");
  c->Print(Form("alpha_jetEnScale_%s_cat%s.pdf", channel.data(), catcut.data()));

  c->cd();
  mZHsgFrame->Draw();
  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.72, 0.80, Form("%s  %s btag", channel.data(), catcut.data()));
  lar->DrawLatexNDC(0.72, 0.75, "signal region");
  c->Print(Form("alpha_jetEnScale_%s_cat%s.pdf)", channel.data(), catcut.data()));

}
