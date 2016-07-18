R__LOAD_LIBRARY(PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(PDFs/PdfDiagonalizer_cc.so)
using namespace RooFit;

void rooFitData(string channel, string catcut){

  // Suppress all the INFO message

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooMsgService::instance().setSilentMode(true);

  // Input files and sum all backgrounds

  TChain* treeZjets = new TChain("tree");

  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-100to200_13TeV_toyMC.root", channel.data()));
  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-200to400_13TeV_toyMC.root", channel.data()));
  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-400to600_13TeV_toyMC.root", channel.data()));
  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-600toInf_13TeV_toyMC.root", channel.data()));

  RooRealVar mzh("mzh", "M_{ZH}", 900., 3000., "GeV");
  RooPlot* alphaFrame = mzh.frame();
  
  for(int js = 2; js >= 0; --js){

    // Define all the variables from the trees 

    RooRealVar cat ("cat", "", 0, 2);
    RooRealVar mJet("prmass", "M_{jet}",  30.,  300., "GeV");
    RooRealVar mZH (Form("mllbb%i",js), "M_{ZH}", 900., 3000., "GeV");
    RooRealVar evWeight("evweight", "", -1.e10, 1.e10);

    // Set the range in zh mass 

    mZH.setRange("fullRange", 900., 3000.);

    RooBinning binsmZH(21, 900, 3000);

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

    RooRealVar a("a", "a", -0.0005, -1.,     1.);
    RooRealVar b("b", "b",    1200,  0., 10000.);
  
    RooGenericPdf model_ZHSB("model_ZHSB", "model_ZHSB", "TMath::Exp(@1*@0+@2/@0)", RooArgSet(mZH,a,b));
    RooGenericPdf model_ZHSG("model_ZHSG", "model_ZHSG", "TMath::Exp(@1*@0+@2/@0)", RooArgSet(mZH,a,b));

    RooExtendPdf ext_model_ZHSB("ext_model_ZHSB", "ext_model_ZHSB", model_ZHSB, nSBMcEvents);
    RooExtendPdf ext_model_ZHSG("ext_model_ZHSG", "ext_model_ZHSG", model_ZHSG, nSGMcEvents);

    // Fit ZH mass in side band  

    RooFitResult* mZHSB_result = ext_model_ZHSB.fitTo(dataSetZjetsSB, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));

    RooAbsReal* nZHSBFit = ext_model_ZHSB.createIntegral(RooArgSet(mZH), NormSet(mZH), Range("fullRange"));

    float p0 = a.getVal();
    float p1 = b.getVal();

    // Fit ZH mass in signal region

    RooFitResult* mZHSG_result = ext_model_ZHSG.fitTo(dataSetZjetsSG, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));

    RooAbsReal* nZHSGFit = ext_model_ZHSG.createIntegral(RooArgSet(mZH), NormSet(mZH), Range("fullRange"));

    float p2 = a.getVal();
    float p3 = b.getVal();

    // Draw the model of alpha ratio

    RooGenericPdf model_alpha("model_alpha", "model_alpha", Form("TMath::Exp(%f*@0+%f/@0)/TMath::Exp(%f*@0+%f/@0)", p2,p3,p0,p1), RooArgSet(mzh));

    // Plot the results to a frame 

    model_alpha.plotOn(alphaFrame, LineColor((js==0)?kBlue:kYellow));

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

  TLatex* lar = new TLatex();

  lar->SetTextSize(0.035);
  lar->SetLineWidth(5);
  
  TCanvas* c = new TCanvas("c","",0,0,1000,800);

  c->cd();
  alphaFrame->Draw();
  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.75, 0.80, Form("%s  %s btag", channel.data(), catcut.data()));
  lar->DrawLatexNDC(0.75, 0.75, "JES");
  c->Print(Form("alpha_jetEnScale_%s_cat%s.pdf", channel.data(), catcut.data()));

}
