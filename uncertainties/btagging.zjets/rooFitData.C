R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/PdfDiagonalizer_cc.so)
using namespace RooFit;

void rooFitData(string channel, string catcut){

  // Suppress all the INFO message

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooMsgService::instance().setSilentMode(true);

  // Define all the variables from the trees

  RooRealVar cat ("cat", "", 0, 2);
  RooRealVar mJet("prmass", "", 30., 300.);
  RooRealVar evWeight("evweight", "", 0., 1.e10);
  RooRealVar mZH("mllbb", "M_{ZH}", 900., 3000., "GeV");
 
  mZH.setRange("fullRange", 900., 3000.);
 
  RooBinning mZHbin(21, 900., 3000.);
  RooArgSet variables(cat, mJet, mZH, evWeight);
 
  TCut catCut = Form("cat==%s", catcut.c_str());
  TCut sbCut  = "prmass>30 && !(prmass>65 && prmass<135) && prmass<300";
  TCut sigCut = "prmass>105 && prmass<135";

  string region[3] = {"central","up","down"};
  float alpha[11][3];

  TF1* f_alpha = new TF1("f_alpha", "TMath::Exp([0]*x+[1]/x)/TMath::Exp([2]*x+[3]/x)", 800, 4000);   

  for(int nw = 2; nw >= 0; --nw){

    // Input files and sum all backgrounds

    TChain* treeZjets = new TChain("tree");
    
    treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-100to200_13TeV_%s_toyMC.root", channel.data(), region[nw].data()));
    treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-200to400_13TeV_%s_toyMC.root", channel.data(), region[nw].data()));
    treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-400to600_13TeV_%s_toyMC.root", channel.data(), region[nw].data()));
    treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-600toInf_13TeV_%s_toyMC.root", channel.data(), region[nw].data()));
    
    // Create a dataset from a tree -> to process an unbinned likelihood fitting

    RooDataSet dataSetZjetsSB("dataSetZjetsSB", "dataSetZjetsSB", variables, Cut(catCut && sbCut),  WeightVar(evWeight), Import(*treeZjets));  
    RooDataSet dataSetZjetsSG("dataSetZjetsSG", "dataSetZjetsSG", variables, Cut(catCut && sigCut), WeightVar(evWeight), Import(*treeZjets));

    // Total events number

    RooRealVar nSBMcEvents("nSBMcEvents", "nSBMcEvents", 0., 1.e10);
    RooRealVar nSGMcEvents("nSGMcEvents", "nSGMcEvents", 0., 1.e10);

    nSBMcEvents.setVal(dataSetZjetsSB.sumEntries());
    nSBMcEvents.setConstant(true);
  
    nSGMcEvents.setVal(dataSetZjetsSG.sumEntries());
    nSGMcEvents.setConstant(true);
  
    // Alpha ratio part

    // Fit ZH mass in side band 

    float bmin, bmax;

    if( channel == "ele" ){
      bmin = (catcut=="1") ?  700. : 1500.;
      bmax = (catcut=="1") ? 1100. : 2500.;
    }

    else if( channel == "mu" ){
      bmin = (catcut=="1") ? 200. : 2100.; 
      bmax = (catcut=="1") ? 700. : 2700.;
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
      dmin = (catcut=="1") ? 3200. : 0.;
      dmax = (catcut=="1") ? 3500. : 1.;
    }
    
    else if( channel == "mu" ){
      dmin = (catcut=="1") ? 0. : 10.;
      dmax = (catcut=="1") ? 1. : 15.;
    }
    
    RooRealVar c("c", "c", -0.002, -0.005, 0.);
    RooRealVar d("d", "d", (dmin+dmax)*0.5, dmin, dmax);

    RooGenericPdf model_ZHSG("model_ZHSG", "model_ZHSG", "TMath::Exp(@1*@0+@2/@0)", RooArgSet(mZH,c,d));
    RooExtendPdf ext_model_ZHSG("ext_model_ZHSG", "ext_model_ZHSG", model_ZHSG, nSGMcEvents);

    RooFitResult* mZHSG_result = ext_model_ZHSG.fitTo(dataSetZjetsSG, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));
    RooAbsReal* nZHSGFit = ext_model_ZHSG.createIntegral(RooArgSet(mZH), NormSet(mZH), Range("fullRange"));

    float p2 = c.getVal();
    float p3 = d.getVal();

    // Set the model of alpha ratio

    f_alpha->SetParameters(p2,p3,p0,p1);

    int mzh = 800;
    for( int im = 0; im < 11; ++im ){
      alpha[im][nw] = f_alpha->Eval(mzh);
      mzh += (mzh<2000) ? 200 : 500;
    }

    fprintf(stdout, "p0=%f\tp1=%f\tp2=%f\tp3=%f\n", p0,p1,p2,p3);

    delete treeZjets;

  } // end of weight loop
   
  // Calculate uncertainty of each mass bin

  float Alpha[11], Unc[11], Mzh[11] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  for( int im = 0; im < 11; ++im ){

    Alpha[im] = alpha[im][0];
    Unc[im] = ( fabs(alpha[im][1]-alpha[im][0]) > fabs(alpha[im][2]-alpha[im][0]) ) ?
      fabs(alpha[im][1]-alpha[im][0]) : fabs(alpha[im][2]-alpha[im][0]);

  } // end of mass points
  
  TGraphErrors *g_alpha = new TGraphErrors(11, Mzh, Alpha, 0, Unc);

  g_alpha->SetTitle("");
  g_alpha->GetXaxis()->SetTitle("m_{ZH} (GeV)");
  g_alpha->GetYaxis()->SetTitle("#alpha Ratio");  
  g_alpha->GetYaxis()->SetTitleOffset(1.3);
  g_alpha->GetXaxis()->SetLimits(800,4000);
  g_alpha->SetMinimum(0);
  g_alpha->SetLineWidth(2);
  g_alpha->SetLineColor(kBlue);
  g_alpha->SetMarkerStyle(8);
  g_alpha->SetMarkerColor(kBlue);
  g_alpha->SetFillStyle(1001);
  g_alpha->SetFillColor(kYellow);  
  
  TLegend* leg = new TLegend(0.15,0.15,0.35,0.25);

  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);

  TLatex* lar = new TLatex();
  
  lar->SetTextSize(0.035);
  lar->SetLineWidth(5);

  TCanvas* cv = new TCanvas("cv", "", 0, 0, 1000, 800);

  cv->cd();
  g_alpha->Draw("apz");
  g_alpha->Draw("3same");
  g_alpha->Draw("cxsame");
  leg->AddEntry(g_alpha, "alpha ratio with uncertainties", "lf");
  leg->Draw();
  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.72, 0.80, Form("%s  %s btag", channel.data(), catcut.data()));
  cv->Print(Form("alpha_bTagScale_%s_cat%s.pdf", channel.data(), catcut.data()));

}
