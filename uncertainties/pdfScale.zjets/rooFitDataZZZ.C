R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/PdfDiagonalizer_cc.so)
using namespace RooFit;

void rooFitData(string channel, string catcut, string type, int first, int last){

  // Suppress all the INFO message

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooMsgService::instance().setSilentMode(true);

  // Input files and sum all backgrounds

  TChain* treeZjets = new TChain("tree");

  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-100to200_13TeV_toyMC.root", channel.data()));
  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-200to400_13TeV_toyMC.root", channel.data()));
  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-400to600_13TeV_toyMC.root", channel.data()));
  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-600toInf_13TeV_toyMC.root", channel.data()));

  // Define all the variables from the trees 

  RooRealVar cat("cat", "", 0, 2);
  RooRealVar mJet("prmass", "",  30.,  300.);
  RooRealVar mZH("mllbb", "M_{ZH}", 800., 4000., "GeV");
  RooBinning mZHbin(32, 800., 4000.);

  mZH.setRange("fullRange", 800., 4000.);

  TCut catCut = Form("cat==%s", catcut.c_str());
  TCut sbCut  = "prmass>30 && !(prmass>65 && prmass<135) && prmass<300";
  TCut sigCut = "prmass>105 && prmass<135";

  RooPlot* mZHsbFrame = mZH.frame();
  RooPlot* mZHsgFrame = mZH.frame();

  int N = last-first;
  int iw = N-1;
  float alphaScale[11][N], alphaCentral[11];

  TF1* f_alpha = new TF1("f_alpha", "TMath::Exp([0]*x+[1]/x)/TMath::Exp([2]*x+[3]/x)", 800, 4000);
   
  f_alpha->SetTitle("");
  f_alpha->GetXaxis()->SetTitle("m_{ZH} (GeV)");
  f_alpha->GetYaxis()->SetTitle("#alpha Ratio");
  f_alpha->GetYaxis()->SetTitleOffset(1.3);
  
  TCanvas* cv = new TCanvas("cv", "", 0, 0, 1000, 800);

  cv->cd();
  
  for( int nw = last; nw >= first; --nw ){
    
    RooRealVar evWeight(Form("evweight%02i",nw), "", -1.e10, 1.e10);
    RooArgSet variables(cat, mJet, mZH, evWeight);

    // Create a dataset from a tree -> to process unbinned likelihood fitting

    RooDataSet dataSetZjetsSB("dataSetZjetsSB", "dataSetZjetsSB", variables, Cut(catCut && sbCut),  WeightVar(evWeight), Import(*treeZjets));  
    RooDataSet dataSetZjetsSG("dataSetZjetsSG", "dataSetZjetsSG", variables, Cut(catCut && sigCut), WeightVar(evWeight), Import(*treeZjets));
  
    // Total event numbers

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

    ext_model_ZHSB.fitTo(dataSetZjetsSB, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));

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
    
    ext_model_ZHSG.fitTo(dataSetZjetsSG, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));

    float p2 = c.getVal();
    float p3 = d.getVal();

    // Set the model of alpha ratio

    f_alpha->SetParameters(p2,p3,p0,p1);
    f_alpha->SetMinimum(0);
    f_alpha->SetMaximum( (channel=="ele"&&catcut=="1") ? 50 : ( (channel=="mu"&&catcut=="1") ? 1.5 : 0.2 ) );
    f_alpha->SetLineColor((nw==first)?kBlue:kCyan);
    f_alpha->DrawCopy((nw==last) ? "" : "same");

    int mzh = 800;

    for( int im = 0; im < 11; ++im ){

      if( nw != first )
	alphaScale[im][iw] = f_alpha->Eval(mzh);
      else
	alphaCentral[im] = f_alpha->Eval(mzh);

      mzh += (mzh<2000) ? 200 : 500;

    }

    --iw;

    fprintf(stdout, "weight=%i\tp0=%f\tp1=%f\tp2=%f\tp3=%f\n", nw, p0, p1, p2, p3);

    // Plot the results to a frame 

    dataSetZjetsSB.plotOn(mZHsbFrame, Binning(mZHbin), MarkerColor((nw==first)?kBlue:kCyan), LineColor((nw==first)?kBlue:kCyan));
    model_ZHSB.plotOn(mZHsbFrame, Range("fullRange"), LineColor((nw==first)?kBlue:kCyan));

    dataSetZjetsSG.plotOn(mZHsgFrame, Binning(mZHbin), MarkerColor((nw==first)?kBlue:kCyan), LineColor((nw==first)?kBlue:kCyan));
    model_ZHSG.plotOn(mZHsgFrame, Range("fullRange"), LineColor((nw==first)?kBlue:kCyan));

  } // end of weight loop

  TLatex* lar = new TLatex();

  lar->SetTextSize(0.035);
  lar->SetLineWidth(5);

  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.72, 0.80, Form("%s  %s btag", channel.data(), catcut.data()));
  cv->Print(Form("alpha_%sScale_%s_cat%s.pdf(", type.data(), channel.data(), catcut.data()));

  // Calculate RMS value of each mass bin

  float Unc[11], Mzh[11] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  for( int im = 0; im < 11; ++im ){

    if( type == "mur1" )
      Unc[im] = ( fabs(alphaScale[im][1]-alphaCentral[im]) > fabs(alphaScale[im][0]-alphaCentral[im]) ) ?
	fabs(alphaScale[im][1]-alphaCentral[im]) : fabs(alphaScale[im][0]-alphaCentral[im]);

    else 
      Unc[im] = TMath::RMS(N, alphaScale[im]);

  } // end of mass points

  TGraphErrors *g_alpha = new TGraphErrors(11, Mzh, alphaCentral, 0, Unc);

  g_alpha->SetTitle("");
  g_alpha->GetXaxis()->SetTitle("m_{ZH} (GeV)");
  g_alpha->GetYaxis()->SetTitle("#alpha Ratio");  
  g_alpha->GetYaxis()->SetTitleOffset(1.3);
  g_alpha->GetXaxis()->SetLimits(800,4000);
  g_alpha->SetMinimum(0);
  g_alpha->SetMaximum( (channel=="ele"&&catcut=="1") ? 50 : ( (channel=="mu"&&catcut=="1") ? 1.5 : 0.2 ) );
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
  
  cv->Clear();
  cv->cd();
  g_alpha->Draw("apz");
  g_alpha->Draw("3same");
  g_alpha->Draw("cxsame");
  leg->AddEntry(g_alpha, "alpha ratio with uncertainties", "lf");
  leg->Draw();
  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.72, 0.80, Form("%s  %s btag", channel.data(), catcut.data()));
  cv->Print(Form("alpha_%sScale_%s_cat%s.pdf", type.data(), channel.data(), catcut.data()));  

  TLegend* leg1 = new TLegend(0.55,0.55,0.85,0.70);

  leg1->AddEntry(mZHsbFrame->findObject(mZHsbFrame->nameOf(0)), "MC side band (central)", "lep");
  leg1->AddEntry(mZHsbFrame->findObject(mZHsbFrame->nameOf(1)), "Fit curve (central)", "l");
  leg1->Draw();

  cv->Clear();
  cv->cd();
  mZHsbFrame->addObject(leg1);
  mZHsbFrame->SetTitle("");
  mZHsbFrame->SetMinimum(0);
  mZHsbFrame->Draw();
  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.72, 0.80, Form("%s  %s btag", channel.data(), catcut.data()));
  cv->Print(Form("alpha_%sScale_%s_cat%s.pdf", type.data(), channel.data(), catcut.data()));

  TLegend* leg2 = new TLegend(0.55,0.55,0.85,0.70);

  leg2->AddEntry(mZHsgFrame->findObject(mZHsgFrame->nameOf(0)), "MC signal region (central)", "lep");
  leg2->AddEntry(mZHsgFrame->findObject(mZHsgFrame->nameOf(1)), "Fit curve (central)", "l");
  leg2->Draw();

  cv->Clear();
  cv->cd();
  mZHsgFrame->addObject(leg2);
  mZHsgFrame->SetTitle("");
  mZHsgFrame->SetMinimum(0);
  mZHsgFrame->Draw();
  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.72, 0.80, Form("%s  %s btag", channel.data(), catcut.data()));
  cv->Print(Form("alpha_%sScale_%s_cat%s.pdf)", type.data(), channel.data(), catcut.data()));

}
