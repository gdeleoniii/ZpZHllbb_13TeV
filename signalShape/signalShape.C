R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/PdfDiagonalizer_cc.so)
#include "zpEleShape.h"
#include "zpMuShape.h"
using namespace RooFit;

void signalShape(string chan, int cat){

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooMsgService::instance().setSilentMode(true);

  float mzh[10] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500};

  RooRealVar mZH("mZH", "m_{ZH} (GeV)", 500, 4000);
  mZH.setRange("fullRange", 500., 4000.);

  RooPlot* myFrame = mZH.frame();
  
  TH1F* f_ = NULL;
  TLegend* leg = new TLegend(0.65,0.40,0.85,0.75);
  
  for( int i = 0; i < 10; ++i ){

    RooRealVar m("m", "mean", mzh[i], mzh[i]-50, mzh[i]+50);
    RooRealVar s("s", "sigma", 25., 80.);
    RooRealVar a("a", "alpha", 1., 0.1, 3.);
    RooRealVar n("n", "n", 80., 0., 160.);

    m.setConstant(true);

    f_ = (chan == "mu") ?
      zpMuShape(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",(int)mzh[i]), cat, (int)mzh[i]) :
      zpEleShape(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",(int)mzh[i]), cat, (int)mzh[i]);

    RooDataHist h_("h_", "", mZH, Import(*f_));
    RooCBShape model("model", "Cystal Ball Function", mZH, m, s, a, n);

    int color = 100-4*i;

    mZH.setRange("range", mzh[i]-400, mzh[i]+300);
    model.fitTo(h_, SumW2Error(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));
    h_.plotOn(myFrame,MarkerColor(color),LineColor(color));
    model.plotOn(myFrame,LineColor(color));
    leg->AddEntry(myFrame->findObject(myFrame->nameOf(1)), Form("M_{ZH} = %d GeV", (int)mzh[i]), "l");
  
    fprintf(stdout, "m=%f\ts=%f\ta=%f\tn=%f\n", m.getVal(), s.getVal(), a.getVal(), n.getVal());
      
  }

  leg->Draw();
  myFrame->addObject(leg);
  myFrame->SetTitle("");
  myFrame->SetMinimum(1.e-3);

  TLatex lar;

  lar.SetTextSize(0.04);
  lar.SetLineWidth(5);

  TCanvas c("c", "", 0, 0, 800, 600);

  c.cd();
  c.SetLogy();
  myFrame->Draw();
  lar.DrawLatexNDC(0.20, 0.83, "CMS");
  lar.DrawLatexNDC(0.20, 0.79, "#it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.72, 0.79, Form("%s  %i btag", chan.data(), cat));
  c.Print(Form("signalShape_%s_%ibtag.pdf", chan.data(), cat));

}
