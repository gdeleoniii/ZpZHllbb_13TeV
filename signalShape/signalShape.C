#include "zpEleShape.h"
#include "zpMuShape.h"
R__LOAD_LIBRARY(PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(PDFs/PdfDiagonalizer_cc.so)
using namespace RooFit;

const int N = 11;

void signalShape(){

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  //RooMsgService::instance().setSilentMode(true);

  float mzh[N] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};

  RooRealVar mZH("mZH", "m_{ZH} (GeV)", 500, 4500);

  RooPlot* mu1btagFrame = mZH.frame();
  RooPlot* mu2btagFrame = mZH.frame();
  RooPlot* ele1btagFrame = mZH.frame();
  RooPlot* ele2btagFrame = mZH.frame();
  
  TLegend* leg0 = new TLegend(0.65,0.40,0.85,0.75);
  TLegend* leg1 = new TLegend(0.65,0.40,0.85,0.75);
  TLegend* leg2 = new TLegend(0.65,0.40,0.85,0.75);
  TLegend* leg3 = new TLegend(0.65,0.40,0.85,0.75);

  for( int i = 0; i < N; ++i ){

    RooRealVar m("m", "mean",   mzh[i], 800., 4000.);
    RooRealVar s("s", "sigma",            0.,   50.);
    RooRealVar a("a", "alpha",     0.5, -50.,   50.);
    RooRealVar n("n", "n",        100.,   0.,  200.);

    TH1F* mu_1btag = zpMuShape(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",(int)mzh[i]),1);
    RooDataHist mu1btag("mu1btag", "", mZH, Import(*mu_1btag));
    RooCBShape cbpdfm1("cbpdfm1", "Cystal Ball Function", mZH, m, s, a, n);
    cbpdfm1.fitTo(mu1btag);
    mu1btag.plotOn(mu1btagFrame);
    cbpdfm1.plotOn(mu1btagFrame,LineColor(100-2*i));
    leg0->AddEntry(mu1btagFrame->findObject(mu1btagFrame->nameOf(1)), Form("M_{ZH} = %d GeV", (int)mzh[i]), "l");
    
    TH1F* mu_2btag = zpMuShape(Form("/data7/htong/skim_NCUGlobalTuples/skim_mu_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",(int)mzh[i]),2);
    RooDataHist mu2btag("mu2btag", "", mZH, Import(*mu_2btag));    
    RooCBShape cbpdfm2("cbpdfm2", "Cystal Ball Function", mZH, m, s, a, n);
    cbpdfm2.fitTo(mu2btag);
    mu2btag.plotOn(mu2btagFrame);
    cbpdfm2.plotOn(mu2btagFrame,LineColor(100-2*i));
    leg1->AddEntry(mu2btagFrame->findObject(mu2btagFrame->nameOf(1)), Form("M_{ZH} = %d GeV", (int)mzh[i]), "l");

    TH1F* ele_1btag = zpEleShape(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",(int)mzh[i]),1);
    RooDataHist ele1btag("ele1btag", "", mZH, Import(*ele_1btag));
    RooCBShape cbpdfe1("cbpdfe1", "Cystal Ball Function", mZH, m, s, a, n);
    cbpdfe1.fitTo(ele1btag);
    ele1btag.plotOn(ele1btagFrame);
    cbpdfe1.plotOn(ele1btagFrame,LineColor(100-2*i));
    leg2->AddEntry(ele1btagFrame->findObject(ele1btagFrame->nameOf(1)), Form("M_{ZH} = %d GeV", (int)mzh[i]), "l");

    TH1F* ele_2btag = zpEleShape(Form("/data7/htong/skim_NCUGlobalTuples/skim_ele_crab_ZprimeToZhToZlephbb_narrow_M-%d_13TeV-madgraph.root",(int)mzh[i]),2);
    RooDataHist ele2btag("ele2btag", "", mZH, Import(*ele_2btag));
    RooCBShape cbpdfe2("cbpdfe2", "Cystal Ball Function", mZH, m, s, a, n);
    cbpdfe2.fitTo(ele2btag);
    ele2btag.plotOn(ele2btagFrame);
    cbpdfe2.plotOn(ele2btagFrame,LineColor(100-2*i));
    leg3->AddEntry(ele2btagFrame->findObject(ele2btagFrame->nameOf(1)), Form("M_{ZH} = %d GeV", (int)mzh[i]), "l");
    
  }

  leg0->Draw();
  leg1->Draw();
  leg2->Draw();
  leg3->Draw();

  mu1btagFrame->addObject(leg0);
  mu2btagFrame->addObject(leg1);
  ele1btagFrame->addObject(leg2);
  ele2btagFrame->addObject(leg3);

  mu1btagFrame->SetTitle("");
  mu2btagFrame->SetTitle("");
  ele1btagFrame->SetTitle("");
  ele2btagFrame->SetTitle("");

  TLatex lar;

  lar.SetTextSize(0.04);
  lar.SetLineWidth(5);

  TCanvas c("c", "", 0, 0, 800, 600);

  c.cd();
  c.SetLogy();
  mu1btagFrame->Draw();
  lar.DrawLatexNDC(0.20, 0.83, "CMS");
  lar.DrawLatexNDC(0.20, 0.79, "#it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.72, 0.79, "mu  1 btag");
  c.Print("signalShape.pdf(");

  c.cd();
  c.SetLogy();
  mu2btagFrame->Draw();
  lar.DrawLatexNDC(0.20, 0.83, "CMS");
  lar.DrawLatexNDC(0.20, 0.79, "#it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.72, 0.79, "mu  2 btag");
  c.Print("signalShape.pdf");

  c.cd();
  c.SetLogy();
  ele1btagFrame->Draw();
  lar.DrawLatexNDC(0.20, 0.83, "CMS");
  lar.DrawLatexNDC(0.20, 0.79, "#it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.72, 0.79, "ele  1 btag");
  c.Print("signalShape.pdf");
  
  c.cd();
  c.SetLogy();
  ele2btagFrame->Draw();
  lar.DrawLatexNDC(0.20, 0.83, "CMS");
  lar.DrawLatexNDC(0.20, 0.79, "#it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.72, 0.79, "ele  2 btag");
  c.Print("signalShape.pdf)");

}
