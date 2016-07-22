#include <iostream>
#include <TFile.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include "../readHists.h"

void bTagEff(string channel, string flavor){

  readHist zjets100(Form("%s/Zjets/DYJetsToLL_M-50_HT-100to200_13TeV_%sMCbtagEff.root", channel.data(), channel.data()));
  readHist zjets200(Form("%s/Zjets/DYJetsToLL_M-50_HT-200to400_13TeV_%sMCbtagEff.root", channel.data(), channel.data()));
  readHist zjets400(Form("%s/Zjets/DYJetsToLL_M-50_HT-400to600_13TeV_%sMCbtagEff.root", channel.data(), channel.data()));
  readHist zjets600(Form("%s/Zjets/DYJetsToLL_M-50_HT-600toInf_13TeV_%sMCbtagEff.root", channel.data(), channel.data()));

  string ptnoCSVName, ptwtCSVName;

  if( flavor == "udsg" ){
    ptnoCSVName = "lJetPtnoCSV";
    ptwtCSVName = "lJetPtwtCSV";
  }

  else if( flavor == "c" ){
    ptnoCSVName = "cJetPtnoCSV";
    ptwtCSVName = "cJetPtwtCSV";
  }

  else if( flavor == "b" ){
    ptnoCSVName = "bJetPtnoCSV";
    ptwtCSVName = "bJetPtwtCSV";
  }

  float varBinslight[] = {30,50,70,100,140,200,300,670,1000,2000};
  float varBins[] = {30,50,70,100,140,200,300,670,2000};
  
  int nvarBinslight = sizeof(varBinslight)/sizeof(varBinslight[0])-1;
  int nvarBins = sizeof(varBins)/sizeof(varBins[0])-1;

  TH1F* h_jetPtnoCSV = new TH1F("h_jetPtnoCSV", "", (flavor=="udsg")?nvarBinslight:nvarBins, (flavor=="udsg")?varBinslight:varBins);
  TH1F* h_jetPtwtCSV = new TH1F("h_jetPtwtCSV", "", (flavor=="udsg")?nvarBinslight:nvarBins, (flavor=="udsg")?varBinslight:varBins);

  h_jetPtnoCSV->Add((TH1F*)(zjets100.getHist(ptnoCSVName.data()))->Clone("h_jetPtnoCSV"));
  h_jetPtnoCSV->Add((TH1F*)(zjets200.getHist(ptnoCSVName.data()))->Clone("h_jetPtnoCSV"));
  h_jetPtnoCSV->Add((TH1F*)(zjets400.getHist(ptnoCSVName.data()))->Clone("h_jetPtnoCSV"));
  h_jetPtnoCSV->Add((TH1F*)(zjets600.getHist(ptnoCSVName.data()))->Clone("h_jetPtnoCSV"));

  h_jetPtwtCSV->Add((TH1F*)(zjets100.getHist(ptwtCSVName.data()))->Clone("h_jetPtwtCSV"));
  h_jetPtwtCSV->Add((TH1F*)(zjets200.getHist(ptwtCSVName.data()))->Clone("h_jetPtwtCSV"));
  h_jetPtwtCSV->Add((TH1F*)(zjets400.getHist(ptwtCSVName.data()))->Clone("h_jetPtwtCSV"));
  h_jetPtwtCSV->Add((TH1F*)(zjets600.getHist(ptwtCSVName.data()))->Clone("h_jetPtwtCSV"));

  TFile f(Form("%s_%sflavor_zjetsBtagEff.root", channel.data(),flavor.data()),"recreate");

  TGraphAsymmErrors *g = new TGraphAsymmErrors();

  g->BayesDivide(h_jetPtwtCSV, h_jetPtnoCSV, "B");
  g->SetMarkerStyle(8);
  g->SetMaximum(1.3);
  g->GetYaxis()->SetTitle("Efficiency");  
  g->GetXaxis()->SetTitle("p_{T subjet} [GeV]");
  g->GetXaxis()->SetLimits(0,2000);

  g->Write(Form("%s_%sflavor",channel.data(),flavor.data()));
  f.Write();

  TLegend leg(0.60, 0.70, 0.90, 0.87);

  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.04);
  leg.AddEntry(g, Form("%s, %s flavor",channel.data(),flavor.data()), "lp");
  
  TLatex lar;

  lar.SetNDC(kTRUE);
  lar.SetTextSize(0.04);
  lar.SetLineWidth(5);

  TCanvas c("c", "", 0, 0, 800, 600);

  c.cd();
  g->Draw("ap");
  leg.Draw();
  lar.DrawLatex(0.15, 0.83, "CMS");
  lar.DrawLatex(0.15, 0.79, "#it{#bf{Simulation}}");
  lar.DrawLatex(0.62, 0.83, "#bf{Z+jets samples}");
  c.Print(Form("%s_%sflavor_zjetsBtagEff.pdf", channel.data(), flavor.data()));
  
}
