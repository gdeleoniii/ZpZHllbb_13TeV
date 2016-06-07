#include <string>
#include <iostream>
#include <TPad.h>
#include <TH1D.h>
#include <TFile.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TLegend.h>
#include <TTree.h>
#include <TKey.h>
#include <TSystemDirectory.h>
#include "../setNCUStyle.h"
#include "../readHists.h"

void myPlot(TH1D* h_Zjets100, TH1D* h_Zjets200, TH1D* h_Zjets400, TH1D* h_Zjets600, 
	    TH1D* h_TT, TH1D* h_WW, TH1D* h_WZ, TH1D* h_ZZ, TH1D* h_ZH,
            TH1D* h_data0, TH1D* h_data1,
	    TH1D* h_data, TH1D* h_bkg){

  h_data->Reset();
  h_data->Add(h_data0);
  h_data->Add(h_data1);

  TH1D* h_Zjets = (TH1D*)h_Zjets100->Clone("h_Zjets");

  h_Zjets->Reset();
  h_Zjets->Add(h_Zjets100);
  h_Zjets->Add(h_Zjets200);
  h_Zjets->Add(h_Zjets400);
  h_Zjets->Add(h_Zjets600);
  h_Zjets->SetFillColor(97);
  h_Zjets->SetLineColor(kBlack);

  h_TT->SetFillColor(94);
  h_TT->SetLineColor(kBlack);

  h_WW->SetFillColor(91);
  h_WW->SetLineColor(kBlack);

  h_WZ->SetFillColor(88);
  h_WZ->SetLineColor(kBlack);

  h_ZZ->SetFillColor(85);
  h_ZZ->SetLineColor(kBlack);

  h_ZH->SetFillColor(79);
  h_ZH->SetLineColor(kBlack);

  h_bkg->Reset();
  h_bkg->Add(h_Zjets);
  h_bkg->Add(h_TT);
  h_bkg->Add(h_WW);
  h_bkg->Add(h_WZ);
  h_bkg->Add(h_ZZ);
  h_bkg->Add(h_ZH);

  THStack *h_stack = new THStack("h_stack", "");

  h_stack->Add(h_Zjets);
  h_stack->Add(h_TT);
  h_stack->Add(h_WW);
  h_stack->Add(h_WZ);
  h_stack->Add(h_ZZ);
  h_stack->Add(h_ZH);

  h_data->SetLineColor(kBlack);
  h_data->SetMarkerStyle(8);
  h_data->SetMarkerSize(1.5);
  h_data->GetYaxis()->SetTitleOffset(1.3);
  h_data->GetXaxis()->SetTitle("");
  h_data->GetXaxis()->SetLabelOffset(999);
  h_data->GetXaxis()->SetLabelSize(0);

  if( h_data->GetMaximum() < h_stack->GetMaximum() ){
  
    h_stack->Draw("histe");
    h_stack->GetHistogram()->GetYaxis()->SetTitle("Event Numbers");
    h_stack->GetHistogram()->GetYaxis()->SetTitleSize(h_data->GetYaxis()->GetTitleSize());
    h_stack->GetHistogram()->GetYaxis()->SetLabelSize(h_data->GetYaxis()->GetLabelSize());
    h_stack->GetHistogram()->GetYaxis()->SetTitleOffset(1.3);
    h_stack->GetHistogram()->GetXaxis()->SetTickLength(0);
    h_stack->GetHistogram()->GetXaxis()->SetLabelOffset(999);
    h_data->Draw("elsame");
  
  }
    
  else{

    h_data->GetYaxis()->SetTitle("Event Numbers");
    h_data->Draw("el");
    h_stack->Draw("histesame");
    h_data->Draw("elsame");

  }
  
  TLegend *leg = new TLegend(0.73, 0.60, 0.90, 0.87);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
    
  leg->AddEntry(h_Zjets, "Z+Jets", "f");
  leg->AddEntry(h_TT, "t#bar{t}", "f");
  leg->AddEntry(h_WW, "WW", "f");
  leg->AddEntry(h_WZ, "WZ", "f");
  leg->AddEntry(h_ZZ, "ZZ", "f");
  leg->AddEntry(h_ZH, "ZH", "f");
  leg->AddEntry(h_data, "Data", "lp");
  leg->Draw();

  TLatex *lar = new TLatex();

  lar->SetNDC(kTRUE);
  lar->SetTextSize(0.04);
  lar->SetLineWidth(5);
  lar->DrawLatex(0.14, 0.94, "CMS #it{#bf{2015}}");
  lar->DrawLatex(0.60, 0.94, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");

}

void myRatio(TH1D* h_data, TH1D *h_bkg){

  TH1D* h_ratio = (TH1D*)h_bkg->Clone("h_ratio");

  h_ratio->Reset();

  Int_t nbin = h_ratio->GetNbinsX();
  Float_t ratio[nbin];
  Float_t error[nbin];
  Float_t numer_nbincontent[nbin];
  Float_t denom_nbincontent[nbin];
  Float_t numer_binerror[nbin];
  Float_t denom_binerror[nbin];

  for(Int_t i = 1; i <= nbin; i++){

    numer_nbincontent[i] = h_data->GetBinContent(i);
    denom_nbincontent[i] = h_bkg ->GetBinContent(i);
    numer_binerror[i]    = h_data->GetBinError(i);
    denom_binerror[i]    = h_bkg ->GetBinError(i);

    if( denom_nbincontent[i] <= 0 || numer_nbincontent[i] <= 0 ) continue;
    if( denom_binerror[i] <= 0 || numer_binerror[i] <= 0 ) continue;

    ratio[i] = (Float_t)numer_nbincontent[i]/denom_nbincontent[i];
    error[i] = (ratio[i])*sqrt(pow(numer_binerror[i]/numer_nbincontent[i],2)+pow(denom_binerror[i]/denom_nbincontent[i],2));

    h_ratio->SetBinContent(i,ratio[i]);
    h_ratio->SetBinError(i,error[i]);

  }

  h_ratio->SetLineColor(kBlack);
  h_ratio->SetMarkerStyle(8);
  h_ratio->SetMarkerSize(1.5);
  h_ratio->SetTitle("");
  h_ratio->GetYaxis()->SetTitle("Data/MC");
  h_ratio->GetYaxis()->SetTitleOffset(0.45);
  h_ratio->GetXaxis()->SetLabelSize(0.1);
  h_ratio->GetXaxis()->SetLabelOffset(0.005);
  h_ratio->GetXaxis()->SetTitleSize(0.125);
  h_ratio->GetXaxis()->SetTitleOffset(0.8);
  h_ratio->GetYaxis()->SetLabelSize(0.1);
  h_ratio->GetYaxis()->SetTitleSize(0.1);
  h_ratio->GetYaxis()->SetNdivisions(505);
  h_ratio->GetYaxis()->SetRangeUser(0,2);
  h_ratio->Draw();

  Float_t x0 = h_bkg->GetXaxis()->GetXmin();
  Float_t x1 = h_bkg->GetXaxis()->GetXmax();
  Float_t y0 = 1.;
  Float_t y1 = 1.;

  TLine* one = new TLine(x0,y0,x1,y1);

  one->SetLineColor(2);
  one->SetLineStyle(1);
  one->SetLineWidth(2);
  one->Draw("same");

  h_ratio->Draw("same");

}

void dataMCplots(std::string channel, std::string outputFolder, std::string pdfName){

  setNCUStyle(true);
 
  readHist data1   (Form("%s/Single%s-Run2015D-v1_%s.root",              outputFolder.data(), channel.data(), pdfName.data()));
  readHist data2   (Form("%s/Single%s-Run2015D-v4_%s.root",              outputFolder.data(), channel.data(), pdfName.data()));
  readHist zjets100(Form("%s/DYJetsToLL_M-50_HT-100to200_13TeV_%s.root", outputFolder.data(), pdfName.data()));
  readHist zjets200(Form("%s/DYJetsToLL_M-50_HT-200to400_13TeV_%s.root", outputFolder.data(), pdfName.data()));
  readHist zjets400(Form("%s/DYJetsToLL_M-50_HT-400to600_13TeV_%s.root", outputFolder.data(), pdfName.data()));
  readHist zjets600(Form("%s/DYJetsToLL_M-50_HT-600toInf_13TeV_%s.root", outputFolder.data(), pdfName.data()));
  readHist tt      (Form("%s/TT_TuneCUETP8M1_13TeV_%s.root",             outputFolder.data(), pdfName.data()));
  readHist ww      (Form("%s/WW_TuneCUETP8M1_13TeV_%s.root",             outputFolder.data(), pdfName.data()));
  readHist wz      (Form("%s/WZ_TuneCUETP8M1_13TeV_%s.root",             outputFolder.data(), pdfName.data()));
  readHist zz      (Form("%s/ZZ_TuneCUETP8M1_13TeV_%s.root",             outputFolder.data(), pdfName.data()));
  readHist zh      (Form("%s/ZH_HToBB_ZToLL_M125_13TeV_%s.root",         outputFolder.data(), pdfName.data()));

  Float_t up_height     = 0.8;
  Float_t dw_correction = 1.455;
  Float_t dw_height     = (1-up_height)*dw_correction;

  TCanvas c("c","",0,0,1000,900);
  c.Divide(1,2);

  TPad* c_up = (TPad*) c.GetListOfPrimitives()->FindObject("c_1");
  TPad* c_dw = (TPad*) c.GetListOfPrimitives()->FindObject("c_2"); 

  c_up->SetPad(0,1-up_height,1,1);
  c_dw->SetPad(0,0,1,dw_height);
  c_dw->SetBottomMargin(0.25);

  // To get the name of histograms
  
  TFile *f_ = TFile::Open(Form("%s/Single%s-Run2015D-v1_%s.root", outputFolder.data(), channel.data(), pdfName.data()));
  f_->cd();
  
  TDirectory *current_sourcedir = gDirectory;
  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key;

  vector<std::string> h_name;

  while( (key = (TKey*)nextkey()) ){

    TObject *obj = key->ReadObj();

    if( obj->IsA()->InheritsFrom("TH1") ) 
      h_name.push_back(obj->GetTitle());

  }

  // Draw and output
  
  for(unsigned int i = 0; i < h_name.size()-1; ++i){

    if( h_name[i]=="eleHoverE")
      c_up->cd()->SetLogy(1);
    else
      c_up->cd()->SetLogy(0);
    
    TH1D *h_data = (TH1D*)(data1.getHist(h_name[i].data()))->Clone("h_data");
    TH1D *h_bkg  = (TH1D*)(data1.getHist(h_name[i].data()))->Clone("h_bkg");

    myPlot(zjets100.getHist(h_name[i].data()),
	   zjets200.getHist(h_name[i].data()),
	   zjets400.getHist(h_name[i].data()),
	   zjets600.getHist(h_name[i].data()),
	   tt.getHist(h_name[i].data()),
	   ww.getHist(h_name[i].data()),
	   wz.getHist(h_name[i].data()),
	   zz.getHist(h_name[i].data()),
	   zh.getHist(h_name[i].data()),
	   data1.getHist(h_name[i].data()),
	   data2.getHist(h_name[i].data()),
	   h_data, h_bkg);

    c_up->RedrawAxis();
    c_dw->cd();

    myRatio(h_data, h_bkg);

    c.Draw();
    
    if( i == 0 ) 
      c.Print(Form("%s.pdf(", pdfName.data()), "pdf");
    else if( i == h_name.size()-2 ) 
      c.Print(Form("%s.pdf)", pdfName.data()), "pdf");
    else 
      c.Print(Form("%s.pdf", pdfName.data()), "pdf");
    
  }

}
