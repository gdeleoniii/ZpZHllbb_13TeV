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
#include "setNCUStyle.h"

void myPlot(TH1D* h_DY100,
            TH1D* h_DY200,
            TH1D* h_DY400,
            TH1D* h_DY600,
            TH1D* h_TTbar,
            TH1D* h_WW,
            TH1D* h_WZ,
            TH1D* h_ZZ,
            TH1D* h_data0,
	    TH1D* h_data1,
            Double_t scaleDY100,
            Double_t scaleDY200,
            Double_t scaleDY400,
            Double_t scaleDY600,
            Double_t scaleTTbar,
            Double_t scaleWW,
            Double_t scaleWZ,
            Double_t scaleZZ){

  TH1D* h_data = (TH1D*)h_data0->Clone("h_data");

  h_data->Reset();
  h_data->Add(h_data0);
  h_data->Add(h_data1);

  TH1D* h_DY = (TH1D*)h_DY100->Clone("h_DY");

  h_DY->Reset();
  h_DY->Add(h_DY100, scaleDY100);
  h_DY->Add(h_DY200, scaleDY200);
  h_DY->Add(h_DY400, scaleDY400);
  h_DY->Add(h_DY600, scaleDY600);
  h_DY->SetFillColor(kOrange-3);
  h_DY->SetLineColor(kBlack);

  h_TTbar->Scale(scaleTTbar);
  h_TTbar->SetFillColor(kGreen);
  h_TTbar->SetLineColor(kBlack);

  h_WW->Scale(scaleWW);
  h_WW->SetFillColor(kYellow);
  h_WW->SetLineColor(kBlack);

  h_WZ->Scale(scaleWZ);
  h_WZ->SetFillColor(kCyan);
  h_WZ->SetLineColor(kBlack);

  h_ZZ->Scale(scaleZZ);
  h_ZZ->SetFillColor(kPink);
  h_ZZ->SetLineColor(kBlack);

  THStack *h_stack = new THStack("h_stack", "");

  h_stack->Add(h_DY);
  h_stack->Add(h_TTbar);
  h_stack->Add(h_WW);
  h_stack->Add(h_WZ);
  h_stack->Add(h_ZZ);

  h_data->SetLineColor(kBlack);
  h_data->GetXaxis()->SetTitle("");
  h_data->GetXaxis()->SetLabelOffset(999);
  h_data->GetXaxis()->SetLabelSize(0);

  if( h_data->GetMaximum() < h_stack->GetMaximum() ){
  
    h_stack->Draw("histe");
    h_stack->GetHistogram()->GetYaxis()->SetTitle("Event Numbers");
    h_stack->GetHistogram()->GetYaxis()->SetTitleSize(h_data->GetYaxis()->GetTitleSize());
    h_stack->GetHistogram()->GetYaxis()->SetLabelSize(h_data->GetYaxis()->GetLabelSize());
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
    
  leg->AddEntry(h_DY, "DY+Jets", "f");
  leg->AddEntry(h_TTbar, "t#bar{t}", "f");
  leg->AddEntry(h_WW, "WW", "f");
  leg->AddEntry(h_WZ, "WZ", "f");
  leg->AddEntry(h_ZZ, "ZZ", "f");
  leg->AddEntry(h_data, "Data", "lp");
  leg->Draw();

  TLatex *lar = new TLatex();

  lar->SetNDC(kTRUE);
  lar->SetTextSize(0.04);
  lar->SetLineWidth(5);
  lar->DrawLatex(0.15, 0.94, "CMS preliminary 2015");
  lar->DrawLatex(0.65, 0.94, "L = 2.08 fb^{-1} at #sqrt{s} = 13 TeV");

}

void myRatio(TH1D* h_DY100,
             TH1D* h_DY200,
             TH1D* h_DY400,
             TH1D* h_DY600,
             TH1D* h_TTbar,
             TH1D* h_WW,
             TH1D* h_WZ,
             TH1D* h_ZZ,
             TH1D* h_data0,
	     TH1D* h_data1,
             Double_t scaleDY100,
             Double_t scaleDY200,
             Double_t scaleDY400,
             Double_t scaleDY600,
             Double_t scaleTTbar,
             Double_t scaleWW,
             Double_t scaleWZ,
             Double_t scaleZZ){

  TH1D* h_data = (TH1D*)h_data0->Clone("h_data");

  h_data->Reset();
  h_data->Add(h_data0);
  h_data->Add(h_data1);

  TH1D *h_bkg = (TH1D*)h_data->Clone("h_bkg");

  h_bkg->Reset();
  h_bkg->Add(h_DY100, scaleDY100);
  h_bkg->Add(h_DY200, scaleDY200);
  h_bkg->Add(h_DY400, scaleDY400);
  h_bkg->Add(h_DY600, scaleDY600);
  h_bkg->Add(h_TTbar, scaleTTbar);
  h_bkg->Add(h_WW, scaleWW);
  h_bkg->Add(h_WZ, scaleWZ);
  h_bkg->Add(h_ZZ, scaleZZ);

  TH1D* h_ratio = (TH1D*)h_data->Clone("h_ratio");

  h_ratio->Reset();

  Int_t nbin = h_ratio->GetNbinsX();
  Double_t ratio[nbin];
  Double_t error[nbin];
  Double_t numer_nbincontent[nbin];
  Double_t denom_nbincontent[nbin];
  Double_t numer_binerror[nbin];
  Double_t denom_binerror[nbin];

  for(Int_t i=1; i<=nbin; i++){

    numer_nbincontent[i] = h_data->GetBinContent(i);
    denom_nbincontent[i] = h_bkg->GetBinContent(i);
    numer_binerror[i] = h_data->GetBinError(i);
    denom_binerror[i] = h_bkg->GetBinError(i);

    if( denom_nbincontent[i] <= 0 || numer_nbincontent[i] <= 0 ) continue;
    if( denom_binerror[i] <= 0 || numer_binerror[i] <= 0 ) continue;

    ratio[i] = (Double_t)numer_nbincontent[i]/denom_nbincontent[i];
    error[i] = (ratio[i])*sqrt(pow(numer_binerror[i]/numer_nbincontent[i],2)+pow(denom_binerror[i]/denom_nbincontent[i],2));

    h_ratio->SetBinContent(i,ratio[i]);
    h_ratio->SetBinError(i,error[i]);

  }

  h_ratio->SetLineColor(kBlack);
  h_ratio->SetTitle("");
  h_ratio->GetYaxis()->SetTitle("data/MC");
  h_ratio->GetYaxis()->SetTitleOffset(0.3);
  h_ratio->GetXaxis()->SetLabelSize(0.1);
  h_ratio->GetXaxis()->SetTitleSize(0.125);
  h_ratio->GetYaxis()->SetLabelSize(0.1);
  h_ratio->GetYaxis()->SetTitleSize(0.1);
  h_ratio->GetYaxis()->SetNdivisions(505);
  h_ratio->GetYaxis()->SetRangeUser(0,2);
  h_ratio->Draw();

  Double_t x0 = h_bkg->GetXaxis()->GetXmin();
  Double_t x1 = h_bkg->GetXaxis()->GetXmax();
  Double_t y0 = 1.;
  Double_t y1 = 1.;

  TLine* one = new TLine(x0,y0,x1,y1);
  one->SetLineColor(2);
  one->SetLineStyle(1);
  one->SetLineWidth(2);
  one->Draw("same");

  h_ratio->Draw("same");

}

TFile* getFile(std::string infiles, Double_t crossSection, Double_t* scale){

  TFile* f = TFile::Open(infiles.data());
  TH1D* h = (TH1D*)(f->Get("eventWeight"));
  Double_t dataLumi = 2080; //pb-1 
  *scale = dataLumi/(h->Integral()/crossSection);

  return f;

}

void dataMCplots(std::string outputFolder, std::string pdfName){

  std::vector<string> infiles;
 
  TSystemDirectory *base = new TSystemDirectory("root","root");
  base->SetDirectory(outputFolder.data());
  TList *listOfFiles = base->GetListOfFiles();
  TIter fileIt(listOfFiles);
  TFile *fileH = new TFile();
  Long64_t nfiles = 0;

  while( (fileH = (TFile*)fileIt()) ){
    
    std::string fileN = fileH->GetName();
    std::string baseString = "root";
    if( fileN.find(baseString) == std::string::npos ) continue;
    infiles.push_back(Form("%s/%s",outputFolder.data(),fileN.data()));
    nfiles++;
    
  }

  TFile *f_DY100 = NULL;
  TFile *f_DY200 = NULL;
  TFile *f_DY400 = NULL;
  TFile *f_DY600 = NULL;
  TFile *f_TTbar = NULL;
  TFile *f_WW    = NULL;
  TFile *f_WZ    = NULL;
  TFile *f_ZZ    = NULL;
  TFile *f_data0 = NULL;
  TFile *f_data1 = NULL;

  Double_t xSecDY100 = 139.4*1.23;
  Double_t xSecDY200 = 42.75*1.23;
  Double_t xSecDY400 = 5.497*1.23;
  Double_t xSecDY600 = 2.21*1.23;
  Double_t xSecTTbar = 831.76;
  Double_t xSecWW    = 118.7;
  Double_t xSecWZ    = 47.13;
  Double_t xSecZZ    = 16.523;

  Double_t scaleDY100 = 0;
  Double_t scaleDY200 = 0;
  Double_t scaleDY400 = 0;
  Double_t scaleDY600 = 0;
  Double_t scaleTTbar = 0;
  Double_t scaleWW    = 0;
  Double_t scaleWZ    = 0;
  Double_t scaleZZ    = 0;

  Double_t dummy = -1;

  for(unsigned int i = 0; i < infiles.size(); i++){

    cout << "Input file: " << infiles[i] << endl;

    if( infiles[i].find("HT-100") != std::string::npos )    
      f_DY100 = getFile(infiles[i].data(), xSecDY100, &scaleDY100);

    if( infiles[i].find("HT-200") != std::string::npos )
      f_DY200 = getFile(infiles[i].data(), xSecDY200, &scaleDY200);

    if( infiles[i].find("HT-400") != std::string::npos )
      f_DY400 = getFile(infiles[i].data(), xSecDY400, &scaleDY400); 

    if( infiles[i].find("HT-600") != std::string::npos )
      f_DY600 = getFile(infiles[i].data(), xSecDY600, &scaleDY600);

    if( infiles[i].find("TT_") != std::string::npos )
      f_TTbar = getFile(infiles[i].data(), xSecTTbar, &scaleTTbar);

    if( infiles[i].find("WW_") != std::string::npos )
      f_WW = getFile(infiles[i].data(), xSecWW, &scaleWW);

    if( infiles[i].find("WZ_") != std::string::npos )
      f_WZ = getFile(infiles[i].data(), xSecWZ, &scaleWZ);

    if( infiles[i].find("ZZ_") != std::string::npos )
      f_ZZ = getFile(infiles[i].data(), xSecZZ, &scaleZZ);

    if( infiles[i].find("V12015") != std::string::npos )
      f_data0 = getFile(infiles[i].data(), dummy, &dummy);

    if( infiles[i].find("V42015") != std::string::npos )
      f_data1 = getFile(infiles[i].data(), dummy, &dummy);

  }

  setNCUStyle(true);
  
  Double_t up_height     = 0.8;
  Double_t dw_correction = 1.455;
  Double_t dw_height     = (1-up_height)*dw_correction;

  TCanvas c("c","",0,0,1000,800);
  c.Divide(1,2);

  TPad* c_up = (TPad*) c.GetListOfPrimitives()->FindObject("c_1");
  TPad* c_dw = (TPad*) c.GetListOfPrimitives()->FindObject("c_2"); 

  c_up->SetPad(0,1-up_height,1,1);
  c_dw->SetPad(0,0,1,dw_height);
  c_dw->SetBottomMargin(0.25);

  // To get the name of histograms
  
  TFile *f_ = TFile::Open(infiles[0].data());
  f_->cd();
  
  TDirectory *current_sourcedir = gDirectory;
  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key;

  vector<std::string> h_name;

  while ( (key = (TKey*)nextkey()) ) {

    TObject *obj = key->ReadObj();

    if ( obj->IsA()->InheritsFrom("TH1") ) 
      h_name.push_back(obj->GetTitle());

  }

  // Draw and output
  
  for(unsigned int i = 0; i < h_name.size()-1; i++){

    if( h_name[i]=="eleHoverE" || h_name[i]=="eleMiniIso" || h_name[i]=="muMiniIso" )
      c_up->cd()->SetLogy(1);
    else
      c_up->cd()->SetLogy(0);
    
    myPlot(((TH1D*)(f_DY100->Get(h_name[i].data()))),
	   ((TH1D*)(f_DY200->Get(h_name[i].data()))),
	   ((TH1D*)(f_DY400->Get(h_name[i].data()))),
	   ((TH1D*)(f_DY600->Get(h_name[i].data()))),
	   ((TH1D*)(f_TTbar->Get(h_name[i].data()))),
	   ((TH1D*)(f_WW->Get(h_name[i].data()))),
	   ((TH1D*)(f_WZ->Get(h_name[i].data()))),
	   ((TH1D*)(f_ZZ->Get(h_name[i].data()))),
	   ((TH1D*)(f_data0->Get(h_name[i].data()))),
	   ((TH1D*)(f_data1->Get(h_name[i].data()))),
	   scaleDY100,
	   scaleDY200,
	   scaleDY400,
	   scaleDY600,
	   scaleTTbar,
	   scaleWW,
	   scaleWZ,
	   scaleZZ);

    c_up->RedrawAxis();

    c_dw->cd();
    myRatio(((TH1D*)(f_DY100->Get(h_name[i].data()))),
	    ((TH1D*)(f_DY200->Get(h_name[i].data()))),
	    ((TH1D*)(f_DY400->Get(h_name[i].data()))),
	    ((TH1D*)(f_DY600->Get(h_name[i].data()))),
	    ((TH1D*)(f_TTbar->Get(h_name[i].data()))),
	    ((TH1D*)(f_WW->Get(h_name[i].data()))),
	    ((TH1D*)(f_WZ->Get(h_name[i].data()))),
	    ((TH1D*)(f_ZZ->Get(h_name[i].data()))),
	    ((TH1D*)(f_data0->Get(h_name[i].data()))),
	    ((TH1D*)(f_data1->Get(h_name[i].data()))),
	    scaleDY100,
	    scaleDY200,
	    scaleDY400,
	    scaleDY600,
	    scaleTTbar,
	    scaleWW,
	    scaleWZ,
	    scaleZZ);

    c.Draw();
    
    if( i == 0 ) c.Print(Form("%s.pdf(", pdfName.data()), "pdf");
    else if( i == h_name.size()-2 ) c.Print(Form("%s.pdf)", pdfName.data()), "pdf");
    else c.Print(Form("%s.pdf", pdfName.data()), "pdf");
    
  }

}
