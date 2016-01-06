#include <string>
#include <vector>
#include <iostream>
#include <TF1.h>
#include <TH1.h>
#include <TPad.h>
#include <TMath.h>
#include <TFile.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMatrixD.h>
#include <TFitResult.h>
#include <TSystemDirectory.h>
#include <TGraphAsymmErrors.h>
#include "../setNCUStyle.h"

const double xmin  = 500;
const double xmax  = 5000;
const int    nBins = (xmax-xmin)/100;

double dataLumi  = 3000; //pb-1
double xSecDY100 = 147.4;
double xSecDY200 = 40.99;
double xSecDY400 = 5.678;
double xSecDY600 = 2.198;

TFile* getFile(std::string infiles, std::string hname, 
	       double crossSection, double* scale){

  TFile* f = TFile::Open(infiles.data());
  TH1D*  h = NULL;

  if( hname.find("pMC") != std::string::npos ) 
    h = (TH1D*)(f->Get("eventWeight_pMC"));
  else if( hname.find("pDA") != std::string::npos )
    h = (TH1D*)(f->Get("eventWeight_pDA"));

  *scale = dataLumi/(h->Integral()/crossSection);

  return f;

}

TH1D* addSamples(std::vector<string>& infiles, std::string hname,
		 TFile* f_DY100, TFile* f_DY200, TFile* f_DY400, TFile* f_DY600){ 

  double scaleDY100 = 0;
  double scaleDY200 = 0;
  double scaleDY400 = 0;
  double scaleDY600 = 0;

  for(unsigned int i = 0; i < infiles.size(); i++){

    if( infiles[i].find("HT-100") != std::string::npos )
      f_DY100 = getFile(infiles[i].data(), hname.data(), xSecDY100, &scaleDY100);

    if( infiles[i].find("HT-200") != std::string::npos )
      f_DY200 = getFile(infiles[i].data(), hname.data(), xSecDY200, &scaleDY200);

    if( infiles[i].find("HT-400") != std::string::npos )
      f_DY400 = getFile(infiles[i].data(), hname.data(), xSecDY400, &scaleDY400);

    if( infiles[i].find("HT-600") != std::string::npos )
      f_DY600 = getFile(infiles[i].data(), hname.data(), xSecDY600, &scaleDY600);

  }

  TH1D* DY100_temp = (TH1D*)(f_DY100->Get(Form("%s",hname.c_str())));
  TH1D* DY200_temp = (TH1D*)(f_DY200->Get(Form("%s",hname.c_str())));
  TH1D* DY400_temp = (TH1D*)(f_DY400->Get(Form("%s",hname.c_str())));
  TH1D* DY600_temp = (TH1D*)(f_DY600->Get(Form("%s",hname.c_str())));

  TH1D* h_Total = (TH1D*)(f_DY100->Get(Form("%s",hname.c_str())))->Clone("h_Total");

  h_Total->Reset();
  h_Total->Add(DY100_temp,scaleDY100);
  h_Total->Add(DY200_temp,scaleDY200);
  h_Total->Add(DY400_temp,scaleDY400);
  h_Total->Add(DY600_temp,scaleDY600);

  return h_Total;

}

void myRatio(TH1D* h_numer, TH1D* h_denom){

  TH1D* h_ratio = (TH1D*)h_numer->Clone("h_ratio");

  h_ratio->Reset();

  int nbin = h_ratio->GetNbinsX();
  double ratio[nbin];
  double error[nbin];
  double numer_nbincontent[nbin];
  double denom_nbincontent[nbin];
  double numer_binerror[nbin];
  double denom_binerror[nbin];

  for( int i = 1; i <= nbin; i++ ){

    numer_nbincontent[i] = h_numer->GetBinContent(i);
    denom_nbincontent[i] = h_denom->GetBinContent(i);
    numer_binerror[i]    = h_numer->GetBinError(i);
    denom_binerror[i]    = h_denom->GetBinError(i);

    if( denom_nbincontent[i] <= 0 || numer_nbincontent[i] <= 0 ) continue;
    if( denom_binerror[i] <= 0 || numer_binerror[i] <= 0 ) continue;

    ratio[i] = (double)numer_nbincontent[i]/denom_nbincontent[i];
    error[i] = (ratio[i])*sqrt(pow(numer_binerror[i]/numer_nbincontent[i],2)+pow(denom_binerror[i]/denom_nbincontent[i],2));

    h_ratio->SetBinContent(i,ratio[i]);
    h_ratio->SetBinError(i,error[i]);

  }

  h_ratio->SetLineColor(kBlack);
  h_ratio->SetTitle("");
  h_ratio->GetYaxis()->SetTitle("Predicted/Truth");
  h_ratio->GetYaxis()->SetTitleOffset(0.3);
  h_ratio->GetXaxis()->SetTitle("ZH mass in signal region of pseudo-data");
  h_ratio->GetXaxis()->SetLabelSize(0.1);
  h_ratio->GetXaxis()->SetTitleSize(0.125);
  h_ratio->GetYaxis()->SetLabelSize(0.1);
  h_ratio->GetYaxis()->SetTitleSize(0.1);
  h_ratio->GetYaxis()->SetNdivisions(505);
  h_ratio->GetYaxis()->SetRangeUser(0,2);
  h_ratio->Draw();

  double x0 = h_denom->GetXaxis()->GetXmin();
  double x1 = h_denom->GetXaxis()->GetXmax();
  double y0 = 1.;
  double y1 = 1.;

  TLine* one = new TLine(x0,y0,x1,y1);

  one->SetLineColor(2);
  one->SetLineStyle(1);
  one->SetLineWidth(2);
  one->Draw("same");

  h_ratio->Draw("same");

}

double fitZpmass(double* v, double* p){

  double x = v[0];
  return p[0]*TMath::Exp(p[1]*x + p[2]/x);

}

double divFunc(double* v, double* p){

  double x = v[0];
  return (p[0]*TMath::Exp(p[1]*x+p[2]/x))/(p[3]*TMath::Exp(p[4]*x+p[5]/x));

}

double ErfExp(double x, double c, double offset, double width){
  
  if( width < 1e-2 ) width = 1e-2;
  if( c == 0 ) c = -1e-7;
  return TMath::Exp(c*x)*(1.+TMath::Erf((x-offset)/width))/2.;

}

double integral_ErfExp(const double c, const double offset, double width, const double xmin, const double xmax){

  double width_tmp = width; 
  double minTerm   = 0;
  double maxTerm   = 0;

  if( width < 1e-2 ) width = 1e-2;

  if( c == 0 ){ 
    
    double delta = -1e-7;

    minTerm = (TMath::Exp(delta*delta*width_tmp*width_tmp/4+delta*offset) * 
	       TMath::Erf((2*xmin-delta*width_tmp*width_tmp-2*offset)/2/width_tmp) - 
	       TMath::Exp(delta*xmin) * 
	       TMath::Erf((xmin-offset)/width_tmp) - 
	       TMath::Exp(delta*xmin))/-2/delta;

    maxTerm = (TMath::Exp(delta*delta*width_tmp*width_tmp/4+delta*offset) * 
	       TMath::Erf((2*xmax-delta*width_tmp*width_tmp-2*offset)/2/width_tmp) - 
	       TMath::Exp(delta*xmax) * 
	       TMath::Erf((xmax-offset)/width_tmp) - 
	       TMath::Exp(delta*xmax))/-2/delta;
        
  }

  else{

    minTerm = (TMath::Exp(c*c*width_tmp*width_tmp/4+c*offset) * 
	       TMath::Erf((2*xmin-c*width_tmp*width_tmp-2*offset)/2/width_tmp) - 
	       TMath::Exp(c*xmin) * 
	       TMath::Erf((xmin-offset)/width_tmp) - 
	       TMath::Exp(c*xmin))/-2/c;

    maxTerm = (TMath::Exp(c*c*width_tmp*width_tmp/4+c*offset) * 
	       TMath::Erf((2*xmax-c*width_tmp*width_tmp-2*offset)/2/width_tmp) - 
	       TMath::Exp(c*xmax) * 
	       TMath::Erf((xmax-offset)/width_tmp) - 
	       TMath::Exp(c*xmax))/-2/c;
  }

  return (maxTerm-minTerm);

}  

double fitPRmass(double* v, double* p){

  double x = v[0];
  Double_t width_tmp = p[3];
  Double_t binwidth  = p[4];

  if( p[3] < 1e-2 ) width_tmp = 1e-2;

  return p[0]*ErfExp(x,p[1],p[2],width_tmp)/integral_ErfExp(p[1],p[2],width_tmp,40,240)*binwidth;

}

double hollow_fitPRmass(double* v, double* p){

  double x = v[0];
  double width_tmp = p[3];
  double binwidth  = p[4];

  if( p[3] < 1e-2 ) width_tmp = 1e-2;

  return p[0]*ErfExp(x,p[1],p[2],width_tmp)/(integral_ErfExp(p[1],p[2],width_tmp,40,65)+integral_ErfExp(p[1],p[2],width_tmp,145,240))*binwidth;

}

TGraphAsymmErrors* fitUncertainty(bool considerBins, const TF1* f, TMatrixD* corrMatrix, double (*fitFunc)(double*, double*),
				  TH1D* h, double nBkgSig, double* posUncEv, double* negUncEv){

  double par[5] = {0};

  for( int i = 0; i < 5; i++ )
    par[i] = f->GetParameter(i);

  TF1* posFit[5];
  TF1* negFit[5];

  for( int i = 0; i < 5; i++ ){

    double partemp[5] = {par[0],par[1],par[2],par[3],par[4]};

    posFit[i] = new TF1(Form("posFit%d",i), fitFunc, 40, 240, 5);
    partemp[i]  = par[i] + f->GetParError(i);
    posFit[i]->SetParameters(partemp[0],partemp[1],partemp[2],partemp[3],partemp[4]);

  }

  for( int i = 0; i < 5; i++ ){

    double partemp[5] = {par[0],par[1],par[2],par[3],par[4]};

    negFit[i] = new TF1(Form("negFit%d",i), fitFunc, 40, 240, 5);
    partemp[i]  = par[i] - f->GetParError(i);
    negFit[i]->SetParameters(partemp[0],partemp[1],partemp[2],partemp[3],partemp[4]);

  }

  TMatrixD posColM(5,1);
  TMatrixD negColM(5,1);
  TMatrixD posRowM(1,5);
  TMatrixD negRowM(1,5);

  TGraphAsymmErrors* g = NULL;

  /// for error band of fit curve ///

  if( considerBins == true ){

    int    NBINS = 40;
    double x     = 40.0;
    double width = (240-x)/NBINS;

    double funcX[NBINS];
    double funcY[NBINS];
    double posUnc[NBINS];
    double negUnc[NBINS];

    for( int n = 0; n < NBINS; n++){

      for(int i = 0; i < 5; i++){
    
	posColM(i,0) = fabs(f->Eval(x) - posFit[i]->Eval(x));
	negColM(i,0) = fabs(f->Eval(x) - negFit[i]->Eval(x));
	posRowM(0,i) = posColM(i,0);
	negRowM(0,i) = negColM(i,0);
    
      }

      TMatrixD posTemp = posRowM*(*corrMatrix*posColM);
      TMatrixD negTemp = negRowM*(*corrMatrix*negColM);
    
      posUnc[n] = TMath::Sqrt(posTemp(0,0));
      negUnc[n] = TMath::Sqrt(negTemp(0,0));

      funcX[n] = x;
      funcY[n] = f->Eval(x);

      x += width;

    }

    g = new TGraphAsymmErrors(NBINS, funcX, funcY, 0, 0, negUnc, posUnc);

  }

  /// for uncertainties of event number ///

  else if( considerBins == false ){

    for(int i = 0; i < 5; i++){

      posColM(i,0) = fabs(nBkgSig - posFit[i]->Integral(105,135)/h->GetBinWidth(1));
      negColM(i,0) = fabs(nBkgSig - negFit[i]->Integral(105,135)/h->GetBinWidth(1));
      posRowM(0,i) = posColM(i,0);
      negRowM(0,i) = negColM(i,0);

    }

    TMatrixD posTemp = posRowM*(*corrMatrix*posColM);
    TMatrixD negTemp = negRowM*(*corrMatrix*negColM);

    *posUncEv = TMath::Sqrt(posTemp(0,0));
    *negUncEv = TMath::Sqrt(negTemp(0,0));

  }

  return g;

}

/// Main function start ///

void anotherfit(std::string outputFolder){

  setNCUStyle();
  gStyle->SetOptFit(0);
  gStyle->SetTitleSize(0.04,"XYZ");
  gStyle->SetLabelSize(0.03,"XYZ");

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

  // Declare prefer histogram and add them together

  TH1D* h_hollow_PRmassBKG = addSamples(infiles,"corrPRmass_pMC",f_DY100,f_DY200,f_DY400,f_DY600);
  TH1D* h_PRmassBKG     = addSamples(infiles,"corrPRmassAll_pMC",f_DY100,f_DY200,f_DY400,f_DY600);
  TH1D* h_hollow_PRmass = addSamples(infiles,"corrPRmass_pDA",f_DY100,f_DY200,f_DY400,f_DY600);
  TH1D* h_PRmass        = addSamples(infiles,"corrPRmassAll_pDA",f_DY100,f_DY200,f_DY400,f_DY600);

  h_hollow_PRmassBKG->SetMarkerStyle(8);
  h_hollow_PRmassBKG->SetMarkerSize(1.5);
  h_hollow_PRmassBKG->SetLineColor(kBlack);
  h_hollow_PRmassBKG->SetXTitle("Side band corrected pruned mass in pseudo-data");
  h_hollow_PRmassBKG->SetYTitle("Event numbers");
  h_hollow_PRmassBKG->SetTitleFont(62);

  h_PRmassBKG->SetMarkerStyle(8);
  h_PRmassBKG->SetMarkerSize(1.5);
  h_PRmassBKG->SetLineColor(kBlack);
  h_PRmassBKG->SetXTitle("Side band corrected pruned mass in MC");
  h_PRmassBKG->SetYTitle("Event numbers");
  h_PRmassBKG->SetTitleFont(62);

  h_hollow_PRmass->SetMarkerStyle(8);
  h_hollow_PRmass->SetMarkerSize(1.5);
  h_hollow_PRmass->SetLineColor(kBlack);
  h_hollow_PRmass->SetXTitle("Side band corrected pruned mass in pseudo-data");
  h_hollow_PRmass->SetYTitle("Event numbers");
  h_hollow_PRmass->SetTitleFont(62);

  h_PRmass->SetMarkerStyle(8);
  h_PRmass->SetMarkerSize(1.5);
  h_PRmass->SetLineColor(kBlack);
  h_PRmass->SetXTitle("Corrected pruned mass in pseudo-data");
  h_PRmass->SetYTitle("Event numbers");
  h_PRmass->SetTitleFont(62);

  // Make the statistics error more like data

  for( int i = 1; i <= h_PRmass->GetNbinsX(); i++ ){

    h_hollow_PRmassBKG->SetBinError(i,TMath::Sqrt(h_hollow_PRmassBKG->GetBinContent(i)));
    h_PRmassBKG->SetBinError(i,TMath::Sqrt(h_PRmassBKG->GetBinContent(i)));
    h_hollow_PRmass->SetBinError(i,TMath::Sqrt(h_hollow_PRmass->GetBinContent(i)));
    h_PRmass->SetBinError(i,TMath::Sqrt(h_PRmass->GetBinContent(i)));

  }

  // Fit pruned mass with signal region in MC (BKG)

  TF1* f_fitPRmassBKG = new TF1("f_fitPRmassBKG", fitPRmass, 40, 240, 5);

  f_fitPRmassBKG->SetLineWidth(2);
  f_fitPRmassBKG->SetLineColor(kBlue);

  double parFitPRm[4] = {1224,-0.107,139.6,107.4};

  f_fitPRmassBKG->SetParameters(parFitPRm[0],parFitPRm[1],parFitPRm[2],parFitPRm[3],h_PRmassBKG->GetBinWidth(1));
  f_fitPRmassBKG->FixParameter(4,h_PRmassBKG->GetBinWidth(1));
  f_fitPRmassBKG->FixParameter(0,h_PRmassBKG->Integral());

  h_PRmassBKG->Fit("f_fitPRmassBKG", "Q", "", 40, 240);
  h_PRmassBKG->SetMaximum(190);

  // Fit pruned mass with signal region  

  TF1* f_fitPRmass = new TF1("f_fitPRmass", fitPRmass, 40, 240, 5);

  f_fitPRmass->SetLineWidth(2);
  f_fitPRmass->SetLineColor(kBlue);

  f_fitPRmass->FixParameter(1,f_fitPRmassBKG->GetParameter(1));
  f_fitPRmass->FixParameter(2,f_fitPRmassBKG->GetParameter(2));
  f_fitPRmass->FixParameter(3,f_fitPRmassBKG->GetParameter(3));
  f_fitPRmass->FixParameter(4,h_PRmass->GetBinWidth(1));

  h_PRmass->Fit("f_fitPRmass", "Q", "", 40, 240);
  h_PRmass->SetMaximum(190);

  double chisqr_cpma = f_fitPRmass->GetChisquare();
  int ndf_cpma = f_fitPRmass->GetNDF();

  TFitResultPtr fitptr = h_PRmass->Fit(f_fitPRmass, "QS");
  TFitResult fitresult = (*fitptr);
  TMatrixD corrMatrix  = fitresult.GetCorrelationMatrix();

  double dummy = 0.0;

  TGraphAsymmErrors* g_errorBands = fitUncertainty(true, f_fitPRmass, &corrMatrix, fitPRmass, NULL, 0, &dummy, &dummy);

  g_errorBands->SetFillStyle(1001);
  g_errorBands->SetFillColor(kYellow);

  // Fit pruned mass without signal region in MC (BKG)

  TF1* f_hollow_fitPRmassBKG = new TF1("f_hollow_fitPRmassBKG", hollow_fitPRmass, 40, 240, 5);

  f_hollow_fitPRmassBKG->SetLineWidth(2);
  f_hollow_fitPRmassBKG->SetLineColor(kBlue);
  f_hollow_fitPRmassBKG->SetParameters(parFitPRm[0],parFitPRm[1],parFitPRm[2],parFitPRm[3],h_hollow_PRmassBKG->GetBinWidth(1));
  f_hollow_fitPRmassBKG->FixParameter(4,h_hollow_PRmassBKG->GetBinWidth(1));
  f_hollow_fitPRmassBKG->FixParameter(0,h_hollow_PRmassBKG->Integral());

  h_hollow_PRmassBKG->Fit("f_hollow_fitPRmassBKG", "Q", "", 40, 240);
  h_hollow_PRmassBKG->SetMaximum(190);

  // Fit pruned mass without signal region 

  TF1* f_hollow_fitPRmass = new TF1("f_hollow_fitPRmass", hollow_fitPRmass, 40, 240, 5);

  f_hollow_fitPRmass->SetLineWidth(2);
  f_hollow_fitPRmass->SetLineColor(kBlue);
 
  f_hollow_fitPRmass->FixParameter(1,f_hollow_fitPRmassBKG->GetParameter(1));
  f_hollow_fitPRmass->FixParameter(2,f_hollow_fitPRmassBKG->GetParameter(2));
  f_hollow_fitPRmass->FixParameter(3,f_hollow_fitPRmassBKG->GetParameter(3));
  f_hollow_fitPRmass->FixParameter(4,h_hollow_PRmass->GetBinWidth(1));
 
  h_hollow_PRmass->Fit("f_hollow_fitPRmass", "Q", "", 40, 240);
  h_hollow_PRmass->SetMaximum(190);

  double chisqr_cpm = f_hollow_fitPRmass->GetChisquare();
  int ndf_cpm = f_hollow_fitPRmass->GetNDF();

  TFitResultPtr hollow_fitptr = h_hollow_PRmass->Fit(f_hollow_fitPRmass, "QS");
  TFitResult hollow_fitresult = (*hollow_fitptr);
  TMatrixD hollow_corrMatrix  = hollow_fitresult.GetCorrelationMatrix();

  TGraphAsymmErrors* g_hollow_errorBands = fitUncertainty(true, f_hollow_fitPRmass, &hollow_corrMatrix, hollow_fitPRmass, NULL, 0, &dummy, &dummy);

  g_hollow_errorBands->SetFillStyle(1001);
  g_hollow_errorBands->SetFillColor(kYellow);


  //// Fluctuate the pruned mass histogram to test the fitting results ////       

  TH1D* h_fluc = (TH1D*)h_PRmass->Clone("h_fluc");
  TH1D* h_hollow_fluc = (TH1D*)h_PRmass->Clone("h_hollow_fluc");
  TH1D* h_bias = new TH1D("h_bias", "", 100, -1, 1);

  TCanvas* ctemp = new TCanvas("ctemp","",0,0,1000,800);
  TFile* outemp = new TFile("anotherfit.root","recreate");

  h_bias->SetLineWidth(1);
  h_bias->SetFillColor(kYellow);
  h_bias->GetXaxis()->SetTitle("Bias ((fit-true)/true)");
  h_bias->GetYaxis()->SetTitle("Counts");

  for( int ntoy = 0; ntoy < 1000; ntoy++ ){
    if(ntoy%10==1)cout<<ntoy<<endl;
    h_fluc->Reset();
    h_hollow_fluc->Reset();

    for( int idata = 0; idata < 1e5/*(int)(h_PRmass->Integral())*/; idata++ ){

      h_fluc->Fill(f_fitPRmass->GetRandom(40,240));
     
    }

    int nd1 = (105 - h_fluc->GetBinLowEdge(1))/h_fluc->GetBinWidth(1);
    int nd2 = (135 - h_fluc->GetBinLowEdge(1))/h_fluc->GetBinWidth(1);

    double nSigHist = h_fluc->Integral(nd1,nd2);

    for( int nbin = 1; nbin <= h_fluc->GetNbinsX(); nbin++ ){

      if( h_hollow_fluc->GetBinLowEdge(nbin) < 65 || h_hollow_fluc->GetBinLowEdge(nbin) > 140 ){

	h_hollow_fluc->SetBinContent(nbin,h_fluc->GetBinContent(nbin));
	h_hollow_fluc->SetBinError(nbin,h_fluc->GetBinError(nbin));

      }

    }

    TF1* f_fluc = new TF1("f_fluc", hollow_fitPRmass, 40, 240, 5);

    f_fluc->FixParameter(1,f_hollow_fitPRmassBKG->GetParameter(1));
    f_fluc->FixParameter(2,f_hollow_fitPRmassBKG->GetParameter(2));
    f_fluc->FixParameter(3,f_hollow_fitPRmassBKG->GetParameter(3));
    f_fluc->FixParameter(4,h_hollow_fluc->GetBinWidth(1));


    h_hollow_fluc->Fit("f_fluc", "Q", "", 40, 240);

    // temp region
    ctemp->cd();
    h_fluc->SetLineColor(kRed);
    h_fluc->Draw();
    h_hollow_fluc->Draw("same");
    ctemp->Write();

    double nSigFit = f_fluc->Integral(105,135)/h_hollow_fluc->GetBinWidth(1);

    h_bias->Fill((nSigFit - nSigHist)/nSigHist);

  }
  outemp->Write();
  //// End of the test ////  

  // Output results

  TLegend* leg = new TLegend(0.35, 0.77, 0.87, 0.87);

  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);

  TLatex* lar = new TLatex();

  lar->SetTextSize(0.035);
  lar->SetLineWidth(5);

  TCanvas* c = new TCanvas("c","",0,0,1000,800);

  c->cd()->SetLogy(0);
  h_PRmass->Draw();
  g_errorBands->Draw("3same");
  h_PRmass->Draw("same");
  leg->Clear();
  leg->AddEntry(h_PRmass, "Error = #sqrt{N_{per bin}}", "lp");
  leg->AddEntry(g_errorBands, "Uncertainty based on fitting errors", "f");
  leg->Draw();
  lar->DrawLatexNDC(0.50, 0.65, Form("#chi^{2} / ndf: %f / %d",chisqr_cpma,ndf_cpma));
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2015");
  lar->DrawLatexNDC(0.65, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  c->Print("anotherfit.pdf(");
  
  c->cd()->SetLogy(0);
  h_hollow_PRmass->Draw();
  g_hollow_errorBands->Draw("3same");
  h_hollow_PRmass->Draw("same");
  leg->Clear();
  leg->AddEntry(h_hollow_PRmass, "Error = #sqrt{N_{per bin}}", "lp");
  leg->AddEntry(g_errorBands, "Uncertainty based on fitting errors", "f");
  leg->Draw();
  lar->DrawLatexNDC(0.50, 0.65, Form("#chi^{2} / ndf: %f / %d",chisqr_cpm,ndf_cpm));
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2015");
  lar->DrawLatexNDC(0.65, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  c->Print("anotherfit.pdf");
  
  c->cd()->SetLogy();
  h_hollow_PRmass->SetMaximum(3e3);
  h_hollow_PRmass->Draw();
  g_hollow_errorBands->Draw("3same");
  h_hollow_PRmass->Draw("same");
  leg->Clear();
  leg->AddEntry(h_hollow_PRmass, "Error = #sqrt{N_{per bin}}", "lp");
  leg->AddEntry(g_hollow_errorBands, "Uncertainty based on fitting errors", "f");
  leg->Draw();
  lar->DrawLatexNDC(0.25, 0.40, Form("#chi^{2} / ndf: %f / %d",chisqr_cpm,ndf_cpm));
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2015");
  lar->DrawLatexNDC(0.65, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  c->Print("anotherfit.pdf");

  c->cd()->SetLogy(0);
  h_bias->Draw();
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2015");
  lar->DrawLatexNDC(0.65, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  c->Print("anotherfit.pdf)");
  
}
