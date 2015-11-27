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
double xSecDY100 = 139.4*1.23;
double xSecDY200 = 42.75*1.23;
double xSecDY400 = 5.497*1.23;
double xSecDY600 = 2.21*1.23;

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
    numer_binerror[i] = h_numer->GetBinError(i);
    denom_binerror[i] = h_denom->GetBinError(i);

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

double fitPRmass(double* v, double* p){

  double x = v[0];
  return p[0]*TMath::Exp(p[1]*x)*0.5*(1+TMath::Erf((x-p[2])/p[3]));

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

double eiko_fitPRmass(double* v, double* p){

  double x = v[0];
  double width_tmp = p[3];
  double binwidth  = p[4];

  if( p[3] < 1e-2 ) width_tmp = 1e-2;

  return p[0]*ErfExp(x,p[1],p[2],width_tmp)/integral_ErfExp(p[1],p[2],width_tmp, 40, 240)*binwidth ; 

}

double hollow_fitPRmass(double* v, double* p){

  double x = v[0];
  double width_tmp = p[3];
  double binwidth  = p[4];

  if( p[3] < 1e-2 ) width_tmp = 1e-2;

  return p[0]*ErfExp(x,p[1],p[2],width_tmp)/(integral_ErfExp(p[1],p[2],width_tmp,40,65)+integral_ErfExp(p[1],p[2],width_tmp,145,240))*binwidth ; 

}

TGraphAsymmErrors* fitUncertainty(const TF1* f, const TMatrixD* corrMatrix){

  double par[4] = {0};

  for( int i = 0; i < 4; i++ )
    par[i] = f->GetParameter(i);

  TF1* posFit[4];
  TF1* negFit[4];

  for( int i = 0; i < 4; i++ ){

    double partemp[4] = {par[0],par[1],par[2],par[3]};

    posFit[i]  = new TF1(Form("posFit%d",i), fitPRmass, 40, 240, 4);
    partemp[i] = par[i] + f->GetParError(i);
    posFit[i]->SetParameters(partemp[0],partemp[1],partemp[2],partemp[3]);

  }

  for( int i = 0; i < 4; i++ ){

    double partemp[4] = {par[0],par[1],par[2],par[3]};

    negFit[i]  = new TF1(Form("negFit%d",i), fitPRmass, 40, 240, 4);
    partemp[i] = par[i] - f->GetParError(i);
    negFit[i]->SetParameters(partemp[0],partemp[1],partemp[2],partemp[3]);

  }

  TMatrixD posColM(4,1);
  TMatrixD negColM(4,1);
  TMatrixD posRowM(1,4);
  TMatrixD negRowM(1,4);

  int    NBINS = 40;
  double x     = 40.0;
  double width = (240-x)/NBINS;

  double funcX[NBINS];
  double funcY[NBINS];
  double posUnc[NBINS];
  double negUnc[NBINS];

  for( int n = 0; n < NBINS; n++){

    for(int i = 0; i < 4; i++){
    
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

  TGraphAsymmErrors* g = new TGraphAsymmErrors(NBINS, funcX, funcY, 0, 0, negUnc, posUnc);

  return g;

}

void fitUncNoBins(const TF1* f, const TMatrixD* corrMatrix, TH1D* h,
		  double nBkgSig, double* posUnc, double* negUnc){

  double par[4] = {0};

  for( int i = 0; i < 4; i++ )
    par[i] = f->GetParameter(i);

  TF1* posFit[4];
  TF1* negFit[4];

  for( int i = 0; i < 4; i++ ){

    double partemp[4] = {par[0],par[1],par[2],par[3]};

    posFit[i]  = new TF1(Form("posFit%d",i), fitPRmass, 40, 240, 4);
    partemp[i] = par[i] + f->GetParError(i);
    posFit[i]->SetParameters(partemp[0],partemp[1],partemp[2],partemp[3]);

  }

  for( int i = 0; i < 4; i++ ){

    double partemp[4] = {par[0],par[1],par[2],par[3]};

    negFit[i]  = new TF1(Form("negFit%d",i), fitPRmass, 40, 240, 4);
    partemp[i] = par[i] - f->GetParError(i);
    negFit[i]->SetParameters(partemp[0],partemp[1],partemp[2],partemp[3]);

  }

  TMatrixD posColM(4,1);
  TMatrixD negColM(4,1);
  TMatrixD posRowM(1,4);
  TMatrixD negRowM(1,4);

  for(int i = 0; i < 4; i++){
    
    posColM(i,0) = fabs(nBkgSig - posFit[i]->Integral(105,135)/h->GetBinWidth(1));
    negColM(i,0) = fabs(nBkgSig - negFit[i]->Integral(105,135)/h->GetBinWidth(1));
    posRowM(0,i) = posColM(i,0);
    negRowM(0,i) = negColM(i,0);
    
  }

  TMatrixD posTemp = posRowM*(*corrMatrix*posColM);
  TMatrixD negTemp = negRowM*(*corrMatrix*negColM);
    
  *posUnc = TMath::Sqrt(posTemp(0,0));
  *negUnc = TMath::Sqrt(negTemp(0,0));

}

void alphaRplots(std::string outputFolder){

  setNCUStyle();
  gStyle->SetOptFit(0);
  gStyle->SetMarkerSize(0);
  gStyle->SetTitleSize(0.04,"XYZ");
  gStyle->SetLabelSize(0.03,"XYZ");
  gStyle->SetHistLineWidth(2);

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

  TH1D* h_sideTotalBKG  = addSamples(infiles,"ZprimeSide_pMC",f_DY100,f_DY200,f_DY400,f_DY600);
  TH1D* h_signTotalBKG  = addSamples(infiles,"ZprimeSign_pMC",f_DY100,f_DY200,f_DY400,f_DY600);
  TH1D* h_sideDATA      = addSamples(infiles,"ZprimeSide_pDA",f_DY100,f_DY200,f_DY400,f_DY600);
  TH1D* h_signDATA      = addSamples(infiles,"ZprimeSign_pDA",f_DY100,f_DY200,f_DY400,f_DY600);
  TH1D* h_corrPRmass    = addSamples(infiles,"corrPRmass_pDA",f_DY100,f_DY200,f_DY400,f_DY600);
  TH1D* h_corrPRmassAll = addSamples(infiles,"corrPRmassAll_pDA",f_DY100,f_DY200,f_DY400,f_DY600);

  h_sideTotalBKG->SetLineWidth(2);
  h_sideTotalBKG->SetLineColor(kBlack);
  h_sideTotalBKG->SetXTitle("ZH mass in side band of pseudo-MC");
  h_sideTotalBKG->SetYTitle("Event numbers");
  h_sideTotalBKG->SetTitleFont(62);

  h_signTotalBKG->SetLineWidth(2);
  h_signTotalBKG->SetLineColor(kBlack);
  h_signTotalBKG->SetXTitle("ZH mass in signal region of pseudo-MC");
  h_signTotalBKG->SetYTitle("Event numbers");
  h_signTotalBKG->SetTitleFont(62);

  h_signDATA->SetLineWidth(2);
  h_signDATA->SetLineColor(kBlue);
  h_signDATA->SetXTitle("ZH mass in signal region of pseudo-data");
  h_signDATA->SetYTitle("Event numbers");
  h_signDATA->SetTitleFont(62);

  h_corrPRmass->SetLineWidth(2);
  h_corrPRmass->SetLineColor(kBlack);
  h_corrPRmass->SetXTitle("Side band corrected pruned mass in pseudo-data");
  h_corrPRmass->SetYTitle("Event numbers");
  h_corrPRmass->SetTitleFont(62);

  h_corrPRmassAll->SetLineWidth(2);
  h_corrPRmassAll->SetLineColor(kBlack);
  h_corrPRmassAll->SetXTitle("Corrected pruned mass in pseudo-data");
  h_corrPRmassAll->SetYTitle("Event numbers");
  h_corrPRmassAll->SetTitleFont(62);

  // Make the statistics error more like data

  for( int i = 1; i <= nBins; i++ ){

    h_sideDATA->SetBinError(i,TMath::Sqrt(h_sideDATA->GetBinContent(i)));
    h_signDATA->SetBinError(i,TMath::Sqrt(h_signDATA->GetBinContent(i)));

  }

  for( int i = 1; i <= h_corrPRmass->GetNbinsX(); i++ ){

    h_corrPRmass   ->SetBinError(i,TMath::Sqrt(h_corrPRmass   ->GetBinContent(i)));
    h_corrPRmassAll->SetBinError(i,TMath::Sqrt(h_corrPRmassAll->GetBinContent(i)));

  }

  // Calculate alpha ratio

  TH1D* h_alphaRatio = new TH1D("h_alphaRatio", "", nBins, xmin, xmax); 

  h_alphaRatio->Sumw2();
  h_alphaRatio->SetXTitle("ZH mass");
  h_alphaRatio->SetYTitle("Alpha ratio");
  h_alphaRatio->Divide(h_signTotalBKG,h_sideTotalBKG);
  h_alphaRatio->SetMinimum(0);

  // Calculate number of backgrounds in signal region

  TH1D* h_numbkgDATA = (TH1D*)h_alphaRatio->Clone("h_numbkgDATA");

  h_numbkgDATA->Reset();

  for( int i = 1; i <= nBins; i++ ){

    double alphaRatio      = h_alphaRatio->GetBinContent(i); 
    double sideDATA        = h_sideDATA->GetBinContent(i);
    double numbkgDATA      = alphaRatio*sideDATA;      
    double alphaRatioError = h_alphaRatio->GetBinError(i);
    double sideDATAError   = h_sideDATA->GetBinError(i);

    if( alphaRatio == 0 || sideDATA == 0 ) continue;

    double numbkgDATAError = numbkgDATA*sqrt(pow((alphaRatioError/alphaRatio),2)+pow((sideDATAError/sideDATA),2));

    h_numbkgDATA->SetBinContent(i,numbkgDATA);
    h_numbkgDATA->SetBinError(i,numbkgDATAError);

  }

  TF1* f_fitPRmass = new TF1("f_fitPRmass", fitPRmass, 40, 240, 4);
  TF1* f_hollow_fitPRmass = new TF1("f_hollow_fitPRmass", hollow_fitPRmass, 40, 240, 5);
  TF1* f_fitZpmass = new TF1("f_fitZpmass", fitZpmass, xmin, xmax, 3);
  TF1* f_fitAlphaR = new TF1("f_fitAlphaR", divFunc, xmin, xmax, 6);

  f_fitPRmass->SetLineWidth(2);
  f_hollow_fitPRmass->SetLineWidth(2);
  f_fitZpmass->SetLineWidth(2);
  f_fitAlphaR->SetLineWidth(2);

  // Fit ZH mass and draw function of alpha ratio    

  double parAR[6] = {0};

  f_fitZpmass->SetParameters(0,0,0);
  h_signTotalBKG->Fit("f_fitZpmass", "", "", xmin, xmax);

  double chisqr_sgb = f_fitZpmass->GetChisquare();
  int ndf_sgb = f_fitZpmass->GetNDF();

  parAR[0] = f_fitZpmass->GetParameter(0);
  parAR[1] = f_fitZpmass->GetParameter(1);
  parAR[2] = f_fitZpmass->GetParameter(2);

  f_fitZpmass->SetParameters(0,0,0);
  h_sideTotalBKG->Fit("f_fitZpmass", "", "", xmin, xmax);

  double chisqr_sdb = f_fitZpmass->GetChisquare();
  int ndf_sdb = f_fitZpmass->GetNDF();

  parAR[3] = f_fitZpmass->GetParameter(0);
  parAR[4] = f_fitZpmass->GetParameter(1);
  parAR[5] = f_fitZpmass->GetParameter(2);

  f_fitAlphaR->SetParameters(parAR[0],parAR[1],parAR[2],parAR[3],parAR[4],parAR[5]);

  // Fit pruned mass

  double parFitPRm[4] = {1224,-0.107,139.6,107.4};

  f_fitPRmass->SetParameters(parFitPRm[0],parFitPRm[1],parFitPRm[2],parFitPRm[3],h_corrPRmassAll->GetBinWidth(1));
  f_fitPRmass->FixParameter(4,h_corrPRmassAll->GetBinWidth(1));
  f_fitPRmass->FixParameter(0,h_corrPRmassAll->Integral());
  
  h_corrPRmassAll->Fit("f_fitPRmass", "", "", 40, 240);
  h_corrPRmassAll->Fit("f_fitPRmass", "", "", 40, 240);

  double chisqr_cpma = f_fitPRmass->GetChisquare();
  int ndf_cpma = f_fitPRmass->GetNDF();

  // Fit pruned mass without signal region

  f_hollow_fitPRmass->SetParameters(parFitPRm[0],parFitPRm[1],parFitPRm[2],parFitPRm[3],h_corrPRmass->GetBinWidth(1));
  f_hollow_fitPRmass->FixParameter(4,h_corrPRmass->GetBinWidth(1));
  f_hollow_fitPRmass->FixParameter(0,h_corrPRmass->Integral());
  h_corrPRmass->Fit("f_hollow_fitPRmass", "", "", 40, 240);
  h_corrPRmass->Fit("f_hollow_fitPRmass", "", "", 40, 240);

  double chisqr_cpm = f_fitPRmass->GetChisquare();
  int ndf_cpm = f_fitPRmass->GetNDF();

  TFitResultPtr fitptr = h_corrPRmass->Fit(f_fitPRmass, "S");
  TFitResult fitresult = (*fitptr);
  TMatrixD corrMatrix  = fitresult.GetCorrelationMatrix();

  TGraphAsymmErrors* g_errorBands = fitUncertainty(f_fitPRmass, &corrMatrix);

  g_errorBands->SetFillStyle(3004);

  double nBkgSig = f_fitPRmass->Integral(105,135)/h_corrPRmass->GetBinWidth(1);
  double posUnc  = 0;
  double negUnc  = 0;

  fitUncNoBins(f_fitPRmass, &corrMatrix, h_corrPRmassAll, nBkgSig, &posUnc, &negUnc);

  cout << "\n*************************************************************************" << endl;
  cout << "** Number of backgrounds in signal region: "
       << nBkgSig << " + " << posUnc << " - " << negUnc << " **" << endl;
  cout << "*************************************************************************\n" << endl;

  h_numbkgDATA->Scale(nBkgSig/h_numbkgDATA->Integral(0,h_numbkgDATA->GetNbinsX()+1));

  // Output results

  TLegend* leg = new TLegend(0.21, 0.77, 0.87, 0.87);

  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);

  TLatex* lar = new TLatex();

  lar->SetTextSize(0.035);
  lar->SetLineWidth(5);

  TCanvas* c = new TCanvas("c","",0,0,1000,800);

  c->cd()->SetLogy(0);
  h_corrPRmassAll->Draw();
  lar->DrawLatexNDC(0.50, 0.70, Form("#chi^{2} / ndf: %f / %d",chisqr_cpma,ndf_cpma));
  lar->DrawLatexNDC(0.50, 0.60, "#font[22]{#color[4]{f(x) = #frac{1}{2} p_{0} e^{p_{1}x} ( 1 + erf ( #frac{x - p_{2}}{p_{3}} ) )}}");
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2015");
  lar->DrawLatexNDC(0.62, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  c->Print("alphaRatio.pdf(");

  c->cd()->SetLogy(0);
  h_corrPRmass->SetMaximum(300);
  h_corrPRmass->Draw();
  g_errorBands->Draw("3same");
  h_corrPRmass->Draw("same");
  leg->Clear();
  leg->AddEntry(h_corrPRmass, "Error = #sqrt{N_{per bin}}", "le");
  leg->AddEntry(g_errorBands, "Uncertainty based on fitting errors", "f");
  leg->Draw();
  lar->DrawLatexNDC(0.50, 0.65, Form("#chi^{2} / ndf: %f / %d",chisqr_cpm,ndf_cpm));
  lar->DrawLatexNDC(0.50, 0.55, "#font[22]{#color[4]{f(x) = #frac{1}{2} p_{0} e^{p_{1}x} ( 1 + erf ( #frac{x - p_{2}}{p_{3}} ) )}}");
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2015");
  lar->DrawLatexNDC(0.62, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  c->Print("alphaRatio.pdf");

  c->cd()->SetLogy();
  h_corrPRmass->SetMaximum(3e3);
  h_corrPRmass->Draw();
  g_errorBands->Draw("3same");
  h_corrPRmass->Draw("same");
  leg->Clear();
  leg->AddEntry(h_corrPRmass, "Error = #sqrt{N_{per bin}}", "le");
  leg->AddEntry(g_errorBands, "Uncertainty based on fitting errors", "f");
  leg->Draw();
  lar->DrawLatexNDC(0.25, 0.40, Form("#chi^{2} / ndf: %f / %d",chisqr_cpm,ndf_cpm));
  lar->DrawLatexNDC(0.25, 0.30, "#font[22]{#color[4]{f(x) = #frac{1}{2} p_{0} e^{p_{1}x} ( 1 + erf ( #frac{x - p_{2}}{p_{3}} ) )}}");
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2015");
  lar->DrawLatexNDC(0.62, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  c->Print("alphaRatio.pdf");
  
  c->cd()->SetLogy(0);
  h_signTotalBKG->Draw();
  lar->DrawLatexNDC(0.50, 0.80, Form("#chi^{2} / ndf: %f / %d",chisqr_sgb,ndf_sgb));
  lar->DrawLatexNDC(0.50, 0.70, "#font[22]{#color[4]{f_{signal}(x) = p_{0} e^{p_{1}x + #frac{p_{2}}{x}}}}");
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2015");
  lar->DrawLatexNDC(0.62, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  c->Print("alphaRatio.pdf");

  c->cd()->SetLogy(0);
  h_sideTotalBKG->Draw();
  lar->DrawLatexNDC(0.50, 0.80, Form("#chi^{2} / ndf: %f / %d",chisqr_sdb,ndf_sdb));
  lar->DrawLatexNDC(0.50, 0.70, "#font[22]{#color[4]{f_{side}(x) = p_{0} e^{p_{1}x + #frac{p_{2}}{x}}}}");
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2015");
  lar->DrawLatexNDC(0.62, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  c->Print("alphaRatio.pdf");

  c->cd()->SetLogy(0);
  h_alphaRatio->Draw();
  f_fitAlphaR->Draw("same");
  leg->Clear();
  leg->AddEntry(f_fitAlphaR, "#frac{f_{signal}(x)}{f_{side}(x)}", "l");
  leg->Draw();
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2015");
  lar->DrawLatexNDC(0.62, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  c->Print("alphaRatio.pdf");  

  c->Divide(1,2);
  TPad* c_up = (TPad*) c->GetListOfPrimitives()->FindObject("c_1");
  TPad* c_dw = (TPad*) c->GetListOfPrimitives()->FindObject("c_2"); 
  double up_height = 0.8;
  double dw_correction = 1.455;
  double dw_height = (1-up_height)*dw_correction;
  c_up->SetPad(0,1-up_height,1,1);
  c_dw->SetPad(0,0,1,dw_height);
  c_dw->SetBottomMargin(0.25);
  c_up->cd();
  h_signDATA->GetXaxis()->SetTitle("");
  h_signDATA->GetXaxis()->SetLabelOffset(999);
  h_signDATA->GetXaxis()->SetLabelSize(0);
  h_signDATA->Draw();
  h_numbkgDATA->Draw("same");
  lar->DrawLatexNDC(0.15, 0.94, "CMS preliminary 2015");
  lar->DrawLatexNDC(0.70, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");
  leg->Clear();
  leg->AddEntry(h_signDATA, "Truth backgrounds", "le");
  leg->AddEntry(h_numbkgDATA, "Predicted backgrounds", "le");
  leg->Draw();
  c_up->RedrawAxis();
  c_dw->cd();
  myRatio(h_numbkgDATA,h_signDATA);
  c->Draw();
  c->Print("alphaRatio.pdf)");

}
