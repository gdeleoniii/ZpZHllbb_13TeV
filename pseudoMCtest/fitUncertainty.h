#include <iostream>
#include <TF1.h>
#include <TH1.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TMatrixD.h>
#include <TGraphAsymmErrors.h>

/// for error band of fit curve /// 

TGraphAsymmErrors* fitErrorBand(const TF1* f, TMatrixD* corrMatrix, double (*fitFunc)(double*, double*)){

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

  TGraphAsymmErrors* g = new TGraphAsymmErrors(NBINS, funcX, funcY, 0, 0, negUnc, posUnc);

  return g;

}

/// for uncertainties of event number ///

void fitUncertainty(const TF1* f, TMatrixD* corrMatrix, double (*fitFunc)(double*, double*), 
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
