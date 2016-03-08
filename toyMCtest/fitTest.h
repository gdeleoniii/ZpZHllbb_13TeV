#include <iostream>
#include <TF1.h>
#include <TH1.h>
#include <TMatrixD.h>
#include <TFitResult.h>
#include <TGraphAsymmErrors.h>
#include "fitFunction.h"
#include "fitUncertainty.h"

/*
  Steps to study fit bias and pull by generating n toy-MCs :

  1. Fefore doing anything, make the stat error of MC more like data.
  2. Fit pruned mass of (signal + sideband) by only free the shape parameters; using equation A.
  3. Fit pruned mass of (side band only) by using the shape parameters from (2), only free the amplitude; using equation B.
  4. Use the fit function from (3) to generate a lot of toy samples (using TF1::GetRandom()).
  5. Remove the signal region of these toy samples, fit them using equation B, only free the amplitude.
  6. Calculate bias and pull. For the pull, the denominator is the  uncertainty based on error matrix.
*/

void fitTest(TH1D* h_prmass, TH1D* h_prmass_hollow,
	     TGraphAsymmErrors** g_errorBands, TGraphAsymmErrors** g_errorBands_hollow, 
	     TH1D* h_bias, TH1D* h_upPull, TH1D* h_dwPull){
  
  /// Make the statistics error of MC more like data
  /*
  for( int i = 1; i <= h_prmass->GetNbinsX(); i++ ){

    h_prmass->SetBinError(i,TMath::Sqrt(h_prmass->GetBinContent(i)));
    h_prmass_hollow->SetBinError(i,TMath::Sqrt(h_prmass_hollow->GetBinContent(i)));

  }
  */
  /// Fit pruned mass of (signal + sideband) by only free the shape parameters

  TF1* f_fitprmass = new TF1("f_fitprmass", fitPRmass, 40, 240, 5);

  f_fitprmass->SetLineWidth(2);
  f_fitprmass->SetLineColor(kBlue);

  f_fitprmass->SetParameters(h_prmass->Integral(), -0.107, 139.6, 107.4, h_prmass->GetBinWidth(1));
  f_fitprmass->FixParameter(0, h_prmass->Integral());
  f_fitprmass->FixParameter(4, h_prmass->GetBinWidth(1));

  h_prmass->Fit("f_fitprmass", "Q", "", 40, 240);
  
  // Calculate the fitting uncertainties by using error matrix

  TFitResultPtr fitptr = h_prmass->Fit(f_fitprmass, "QS");
  TFitResult fitresult = (*fitptr);
  TMatrixD corrMatrix  = fitresult.GetCorrelationMatrix();  
  
  *g_errorBands = fitErrorBand(f_fitprmass, &corrMatrix, fitPRmass);

  (*g_errorBands)->SetFillStyle(1001);
  (*g_errorBands)->SetFillColor(kYellow);

  /// Fit pruned mass of (side band only) by using the shape parameters from above
  
  TF1* f_fitprmass_hollow = new TF1("f_fitprmass_hollow", hollow_fitPRmass, 40, 240, 5);

  f_fitprmass_hollow->SetLineWidth(2);
  f_fitprmass_hollow->SetLineColor(kBlue);

  f_fitprmass_hollow->FixParameter(1,f_fitprmass->GetParameter(1));
  f_fitprmass_hollow->FixParameter(2,f_fitprmass->GetParameter(2));
  f_fitprmass_hollow->FixParameter(3,f_fitprmass->GetParameter(3));
  f_fitprmass_hollow->FixParameter(4,h_prmass_hollow->GetBinWidth(1));

  h_prmass_hollow->Fit("f_fitprmass_hollow", "Q", "", 40, 240);

  // Calculate the fitting uncertainties by using error matrix

  TFitResultPtr fitptr_hollow = h_prmass_hollow->Fit(f_fitprmass_hollow, "QS");
  TFitResult fitresult_hollow = (*fitptr_hollow);
  TMatrixD corrMatrix_hollow  = fitresult_hollow.GetCorrelationMatrix();
  
  *g_errorBands_hollow = fitErrorBand(f_fitprmass_hollow, &corrMatrix_hollow, hollow_fitPRmass);

  (*g_errorBands_hollow)->SetFillStyle(1001);
  (*g_errorBands_hollow)->SetFillColor(kYellow);

  /// Fluctuate the pruned mass histogram to test the fitting results

  TH1D* h_fluc = (TH1D*)h_prmass->Clone("h_fluc");
  TH1D* h_fluc_hollow = (TH1D*)h_prmass->Clone("h_fluc_hollow");
  
  for( int ntoy = 0; ntoy < 1000; ntoy++ ){

    if( ntoy % 10 == 0 ) std::cout << "Process ntoy = " << ntoy+1 << std::endl;
  
    h_fluc->Reset();
    h_fluc_hollow->Reset();

    for( int idata = 0; idata < (int)(h_prmass->Integral()); idata++ )
      h_fluc->Fill(f_fitprmass->GetRandom(40,240));
     
    int nd1 = (105 - h_fluc->GetBinLowEdge(1))/h_fluc->GetBinWidth(1)+1;
    int nd2 = (135 - h_fluc->GetBinLowEdge(1))/h_fluc->GetBinWidth(1);

    double nSigHist = h_fluc->Integral(nd1,nd2);

    for( int nbin = 1; nbin <= h_fluc->GetNbinsX(); nbin++ ){

      if( h_fluc_hollow->GetBinLowEdge(nbin) < 65 || h_fluc_hollow->GetBinLowEdge(nbin) > 140 ){

	h_fluc_hollow->SetBinContent(nbin,h_fluc->GetBinContent(nbin));
	h_fluc_hollow->SetBinError(nbin,h_fluc->GetBinError(nbin));

      }

    }

    TF1* f_fluc = new TF1("f_fluc", hollow_fitPRmass, 40, 240, 5);

    f_fluc->FixParameter(1,f_fitprmass->GetParameter(1));
    f_fluc->FixParameter(2,f_fitprmass->GetParameter(2));
    f_fluc->FixParameter(3,f_fitprmass->GetParameter(3));
    f_fluc->FixParameter(4,h_fluc_hollow->GetBinWidth(1));

    h_fluc_hollow->Fit("f_fluc", "Q", "", 40, 240);

    double nSigFit = f_fluc->Integral(105.0,135.0)/h_fluc_hollow->GetBinWidth(1);

    h_bias->Fill((nSigFit - nSigHist)/nSigHist);

    double upfitUnc, dwfitUnc;

    TFitResultPtr fitptr_fluc = h_fluc_hollow->Fit(f_fluc, "QS");
    TFitResult fitresult_fluc = (*fitptr_fluc);
    TMatrixD corrMatrix_fluc  = fitresult_fluc.GetCorrelationMatrix();

    fitUncertainty(f_fluc, &corrMatrix_fluc, hollow_fitPRmass, h_fluc_hollow, nSigHist, &upfitUnc, &dwfitUnc);

    h_upPull->Fill((nSigFit - nSigHist)/upfitUnc);
    h_dwPull->Fill((nSigFit - nSigHist)/dwfitUnc);

  } // end of ntoy loop
  
}
