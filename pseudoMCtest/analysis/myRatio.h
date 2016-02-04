#include <iostream>
#include <TH1.h>
#include <TLine.h>

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
