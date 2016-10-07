#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/uncertainties/rooFitUnc.h"

void getAlphaUnc(string channel, string catcut){

  gStyle->SetOptStat(0);

  string region[3] = {"central","up","down"};

  TF1 *f_alpha[3];
  TH1 *h_shape[3];

  for(int nw = 2; nw >= 0; --nw){

    rooFitUnc(channel.data(), catcut.data(), region[nw].data(), &f_alpha[nw], &h_shape[nw], nw);

  }

  // Calculate uncertainty of each mass bin

  float Mzh[13] = {750,800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4300};
  float Alpha[13], UncUp[13], UncDw[13];

  for( int im = 0; im < 13; ++im ){

    Alpha[im] = f_alpha[0]->Eval(Mzh[im]);
    UncUp[im] = fabs(f_alpha[1]->Eval(Mzh[im]) - Alpha[im]);
    UncDw[im] = fabs(f_alpha[2]->Eval(Mzh[im]) - Alpha[im]);
        
  } // end of mass points

  TGraphAsymmErrors *g_alpha = new TGraphAsymmErrors(13, Mzh, Alpha, 0, 0, UncDw, UncUp);

  g_alpha->SetTitle("");
  g_alpha->GetXaxis()->SetLimits(750,4300);
  g_alpha->GetXaxis()->SetTitle("m_{ZH}(GeV)");
  g_alpha->GetYaxis()->SetTitle("#alpha Ratio");  
  g_alpha->GetYaxis()->SetTitleOffset(1.3);
  g_alpha->SetMinimum(1e-2);
  g_alpha->SetMaximum(10);
  g_alpha->SetLineWidth(2);
  g_alpha->SetLineColor(kBlue);
  g_alpha->SetMarkerStyle(8);
  g_alpha->SetMarkerColor(kBlue);
  g_alpha->SetFillStyle(3002);

  // Store histograms in root file (for shape analysis)

  TFile f_shape(Form("background_bTag_%s_cat%s.root", channel.data(), catcut.data()), "recreate");

  h_shape[0]->Write("background_bTag");
  h_shape[1]->Write("background_bTagUp");
  h_shape[2]->Write("background_bTagDown");

  // Output the results

  TLatex lar;

  lar.SetTextSize(0.03);
  lar.SetLineWidth(5);

  TCanvas cv("cv","",0,0,1000,900);

  cv.cd()->SetLogy();

  g_alpha->Draw("Xac");
  g_alpha->Draw("3same");

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar.DrawLatexNDC(0.60, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s b-tag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, "b-tagging scale factor");

  cv.Draw();
  cv.Print(Form("alpha_bTagScale_%s_cat%s.pdf(", channel.data(), catcut.data()));

  cv.Clear();
  cv.cd()->SetLogy();

  h_shape[0]->SetTitle("");
  h_shape[0]->SetMinimum(1e-2);
  h_shape[0]->SetMaximum(10);
  h_shape[0]->Draw();
  h_shape[1]->SetLineColor(kRed);
  h_shape[1]->Draw("same");
  h_shape[2]->SetLineColor(kRed);
  h_shape[2]->Draw("same");

  cv.Draw();
  cv.Print(Form("alpha_bTagScale_%s_cat%s.pdf)", channel.data(), catcut.data()));

}
