#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/uncertainties/rooFitUnc.h"

void getAlphaUnc(string channel, string catcut, string type, int first, int last){

  gStyle->SetOptStat(0);

  float Mzh[13] = {750,800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4300};

  int N = last-first;
  int iw = N-1;
  float alphaScale[13][N], alphaCentral;

  TF1 *f_alpha[N], *f_predict[N];
  TH1 *h_shape[N];

  for( int nw = last; nw >= first; --nw ){
    
    rooFitUnc(channel.data(), catcut.data(), "", &f_alpha[iw], &f_predict[iw], &h_shape[iw], nw, true);

     for( int im = 0; im < 13; ++im ){

      if( nw != first ) alphaScale[im][iw] = f_alpha[iw]->Eval(Mzh[im]);
      if( nw == first )	alphaCentral[im]   = f_alpha[iw]->Eval(Mzh[im]);

    }

    --iw;
  }

  // Calculate RMS value of each mass bin

  float Unc[13], relativeUnc[13];

  for( int im = 0; im < 13; ++im ){

    if( type == "mur1" )
      Unc[im] = ( fabs(alphaScale[im][1]-alphaCentral[im]) > fabs(alphaScale[im][0]-alphaCentral[im]) ) ?
	fabs(alphaScale[im][1]-alphaCentral[im]) : fabs(alphaScale[im][0]-alphaCentral[im]);

    else 
      Unc[im] = TMath::RMS(N, alphaScale[im]);

    relativeUnc[im] = Unc[im]/alphaCentral[im];

    fprintf(stdout, "massPoint=%i\trelativeUnc=%f\n", (int)Mzh[im], relativeUnc[im]);
 
  } // end of mass points

  TGraphErrors *g_alpha = new TGraphErrors(13, Mzh, alphaCentral, 0, Unc);

  g_alpha->SetTitle("");
  g_alpha->GetXaxis()->SetLimits(750,4300);
  g_alpha->GetXaxis()->SetTitle("m_{ZH}(GeV)");
  g_alpha->GetYaxis()->SetTitle("#alpha Ratio");  
  g_alpha->GetYaxis()->SetTitleOffset(1.3);
  g_alpha->SetMinimum(0.05);
  g_alpha->SetMaximum(50);
  g_alpha->SetLineWidth(2);
  g_alpha->SetLineColor(kBlue);
  g_alpha->SetMarkerStyle(8);
  g_alpha->SetMarkerColor(kBlue);
  g_alpha->SetFillStyle(3002);
  
  TLatex lar;

  lar.SetTextSize(0.03);
  lar.SetLineWidth(5);

  TCanvas cv("cv","",0,0,1000,900);

  cv.cd()->SetLogy();

  g_alpha->Draw("Xac");
  g_alpha->Draw("3same");

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.57, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s b-tag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, type=="pdf"?"PDF weight":"QCD scale factor");

  cv.Draw();
  cv.Print(Form("alpha_%sScale_%s_cat%s.pdf(", type.data(), channel.data(), catcut.data()));  

  cv.Clear();
  cv.cd()->SetLogy();

  f_predict[0]->SetTitle("");
  f_predict[0]->SetLineColor(kBlue);
  f_predict[0]->Draw();
  f_predict[1]->Draw("same");
  f_predict[2]->Draw("same");

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar.DrawLatexNDC(0.60, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s b-tag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, type=="pdf"?"PDF weight":"QCD scale factor");

  cv.Draw();
  cv.Print(Form("alpha_%sScale_%s_cat%s.pdf", type.data(), channel.data(), catcut.data()));

  cv.Clear();
  cv.cd()->SetLogy();

  h_shape[0]->SetTitle("");
  h_shape[0]->Draw();
  h_shape[1]->SetLineColor(kRed);
  h_shape[1]->Draw("same");
  h_shape[2]->SetLineColor(kRed);
  h_shape[2]->Draw("same");

  cv.Draw();
  cv.Print(Form("alpha_%sScale_%s_cat%s.pdf)", type.data(), channel.data(), catcut.data()));

}
