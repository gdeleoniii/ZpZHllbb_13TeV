#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <TTree.h>
#include <TH1.h>
#include <TFile.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include "../../setNCUStyle.h"

TGraphAsymmErrors* Asymptotic(std::string foldername){

  setNCUStyle();
  gStyle->SetTitleSize(0.04,"XYZ");
  gStyle->SetLabelSize(0.03,"XYZ");

  ifstream xsect_file("../13TeV_xsec_Zh.txt", ios::in);

  float mH, CS;
  vector<int> v_mhxs;
  vector<float> v_xs, v_toterrh, v_toterrl;

  if( xsect_file.is_open() ){

    while( !xsect_file.eof() ){
      xsect_file >> mH >> CS;
      v_mhxs.push_back(mH);
      v_xs.push_back(CS);

    }

    xsect_file.close();

  }

  else{

    cout << "Failed to open text file: " << xsect_file << endl;
    return NULL;

  }


  int nmZH = v_mhxs.size();
  vector<double> v_mh, v_median, v_68l, v_68h, v_95l, v_95h, v_obs;
  TFile *fFREQ[nmZH];
  TTree *t[nmZH];

  for(int n = 0; n < nmZH; n++){

    char limitfile[100];

    sprintf(limitfile,"%s/higgsCombineCounting.Asymptotic.mZH%d.root",foldername.data(),v_mhxs[n]);

    fFREQ[n] = new TFile(limitfile, "READ");
    t[n] = (TTree*)fFREQ[n]->Get("limit");

    double mh, limit;
    float quant;

    t[n]->SetBranchAddress("mh", &mh);
    t[n]->SetBranchAddress("limit", &limit);
    t[n]->SetBranchAddress("quantileExpected", &quant);

    for(int i = 0; i < t[n]->GetEntries(); i++){

      t[n]->GetEntry(i);
 
      if (quant > -1.01 && quant < -0.99)
	v_obs.push_back(limit);

      else if (quant > 0.02 && quant < 0.03)
	v_95l.push_back(limit);
 
      else if (quant > 0.15 && quant < 0.17)
	v_68l.push_back(limit);
 
      else if (quant > 0.49 && quant < 0.51){
	v_median.push_back(limit);
        v_mh.push_back(mh);
      }
      else if (quant > 0.83 && quant < 0.85)
	v_68h.push_back(limit);
 
      else if (quant > 0.965 && quant < 0.98)
        v_95h.push_back(limit);
 
      else
	cout << "Error! Quantile =  " << quant << endl;
 
    }

  }                                                                                                                          

  double mass[nmZH];
  double medianD[nmZH];
  double up95err[nmZH], down95err[nmZH];
  int nMassEff = 0;

  for(int im = 0; im < nmZH; im++){

    double fl_xs = double(v_xs.at(im));
    fl_xs = (fl_xs);

    mass[nMassEff] = v_mhxs[im];
    medianD[nMassEff] = v_median.at(im) * fl_xs;

    up95err[nMassEff]   = (v_95h.at(im) - v_median.at(im)) * fl_xs;
    down95err[nMassEff] = (v_median.at(im) - v_95l.at(im)) * fl_xs;

    nMassEff++;

  } //end loop over im (mass points)

  TGraphAsymmErrors *gr95_cls = new TGraphAsymmErrors(nMassEff, mass, medianD, 0, 0, down95err, up95err);

  return gr95_cls;

}

void combineExpLimit(){

  TGraphAsymmErrors* g0 = Asymptotic("v1");
  TGraphAsymmErrors* g1 = Asymptotic("v2");

  TCanvas *cMCMC = new TCanvas("cMCMC", "", 1000, 800);

  cMCMC->cd();
  cMCMC->SetGridx(1);
  cMCMC->SetGridy(1);

  // draw a frame to define the range                                                                                                                                                  

  double fr_left = 0.0, fr_down = 1E-4, fr_right = 4500.0, fr_up = 10;

  TH1F *hr = cMCMC->DrawFrame(fr_left, fr_down, fr_right, fr_up, "");

  hr->SetXTitle("M_{ZH} [GeV]");
  hr->SetYTitle("#sigma_{95%} #times BR(Z'#rightarrow ZH) [pb]"); // #rightarrow 2l2q                                                                                                  
  hr->GetYaxis()->SetTitleSize(0.03);
  hr->GetYaxis()->SetTitleOffset(1.6);

  g0->SetLineColor(kRed);
  g0->SetLineWidth(3);
  g0->GetXaxis()->SetRangeUser(fr_left, fr_right);
  g0->Draw("LX");

  g1->SetLineColor(kBlue);
  g1->SetLineWidth(3);
  g1->GetXaxis()->SetRangeUser(fr_left, fr_right);
  g1->Draw("LXsame");


  TH1D* postGrid = new TH1D("postGrid", "", 1, fr_left, fr_right);
  postGrid->GetYaxis()->SetRangeUser(fr_down, fr_up);
  postGrid->Draw("AXIGSAME");

  TLegend *leg = new TLegend(.13, .20, .70, .35);

  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  leg->AddEntry(g0, "#splitline{1 subjet with csv>0.605(1 subjet if both csv>0.605 & dR<0.3) +}{2 subjet with csv>0.605 & dR>0.3}", "L");
  leg->AddEntry(g1, "1 subjet with csv>0.605 + 2 subjet with csv>0.605", "L");
  leg->SetTextSize(0.025);
  leg->Draw();

  TLatex * latex = new TLatex();

  latex->SetNDC(kTRUE);
  latex->SetTextSize(0.035);
  latex->DrawLatex(0.15, 0.94, "CMS preliminary 2015");
  latex->DrawLatex(0.65, 0.94, "L = 3 fb^{-1} at #sqrt{s} = 13 TeV");

  gPad->RedrawAxis("");
  cMCMC->Update();
  gPad->SetLogy();

  cMCMC->Print("use.pdf");

}
