#include <iostream>
#include <TF1.h>
#include <TFile.h>
#include <TAxis.h>
#include <TCanvas.h>

void drawEWfactor(){

  TFile* f = TFile::Open("skimSamples/scalefactors_v4.root");
  TF1* fewk_z = (TF1*)(f->Get("z_ewkcorr/z_ewkcorr_func"));
  TCanvas c("c","",0,0,1000,900);
  c.cd();
  fewk_z->SetTitle("");
  fewk_z->GetXaxis()->SetTitle("Z p_{T} (GeV)");
  fewk_z->GetYaxis()->SetTitle("EWK correction");
  fewk_z->GetYaxis()->SetTitleOffset(1.4);
  fewk_z->SetLineColor(61);
  fewk_z->Draw();
  c.Print("EWfactor.pdf");

}
