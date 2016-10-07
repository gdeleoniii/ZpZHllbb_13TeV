#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/readHists.h"

void nZHplots(string chan, string btag, string rootfilename, string textfilename){

  readHist data1(Form("output_%s_%sbtag/Single%s-Run2015D-v1_mZHLimit.root",              chan.data(),btag.data(),(chan=="ele") ? "Electron" : "Muon"));
  readHist data2(Form("output_%s_%sbtag/Single%s-Run2015D-v4_mZHLimit.root",              chan.data(),btag.data(),(chan=="ele") ? "Electron" : "Muon"));
  readHist dy100(Form("output_%s_%sbtag/DYJetsToLL_M-50_HT-100to200_13TeV_mZHLimit.root", chan.data(),btag.data()));
  readHist dy200(Form("output_%s_%sbtag/DYJetsToLL_M-50_HT-200to400_13TeV_mZHLimit.root", chan.data(),btag.data()));
  readHist dy400(Form("output_%s_%sbtag/DYJetsToLL_M-50_HT-400to600_13TeV_mZHLimit.root", chan.data(),btag.data()));
  readHist dy600(Form("output_%s_%sbtag/DYJetsToLL_M-50_HT-600toInf_13TeV_mZHLimit.root", chan.data(),btag.data()));
  readHist tt   (Form("output_%s_%sbtag/TT_TuneCUETP8M1_13TeV_mZHLimit.root",             chan.data(),btag.data()));
  readHist ww   (Form("output_%s_%sbtag/WW_TuneCUETP8M1_13TeV_mZHLimit.root",             chan.data(),btag.data()));
  readHist wz   (Form("output_%s_%sbtag/WZ_TuneCUETP8M1_13TeV_mZHLimit.root",             chan.data(),btag.data()));
  readHist zz   (Form("output_%s_%sbtag/ZZ_TuneCUETP8M1_13TeV_mZHLimit.root",             chan.data(),btag.data()));
  readHist zh   (Form("output_%s_%sbtag/ZH_HToBB_ZToLL_M125_13TeV_mZHLimit.root",         chan.data(),btag.data()));
  readHist m800 (Form("output_%s_%sbtag/ZprimeToZhToZlephbb_M-800_13TeV_mZHLimit.root",   chan.data(),btag.data()));
  readHist m1000(Form("output_%s_%sbtag/ZprimeToZhToZlephbb_M-1000_13TeV_mZHLimit.root",  chan.data(),btag.data()));
  readHist m1200(Form("output_%s_%sbtag/ZprimeToZhToZlephbb_M-1200_13TeV_mZHLimit.root",  chan.data(),btag.data()));
  readHist m1400(Form("output_%s_%sbtag/ZprimeToZhToZlephbb_M-1400_13TeV_mZHLimit.root",  chan.data(),btag.data()));
  readHist m1600(Form("output_%s_%sbtag/ZprimeToZhToZlephbb_M-1600_13TeV_mZHLimit.root",  chan.data(),btag.data()));
  readHist m1800(Form("output_%s_%sbtag/ZprimeToZhToZlephbb_M-1800_13TeV_mZHLimit.root",  chan.data(),btag.data()));
  readHist m2000(Form("output_%s_%sbtag/ZprimeToZhToZlephbb_M-2000_13TeV_mZHLimit.root",  chan.data(),btag.data()));
  readHist m2500(Form("output_%s_%sbtag/ZprimeToZhToZlephbb_M-2500_13TeV_mZHLimit.root",  chan.data(),btag.data()));
  readHist m3000(Form("output_%s_%sbtag/ZprimeToZhToZlephbb_M-3000_13TeV_mZHLimit.root",  chan.data(),btag.data()));
  readHist m3500(Form("output_%s_%sbtag/ZprimeToZhToZlephbb_M-3500_13TeV_mZHLimit.root",  chan.data(),btag.data()));
  readHist m4000(Form("output_%s_%sbtag/ZprimeToZhToZlephbb_M-4000_13TeV_mZHLimit.root",  chan.data(),btag.data()));

  string histName = "mZH";

  TH1D* h_Data   = (TH1D*)(data1.getHist(histName.data()))->Clone("h_Data");
  TH1D* h_Dom    = (TH1D*)(dy100.getHist(histName.data()))->Clone("h_Dom");
  TH1D* h_SubDom = (TH1D*)(tt.getHist(histName.data()))->Clone("h_SubDom");

  h_Data->Reset();
  h_Dom->Reset();
  h_SubDom->Reset();

  h_Data->Add(data1.getHist(histName.data()));
  h_Data->Add(data2.getHist(histName.data()));
    
  h_Dom->Add(dy100.getHist(histName.data()));
  h_Dom->Add(dy200.getHist(histName.data()));
  h_Dom->Add(dy400.getHist(histName.data()));
  h_Dom->Add(dy600.getHist(histName.data()));
      
  h_SubDom->Add(tt.getHist(histName.data()));
  h_SubDom->Add(ww.getHist(histName.data()));
  h_SubDom->Add(wz.getHist(histName.data()));
  h_SubDom->Add(zz.getHist(histName.data()));
  h_SubDom->Add(zh.getHist(histName.data()));

  TH1D* h_M800  = (TH1D*)(m800 .getHist(histName.data()));
  TH1D* h_M1000 = (TH1D*)(m1000.getHist(histName.data()));
  TH1D* h_M1200 = (TH1D*)(m1200.getHist(histName.data()));
  TH1D* h_M1400 = (TH1D*)(m1400.getHist(histName.data()));
  TH1D* h_M1600 = (TH1D*)(m1600.getHist(histName.data()));
  TH1D* h_M1800 = (TH1D*)(m1800.getHist(histName.data()));
  TH1D* h_M2000 = (TH1D*)(m2000.getHist(histName.data()));
  TH1D* h_M2500 = (TH1D*)(m2500.getHist(histName.data()));
  TH1D* h_M3000 = (TH1D*)(m3000.getHist(histName.data()));
  TH1D* h_M3500 = (TH1D*)(m3500.getHist(histName.data()));
  TH1D* h_M4000 = (TH1D*)(m4000.getHist(histName.data()));
        
  TFile* outFile = new TFile(rootfilename.data(), "recreate");

  h_Data  ->Write("data_obs");		 
  h_Dom   ->Write("DYJETS");
  h_SubDom->Write("SUBDOM");
  h_M800  ->Write("SIGM800");
  h_M1000 ->Write("SIGM1000");
  h_M1200 ->Write("SIGM1200");
  h_M1400 ->Write("SIGM1400");
  h_M1600 ->Write("SIGM1600");
  h_M1800 ->Write("SIGM1800");
  h_M2000 ->Write("SIGM2000");
  h_M2500 ->Write("SIGM2500");
  h_M3000 ->Write("SIGM3000");
  h_M3500 ->Write("SIGM3500");
  h_M4000 ->Write("SIGM4000");

  outFile->Write();
  
  fstream ftext;
  ftext.open(textfilename.data(), ios::out);

  ftext << "DATA\t"   << h_Data  ->Integral() << "\n"; 
  ftext << "DYJETS\t" << h_Dom   ->Integral() << "\n";
  ftext << "SUBDOM\t" << h_SubDom->Integral() << "\n";
  ftext << "M800\t"   << h_M800  ->Integral() << "\n";
  ftext << "M1000\t"  << h_M1000 ->Integral() << "\n";
  ftext << "M1200\t"  << h_M1200 ->Integral() << "\n";
  ftext << "M1400\t"  << h_M1400 ->Integral() << "\n";
  ftext << "M1600\t"  << h_M1600 ->Integral() << "\n";
  ftext << "M1800\t"  << h_M1800 ->Integral() << "\n";
  ftext << "M2000\t"  << h_M2000 ->Integral() << "\n";
  ftext << "M2500\t"  << h_M2500 ->Integral() << "\n";
  ftext << "M3000\t"  << h_M3000 ->Integral() << "\n";
  ftext << "M3500\t"  << h_M3500 ->Integral() << "\n";
  ftext << "M4000\t"  << h_M4000 ->Integral() << "\n";
  
  ftext.close();

}
