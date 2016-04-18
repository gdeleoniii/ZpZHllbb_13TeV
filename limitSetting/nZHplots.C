#include <string>
#include <iostream>
#include <fstream>
#include <TH1.h>
#include <TFile.h>
#include "../setNCUStyle.h"
#include "../readHists.h"

void nZHplots(std::string channel, std::string rootfilename, std::string textfilename){

  setNCUStyle();

  std::string data1name, data2name;
 
  if( channel == "mu" ){
    data1name = "output/SingleElectron-Run2015D-V120151117_mZHmuSignal.root";
    data2name = "output/SingleElectron-Run2015D-V420151117_mZHmuSignal.root";
  }
 
  else if( channel == "ele" ){
    data1name = "output/SingleMuon-Run2015D-V120151117_mZHmuSignal.root";
    data2name = "output/SingleMuon-Run2015D-V420151117_mZHmuSignal.root";
  }
 
  readHist data1(data1name.data());
  readHist data2(data2name.data());
  readHist dy100("output/DYJetsToLL_M-50_HT-100to200_13TeV_mZHmuSignal.root");
  readHist dy200("output/DYJetsToLL_M-50_HT-200to400_13TeV_mZHmuSignal.root");
  readHist dy400("output/DYJetsToLL_M-50_HT-400to600_13TeV_mZHmuSignal.root");
  readHist dy600("output/DYJetsToLL_M-50_HT-600toInf_13TeV_mZHmuSignal.root");
  readHist tt   ("output/TT_TuneCUETP8M1_13TeV_mZHmuSignal.root");
  readHist ww   ("output/WW_TuneCUETP8M1_13TeV_mZHmuSignal.root");
  readHist wz   ("output/WZ_TuneCUETP8M1_13TeV_mZHmuSignal.root");
  readHist zz   ("output/ZZ_TuneCUETP8M1_13TeV_mZHmuSignal.root");
  readHist m800 ("output/ZprimeToZhToZlephbb_M-800_13TeV_mZHmuSignal.root");
  readHist m1000("output/ZprimeToZhToZlephbb_M-1000_13TeV_mZHmuSignal.root");
  readHist m1200("output/ZprimeToZhToZlephbb_M-1200_13TeV_mZHmuSignal.root");
  readHist m1400("output/ZprimeToZhToZlephbb_M-1400_13TeV_mZHmuSignal.root");
  readHist m1600("output/ZprimeToZhToZlephbb_M-1600_13TeV_mZHmuSignal.root");
  readHist m1800("output/ZprimeToZhToZlephbb_M-1800_13TeV_mZHmuSignal.root");
  readHist m2000("output/ZprimeToZhToZlephbb_M-2000_13TeV_mZHmuSignal.root");
  readHist m2500("output/ZprimeToZhToZlephbb_M-2500_13TeV_mZHmuSignal.root");
  readHist m3000("output/ZprimeToZhToZlephbb_M-3000_13TeV_mZHmuSignal.root");
  readHist m3500("output/ZprimeToZhToZlephbb_M-3500_13TeV_mZHmuSignal.root");
  readHist m4000("output/ZprimeToZhToZlephbb_M-4000_13TeV_mZHmuSignal.root");

  TH1D* h_Data = (TH1D*)(data1.getHist("mZprime"))->Clone("h_Data");
  h_Data->Reset();
  h_Data->Add(data1.getHist("mZprime"));
  h_Data->Add(data2.getHist("mZprime"));
    
  TH1D* h_DY = (TH1D*)(dy100.getHist("mZprime"))->Clone("h_DY");
  h_DY->Reset();
  h_DY->Add(dy100.getHist("mZprime"));
  h_DY->Add(dy200.getHist("mZprime"));
  h_DY->Add(dy400.getHist("mZprime"));
  h_DY->Add(dy600.getHist("mZprime"));
      
  TH1D* h_TT    = (TH1D*)(tt   .getHist("mZprime"));
  TH1D* h_WW    = (TH1D*)(ww   .getHist("mZprime"));
  TH1D* h_WZ    = (TH1D*)(wz   .getHist("mZprime"));
  TH1D* h_ZZ    = (TH1D*)(zz   .getHist("mZprime"));
  TH1D* h_M800  = (TH1D*)(m800 .getHist("mZprime"));
  TH1D* h_M1000 = (TH1D*)(m1000.getHist("mZprime"));
  TH1D* h_M1200 = (TH1D*)(m1200.getHist("mZprime"));
  TH1D* h_M1400 = (TH1D*)(m1400.getHist("mZprime"));
  TH1D* h_M1600 = (TH1D*)(m1600.getHist("mZprime"));
  TH1D* h_M1800 = (TH1D*)(m1800.getHist("mZprime"));
  TH1D* h_M2000 = (TH1D*)(m2000.getHist("mZprime"));
  TH1D* h_M2500 = (TH1D*)(m2500.getHist("mZprime"));
  TH1D* h_M3000 = (TH1D*)(m3000.getHist("mZprime"));
  TH1D* h_M3500 = (TH1D*)(m3500.getHist("mZprime"));
  TH1D* h_M4000 = (TH1D*)(m4000.getHist("mZprime"));
        
  TFile* outFile = new TFile(rootfilename.data(), "recreate");

  h_Data ->Write("data_obs");		 
  h_DY   ->Write("DYJETS");
  h_TT   ->Write("TT");
  h_WW   ->Write("WW");
  h_WZ   ->Write("WZ");
  h_ZZ   ->Write("ZZ");
  h_M800 ->Write("SIGM800");
  h_M1000->Write("SIGM1000");
  h_M1200->Write("SIGM1200");
  h_M1400->Write("SIGM1400");
  h_M1600->Write("SIGM1600");
  h_M1800->Write("SIGM1800");
  h_M2000->Write("SIGM2000");
  h_M2500->Write("SIGM2500");
  h_M3000->Write("SIGM3000");
  h_M3500->Write("SIGM3500");
  h_M4000->Write("SIGM4000");

  outFile->Write();
  
  fstream ftext;
  ftext.open(textfilename.data(), ios::out);

  ftext << "DATA\t"   << h_Data ->Integral() << "\n"; 
  ftext << "DYJETS\t" << h_DY   ->Integral() << "\n";
  ftext << "TTBAR\t"  << h_TT   ->Integral() << "\n";
  ftext << "WW\t"     << h_WW   ->Integral() << "\n";
  ftext << "WZ\t"     << h_WZ   ->Integral() << "\n";
  ftext << "ZZ\t"     << h_ZZ   ->Integral() << "\n";
  ftext << "M800\t"   << h_M800 ->Integral() << "\n";
  ftext << "M1000\t"  << h_M1000->Integral() << "\n";
  ftext << "M1200\t"  << h_M1200->Integral() << "\n";
  ftext << "M1400\t"  << h_M1400->Integral() << "\n";
  ftext << "M1600\t"  << h_M1600->Integral() << "\n";
  ftext << "M1800\t"  << h_M1800->Integral() << "\n";
  ftext << "M2000\t"  << h_M2000->Integral() << "\n";
  ftext << "M2500\t"  << h_M2500->Integral() << "\n";
  ftext << "M3000\t"  << h_M3000->Integral() << "\n";
  ftext << "M3500\t"  << h_M3500->Integral() << "\n";
  ftext << "M4000\t"  << h_M4000->Integral() << "\n";
  
  ftext.close();

}
