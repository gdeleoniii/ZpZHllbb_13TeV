#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <TH1D.h>
#include <TFile.h>
#include <TSystemDirectory.h>
#include "../setNCUStyle.h"

TFile* getFile(std::string infiles, Double_t crossSection, Double_t* scale){

  TFile* f = TFile::Open(infiles.data());
  TH1D* h = (TH1D*)(f->Get("eventWeight"));
  Double_t dataLumi = 3000; //pb-1 
  *scale = dataLumi/(h->Integral()/crossSection);

  return f;

}

void mZHSIGplots(std::string inputFolder, std::string rootfilename, std::string textfilename){

  setNCUStyle();

  std::vector<string> infiles;
 
  TSystemDirectory *base = new TSystemDirectory("root","root");
  base->SetDirectory(inputFolder.data());
  TList *listOfFiles = base->GetListOfFiles();
  TIter fileIt(listOfFiles);
  TFile *fileH = new TFile();
  Long64_t nfiles = 0;

  while( (fileH = (TFile*)fileIt()) ){
    
    std::string fileN = fileH->GetName();
    std::string baseString = "root";
    if( fileN.find(baseString) == std::string::npos ) continue;
    infiles.push_back(Form("%s/%s",inputFolder.data(),fileN.data()));
    nfiles++;
    
  }

  TFile *f_DY100 = NULL;
  TFile *f_DY200 = NULL;
  TFile *f_DY400 = NULL;
  TFile *f_DY600 = NULL;
  TFile *f_TTbar = NULL;
  TFile *f_WW    = NULL;
  TFile *f_WZ    = NULL;
  TFile *f_ZZ    = NULL;
  TFile *f_M800  = NULL;
  TFile *f_M1000 = NULL;
  TFile *f_M1200 = NULL;
  TFile *f_M1400 = NULL;
  TFile *f_M1600 = NULL;
  TFile *f_M1800 = NULL;
  TFile *f_M2000 = NULL;
  TFile *f_M2500 = NULL;
  TFile *f_M3000 = NULL;
  TFile *f_M3500 = NULL;
  TFile *f_M4000 = NULL;
  TFile *f_data0 = NULL;
  TFile *f_data1 = NULL;

  Double_t xSecDY100 = 147.4*1.23;
  Double_t xSecDY200 = 40.99*1.23;
  Double_t xSecDY400 = 5.678*1.23;
  Double_t xSecDY600 = 2.198*1.23;
  Double_t xSecTTbar = 831.76;
  Double_t xSecWW    = 118.7;
  Double_t xSecWZ    = 47.13;
  Double_t xSecZZ    = 16.523;

  std::fstream inFile("13TeV_xsec_Zhllbb.txt");

  std::vector<double> mass;
  std::vector<double> xSec;

  if( inFile.is_open() ){

    double mass_temp;
    double xSec_temp;
    
    while( !inFile.eof() ){
      
      inFile >> mass_temp >> xSec_temp;
      mass.push_back(mass_temp);
      xSec.push_back(xSec_temp);

    }

    inFile.close(); 
    
  }
   
  else cout << "Unable to open the file!" << endl;

  Double_t xSecM800  = 0;
  Double_t xSecM1000 = 0;
  Double_t xSecM1200 = 0;
  Double_t xSecM1400 = 0;
  Double_t xSecM1600 = 0;
  Double_t xSecM1800 = 0;
  Double_t xSecM2000 = 0;
  Double_t xSecM2500 = 0;
  Double_t xSecM3000 = 0;
  Double_t xSecM3500 = 0;
  Double_t xSecM4000 = 0;

  for(unsigned int i = 0; i < mass.size(); i++){

    if( mass[i] == 800 )  xSecM800 = xSec[i];
    if( mass[i] == 1000 ) xSecM1000 = xSec[i];
    if( mass[i] == 1200 ) xSecM1200 = xSec[i];
    if( mass[i] == 1400 ) xSecM1400 = xSec[i];
    if( mass[i] == 1600 ) xSecM1600 = xSec[i];
    if( mass[i] == 1800 ) xSecM1800 = xSec[i];
    if( mass[i] == 2000 ) xSecM2000 = xSec[i];
    if( mass[i] == 2500 ) xSecM2500 = xSec[i];
    if( mass[i] == 3000 ) xSecM3000 = xSec[i];
    if( mass[i] == 3500 ) xSecM3500 = xSec[i];
    if( mass[i] == 4000 ) xSecM4000 = xSec[i];

  }

  Double_t scaleDY100 = 0;
  Double_t scaleDY200 = 0;
  Double_t scaleDY400 = 0;
  Double_t scaleDY600 = 0;
  Double_t scaleTTbar = 0;
  Double_t scaleWW    = 0;
  Double_t scaleWZ    = 0;
  Double_t scaleZZ    = 0;
  Double_t scaleM800  = 0;
  Double_t scaleM1000 = 0;
  Double_t scaleM1200 = 0;
  Double_t scaleM1400 = 0;
  Double_t scaleM1600 = 0;
  Double_t scaleM1800 = 0;
  Double_t scaleM2000 = 0;
  Double_t scaleM2500 = 0;
  Double_t scaleM3000 = 0;
  Double_t scaleM3500 = 0;
  Double_t scaleM4000 = 0;
 
  Double_t dummy = -1;

  for(unsigned int i = 0; i < infiles.size(); i++){

    if( infiles[i].find("HT-100") != std::string::npos )    
      f_DY100 = getFile(infiles[i].data(), xSecDY100, &scaleDY100);

    if( infiles[i].find("HT-200") != std::string::npos )
      f_DY200 = getFile(infiles[i].data(), xSecDY200, &scaleDY200);

    if( infiles[i].find("HT-400") != std::string::npos )
      f_DY400 = getFile(infiles[i].data(), xSecDY400, &scaleDY400); 

    if( infiles[i].find("HT-600") != std::string::npos )
      f_DY600 = getFile(infiles[i].data(), xSecDY600, &scaleDY600);

    if( infiles[i].find("TT_") != std::string::npos )
      f_TTbar = getFile(infiles[i].data(), xSecTTbar, &scaleTTbar);

    if( infiles[i].find("WW_") != std::string::npos )
      f_WW    = getFile(infiles[i].data(), xSecWW, &scaleWW);

    if( infiles[i].find("WZ_") != std::string::npos )
      f_WZ    = getFile(infiles[i].data(), xSecWZ, &scaleWZ);

    if( infiles[i].find("ZZ_") != std::string::npos )
      f_ZZ    = getFile(infiles[i].data(), xSecZZ, &scaleZZ);

    if( infiles[i].find("M-800") != std::string::npos )
      f_M800  = getFile(infiles[i].data(), xSecM800 , &scaleM800 ); 

    if( infiles[i].find("M-1000") != std::string::npos )
      f_M1000 = getFile(infiles[i].data(), xSecM1000, &scaleM1000); 

    if( infiles[i].find("M-1200") != std::string::npos )
      f_M1200 = getFile(infiles[i].data(), xSecM1200, &scaleM1200); 

    if( infiles[i].find("M-1400") != std::string::npos )
      f_M1400 = getFile(infiles[i].data(), xSecM1400, &scaleM1400); 

    if( infiles[i].find("M-1600") != std::string::npos )
      f_M1600 = getFile(infiles[i].data(), xSecM1600, &scaleM1600); 

    if( infiles[i].find("M-1800") != std::string::npos )
      f_M1800 = getFile(infiles[i].data(), xSecM1800, &scaleM1800);

    if( infiles[i].find("M-2000") != std::string::npos )
      f_M2000 = getFile(infiles[i].data(), xSecM2000, &scaleM2000);

    if( infiles[i].find("M-2500") != std::string::npos )
      f_M2500 = getFile(infiles[i].data(), xSecM2500, &scaleM2500);

    if( infiles[i].find("M-3000") != std::string::npos )
      f_M3000 = getFile(infiles[i].data(), xSecM3000, &scaleM3000);

    if( infiles[i].find("M-3500") != std::string::npos )
      f_M3500 = getFile(infiles[i].data(), xSecM3500, &scaleM3500);

    if( infiles[i].find("M-4000") != std::string::npos )
      f_M4000 = getFile(infiles[i].data(), xSecM4000, &scaleM4000);

    if( infiles[i].find("V12015") != std::string::npos )
      f_data0 = getFile(infiles[i].data(), dummy, &dummy);

    if( infiles[i].find("V42015") != std::string::npos )
      f_data1 = getFile(infiles[i].data(), dummy, &dummy);

  }

  std::string h_name = "mZprime";

  TH1D* h_Data = (TH1D*)(f_data0->Get(h_name.data()))->Clone("h_Data");
  h_Data->Reset();
  h_Data->Add((TH1D*)(f_data0->Get(h_name.data())));
  h_Data->Add((TH1D*)(f_data1->Get(h_name.data())));
    
  TH1D* h_DY = (TH1D*)(f_DY100->Get(h_name.data()))->Clone("h_DY");
  h_DY->Reset();
  h_DY->Add((TH1D*)(f_DY100->Get(h_name.data())),scaleDY100);
  h_DY->Add((TH1D*)(f_DY200->Get(h_name.data())),scaleDY200);
  h_DY->Add((TH1D*)(f_DY400->Get(h_name.data())),scaleDY400);
  h_DY->Add((TH1D*)(f_DY600->Get(h_name.data())),scaleDY600);
      
  TH1D* h_TTbar = (TH1D*)(f_TTbar->Get(h_name.data()));
  TH1D* h_WW    = (TH1D*)(f_WW   ->Get(h_name.data()));
  TH1D* h_WZ    = (TH1D*)(f_WZ   ->Get(h_name.data()));
  TH1D* h_ZZ    = (TH1D*)(f_ZZ   ->Get(h_name.data()));
  TH1D* h_M800  = (TH1D*)(f_M800 ->Get(h_name.data())); 
  TH1D* h_M1000 = (TH1D*)(f_M1000->Get(h_name.data()));
  TH1D* h_M1200 = (TH1D*)(f_M1200->Get(h_name.data()));
  TH1D* h_M1400 = (TH1D*)(f_M1400->Get(h_name.data()));
  TH1D* h_M1600 = (TH1D*)(f_M1600->Get(h_name.data()));
  TH1D* h_M1800 = (TH1D*)(f_M1800->Get(h_name.data()));
  TH1D* h_M2000 = (TH1D*)(f_M2000->Get(h_name.data()));
  TH1D* h_M2500 = (TH1D*)(f_M2500->Get(h_name.data())); 
  TH1D* h_M3000 = (TH1D*)(f_M3000->Get(h_name.data()));
  TH1D* h_M3500 = (TH1D*)(f_M3500->Get(h_name.data()));
  TH1D* h_M4000 = (TH1D*)(f_M4000->Get(h_name.data()));
  
  h_TTbar->Scale(scaleTTbar);
  h_WW   ->Scale(scaleWW);
  h_WZ   ->Scale(scaleWZ);
  h_ZZ   ->Scale(scaleZZ);
  h_M800 ->Scale(scaleM800 );
  h_M1000->Scale(scaleM1000);
  h_M1200->Scale(scaleM1200);
  h_M1400->Scale(scaleM1400);
  h_M1600->Scale(scaleM1600);
  h_M1800->Scale(scaleM1800);
  h_M2000->Scale(scaleM2000);
  h_M2500->Scale(scaleM2500);
  h_M3000->Scale(scaleM3000);
  h_M3500->Scale(scaleM3500);
  h_M4000->Scale(scaleM4000);
      
  TFile* outFile = new TFile(rootfilename.data(), "recreate");

  h_Data ->Write("data_obs");		 
  h_DY   ->Write("DYJETS");
  h_TTbar->Write("TTBAR");
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
  ftext << "TTBAR\t"  << h_TTbar->Integral() << "\n";
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
