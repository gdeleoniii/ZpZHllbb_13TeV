#include <string>
#include <iostream>
#include <TH1.h>
#include <TFile.h>

using namespace std;

double dataLumi      = 3000; //pb-1
double xSec_DY100    = 147.4*1.23;
double xSec_DY200    = 40.99*1.23;
double xSec_DY400    = 5.678*1.23;
double xSec_DY600    = 2.198*1.23;
double xSec_TT       = 831.76;
double xSec_WW       = 118.7;
double xSec_WZ       = 47.13;
double xSec_ZZ       = 16.523;
double xSec_ST_s     = 3.38;
double xSec_ST_t_at  = 26.49;
double xSec_ST_t_t   = 44.51;
double xSec_ST_tw_at = 35.85;
double xSec_ST_tw_t  = 35.85;

class addBackgrounds{

 public:
  TH1D* add(vector<string>& infiles, string hname, TFile* f_1);
  TH1D* add(vector<string>& infiles, string hname, TFile* f_1, TFile* f_2, TFile* f_3);
  TH1D* add(vector<string>& infiles, string hname, TFile* f_1, TFile* f_2, TFile* f_3, TFile* f_4);
  TH1D* add(vector<string>& infiles, string hname, TFile* f_1, TFile* f_2, TFile* f_3, TFile* f_4, TFile* f_5);
  TFile* getFile(string infiles, string hname, double crossSection, double* scale);
  
};

TFile* addBackgrounds::getFile(string infiles, string hname,
			       double crossSection, double* scale){

  TFile* f = TFile::Open(infiles.data());
  TH1D*  h = NULL;

  if( hname.find("pMC") != std::string::npos ) 
    h = (TH1D*)(f->Get("eventWeight_pMC"));
  else if( hname.find("pDA") != std::string::npos )
    h = (TH1D*)(f->Get("eventWeight_pDA"));

  *scale = dataLumi/(h->Integral()/crossSection);

  return f;

}

TH1D* addBackgrounds::add( vector<string>& infiles, 
			   string hname,
			   TFile* f_TT ){ 
			
  double scale_TT = 0;

  for(unsigned int i = 0; i < infiles.size(); i++){

    if( infiles[i].find("TT_") != std::string::npos ) f_TT = getFile(infiles[i].data(), hname.data(), xSec_TT, &scale_TT);

  }

  TH1D* TT_temp = (TH1D*)(f_TT->Get(Form("%s", hname.c_str())));

  TH1D* h_Total = (TH1D*)(f_TT->Get(Form("%s", hname.c_str())))->Clone("h_Total");

  h_Total->Reset();

  h_Total->Add(TT_temp, scale_TT);

  return h_Total;

}

TH1D* addBackgrounds::add( vector<string>& infiles, 
			   string hname,
			   TFile* f_WW,
			   TFile* f_WZ, 
			   TFile* f_ZZ ){

  double scale_WW = 0;
  double scale_WZ = 0;
  double scale_ZZ = 0;

  for(unsigned int i = 0; i < infiles.size(); i++){

    if( infiles[i].find("WW_") != std::string::npos ) f_WW = getFile(infiles[i].data(), hname.data(), xSec_WW, &scale_WW);
    if( infiles[i].find("WZ_") != std::string::npos ) f_WZ = getFile(infiles[i].data(), hname.data(), xSec_WZ, &scale_WZ);
    if( infiles[i].find("ZZ_") != std::string::npos ) f_ZZ = getFile(infiles[i].data(), hname.data(), xSec_ZZ, &scale_ZZ);

  }

  TH1D* WW_temp  = (TH1D*)(f_WW->Get(Form("%s", hname.c_str())));
  TH1D* WZ_temp  = (TH1D*)(f_WZ->Get(Form("%s", hname.c_str())));
  TH1D* ZZ_temp  = (TH1D*)(f_ZZ->Get(Form("%s", hname.c_str())));

  TH1D* h_Total = (TH1D*)(f_WW->Get(Form("%s", hname.c_str())))->Clone("h_Total");

  h_Total->Reset();

  h_Total->Add(WW_temp, scale_WW);
  h_Total->Add(WZ_temp, scale_WZ);
  h_Total->Add(ZZ_temp, scale_ZZ);

  return h_Total;

}

TH1D* addBackgrounds::add( vector<string>& infiles, 
			   string hname,
			   TFile* f_DY100, 
			   TFile* f_DY200, 
			   TFile* f_DY400, 
			   TFile* f_DY600 ){

  double scale_DY100 = 0;
  double scale_DY200 = 0;
  double scale_DY400 = 0;
  double scale_DY600 = 0;

  for(unsigned int i = 0; i < infiles.size(); i++){

    if( infiles[i].find("HT-100") != std::string::npos ) f_DY100 = getFile(infiles[i].data(), hname.data(), xSec_DY100, &scale_DY100);
    if( infiles[i].find("HT-200") != std::string::npos ) f_DY200 = getFile(infiles[i].data(), hname.data(), xSec_DY200, &scale_DY200);
    if( infiles[i].find("HT-400") != std::string::npos ) f_DY400 = getFile(infiles[i].data(), hname.data(), xSec_DY400, &scale_DY400);
    if( infiles[i].find("HT-600") != std::string::npos ) f_DY600 = getFile(infiles[i].data(), hname.data(), xSec_DY600, &scale_DY600);

  }

  TH1D* DY100_temp = (TH1D*)(f_DY100->Get(Form("%s", hname.c_str())));
  TH1D* DY200_temp = (TH1D*)(f_DY200->Get(Form("%s", hname.c_str())));
  TH1D* DY400_temp = (TH1D*)(f_DY400->Get(Form("%s", hname.c_str())));
  TH1D* DY600_temp = (TH1D*)(f_DY600->Get(Form("%s", hname.c_str())));

  TH1D* h_Total = (TH1D*)(f_DY100->Get(Form("%s", hname.c_str())))->Clone("h_Total");

  h_Total->Reset();

  h_Total->Add(DY100_temp, scale_DY100);
  h_Total->Add(DY200_temp, scale_DY200);
  h_Total->Add(DY400_temp, scale_DY400);
  h_Total->Add(DY600_temp, scale_DY600);

  return h_Total;

}

TH1D* addBackgrounds::add( vector<string>& infiles, 
			   string hname,
			   TFile* f_ST_s,
			   TFile* f_ST_t_at, 
			   TFile* f_ST_t_t,
			   TFile* f_ST_tw_at,
			   TFile* f_ST_tw_t ){ 

  double scale_ST_s     = 0;
  double scale_ST_t_at  = 0;
  double scale_ST_t_t   = 0;
  double scale_ST_tw_at = 0;
  double scale_ST_tw_t  = 0;

  for(unsigned int i = 0; i < infiles.size(); i++){

    if( infiles[i].find("s_4f")          != std::string::npos ) f_ST_s     = getFile(infiles[i].data(), hname.data(), xSec_ST_s,     &scale_ST_s);
    if( infiles[i].find("t_antitop_4f")  != std::string::npos ) f_ST_t_at  = getFile(infiles[i].data(), hname.data(), xSec_ST_t_at,  &scale_ST_t_at);
    if( infiles[i].find("t_top_4f")      != std::string::npos ) f_ST_t_t   = getFile(infiles[i].data(), hname.data(), xSec_ST_t_t,   &scale_ST_t_t);
    if( infiles[i].find("tW_antitop_5f") != std::string::npos ) f_ST_tw_at = getFile(infiles[i].data(), hname.data(), xSec_ST_tw_at, &scale_ST_tw_at);
    if( infiles[i].find("tW_top_5f")     != std::string::npos ) f_ST_tw_t  = getFile(infiles[i].data(), hname.data(), xSec_ST_tw_t,  &scale_ST_tw_t);

  }

  TH1D* s_4f_temp          = (TH1D*)(f_ST_s    ->Get(Form("%s", hname.c_str())));
  TH1D* t_antitop_4f_temp  = (TH1D*)(f_ST_t_at ->Get(Form("%s", hname.c_str())));
  TH1D* t_top_4f_temp      = (TH1D*)(f_ST_t_t  ->Get(Form("%s", hname.c_str())));
  TH1D* tW_antitop_5f_temp = (TH1D*)(f_ST_tw_at->Get(Form("%s", hname.c_str())));
  TH1D* tW_top_5f_temp     = (TH1D*)(f_ST_tw_t ->Get(Form("%s", hname.c_str())));

  TH1D* h_Total = (TH1D*)(f_ST_s->Get(Form("%s", hname.c_str())))->Clone("h_Total");

  h_Total->Reset();

  h_Total->Add(s_4f_temp,          scale_ST_s);
  h_Total->Add(t_antitop_4f_temp,  scale_ST_t_at);
  h_Total->Add(t_top_4f_temp,      scale_ST_t_t);
  h_Total->Add(tW_antitop_5f_temp, scale_ST_tw_at);
  h_Total->Add(tW_top_5f_temp,     scale_ST_tw_t);

  return h_Total;

}
