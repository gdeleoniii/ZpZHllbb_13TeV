#include <string>
#include <iostream>
#include <TH1.h>
#include <TFile.h>
/*
  double data_Lumi_     = 3000; //pb-1
  double xSec_DY100_    = 147.4*1.23;
  double xSec_DY200_    = 40.99*1.23;
  double xSec_DY400_    = 5.678*1.23;
  double xSec_DY600_    = 2.198*1.23;
  double xSec_TT_       = 831.76;
  double xSec_WW_       = 118.7;
  double xSec_WZ_       = 47.13;
  double xSec_ZZ_       = 16.523;
  double xSec_ST_s_     = 3.38;
  double xSec_ST_t_at_  = 26.49;
  double xSec_ST_t_t_   = 44.51;
  double xSec_ST_tw_at_ = 35.85;
  double xSec_ST_tw_t_  = 35.85;

  string root_Name_DY100_    = "DYJetsToLL_M-50_HT-100to200_13TeV_pseudoTest.root";
  string root_Name_DY200_    = "DYJetsToLL_M-50_HT-200to400_13TeV_pseudoTest.root";
  string root_Name_DY400_    = "DYJetsToLL_M-50_HT-400to600_13TeV_pseudoTest.root";
  string root_Name_DY600_    = "DYJetsToLL_M-50_HT-600toInf_13TeV_pseudoTest.root";
  string root_Name_TT_       = "TT_TuneCUETP8M1_13TeV_pseudoTest.root";
  string root_Name_WW_       = "WW_TuneCUETP8M1_13TeV_pseudoTest.root";
  string root_Name_WZ_       = "WZ_TuneCUETP8M1_13TeV_pseudoTest.root";
  string root_Name_ZZ_       = "ZZ_TuneCUETP8M1_13TeV_pseudoTest.root";
  string root_Name_ST_s_     = "ST_s_4f_leptonDecays_13TeV_pseudoTest.root";
  string root_Name_ST_t_at_  = "ST_t_antitop_4f_leptonDecays_13TeV_pseudoTest.root";
  string root_Name_ST_t_t_   = "ST_t_top_4f_leptonDecays_13TeV_pseudoTest.root";
  string root_Name_ST_tw_at_ = "ST_tW_antitop_5f_inclusiveDecays_13TeV_pseudoTest.root";
  string root_Name_ST_tw_t_  = "ST_tW_top_5f_inclusiveDecays_13TeV_pseudoTest.root";
*/

class addBkgs{

 private:

  std::string rootFile_;
  std::string hname_;
  double dataLumi_ = 3000.0;   //pb^-1
  double xSec_;
  double eventNum_;
  TFile* f_;
  TH1D*  h_this_;
  TH1D*  h_evWeight_;
  TH1D*  h_total_;
  
 public:

  addBkgs(std::string, double);
  double scale();
  TH1D* getHist(std::string);
  addBkgs operator + (const addBkgs &);

};

addBkgs::addBkgs(std::string rootFile, double xSec){

  rootFile_ = rootFile;
  xSec_ = xSec;
  f_ = TFile::Open(rootFile_.data());

}

double addBkgs::scale(){

  if( hname_.find("MC") != std::string::npos )
    h_evWeight_ = (TH1D*)(f_->Get("eventWeight_pMC"));
  
  else if( hname_.find("pDA") != std::string::npos )
    h_evWeight_ = (TH1D*)(f_->Get("eventWeight_pDA"));
  
  eventNum_ = h_evWeight_->Integral();

  return dataLumi_/(eventNum_/xSec_);

}

TH1D* addBkgs::getHist(std::string hname){

  hname_ = hname;
  
  h_this_  = (TH1D*)(f_->Get(Form("%s", hname_.c_str())));
  h_total_ = (TH1D*)(h_this_)->Clone("h_total");
  h_total_->Reset();
  h_total_->Add(h_this_, scale());

  return h_total_;
  
}

addBkgs addBkgs::operator + (const addBkgs &tmpBkg){

  addBkgs tmpObj(rootFile_, xSec_);

  TH1D* _this_h = tmpObj.getHist(hname_);
  _this_h->Add(tmpBkg.h_this_, scale());
  
  return tmpObj;

}
