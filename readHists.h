#include <string>
#include <iostream>
#include <fstream>
#include <TH1.h>
#include <TFile.h>

class readHist{
    
 public:
  
  readHist(string);
  TH1D* getHist(string);

 private:

  TFile* thisFile;
  string thisFileName;
  float crossSection(string);
   
};

readHist::readHist(string rootFile){

  string tmpStr = rootFile.erase(rootFile.find_last_not_of("/")+1);
  thisFileName = tmpStr.substr(tmpStr.find_last_of("/")+1);
  thisFile = TFile::Open(rootFile.data());

}

float readHist::crossSection(string thisPath){

  ifstream textFile("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/dataMCcompare/textFile.txt");
  string token;
  float crosssection = 0., thisNum = 0.;

  while( textFile >> token >> thisNum ){

    if( thisPath.find(token) != string::npos )
      crosssection = thisNum;

  }
  
  return crosssection;

}

TH1D* readHist::getHist(string hname){

  TH1D* thisHist = (TH1D*)(thisFile->Get(Form("%s", hname.c_str())));  

  thisHist->Scale(thisFileName.find("Run2015") != string::npos ? 1 : 2080.*crossSection(thisFileName.data())/((TH1D*)(thisFile->Get("totalEvents")))->Integral());

  return thisHist;

}
