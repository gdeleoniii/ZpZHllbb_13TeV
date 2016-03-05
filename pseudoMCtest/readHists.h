#include <string>
#include <iostream>
#include <fstream>
#include <TH1.h>
#include <TFile.h>

class readHist{
    
 public:
  
  readHist(std::string);
  TH1D* getHist(std::string);

 private:

  TFile* thisFile;
  std::string thisFileName;
  double crossSection(std::string);
   
};

readHist::readHist(std::string rootFile){

  std::string tmpStr = rootFile.erase(rootFile.find_last_not_of("/")+1);
  thisFileName = tmpStr.substr(tmpStr.find_last_of("/")+1);
  thisFile = TFile::Open(rootFile.data());

}

double readHist::crossSection(std::string token){

  std::ifstream textFile("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/pseudoMCtest/textFile.txt");
  std::string thisPath;

  double crosssection = 0.;
  double thisNum = 0.;

  while( textFile >> thisPath >> thisNum ){

    if( thisPath.find(token) != std::string::npos )
      crosssection = thisNum;

  }

  return crosssection;
  
}

TH1D* readHist::getHist(std::string hname){

  double totalEvents      = ((TH1D*)(thisFile->Get("totalEvents")))->Integral();
  double thisCrossSection = crossSection(thisFileName.data());
  double thisScale        = 3000./(totalEvents/thisCrossSection); // dataLumi = 3000/pb

  TH1D* thisHist = (TH1D*)(thisFile->Get(Form("%s", hname.c_str())));  
  thisHist->Scale(thisScale);
  
  return thisHist;

}
