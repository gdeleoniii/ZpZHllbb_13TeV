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
  float crossSection(std::string);
   
};

readHist::readHist(std::string rootFile){

  std::string tmpStr = rootFile.erase(rootFile.find_last_not_of("/")+1);
  thisFileName = tmpStr.substr(tmpStr.find_last_of("/")+1);
  thisFile = TFile::Open(rootFile.data());

}

float readHist::crossSection(std::string token){

  std::ifstream textFile("textFile.txt");
  std::string thisPath;

  float crosssection = 0.;
  float thisNum = 0.;

  while( textFile >> thisPath >> thisNum ){

    if( thisPath.find(token) != std::string::npos )
      crosssection = thisNum;

  }

  return crosssection;
  
}

TH1D* readHist::getHist(std::string hname){

  float totalEvents      = ((TH1D*)(thisFile->Get("totalEvents")))->Integral();
  float thisCrossSection = crossSection(thisFileName.data());
  float thisScale        = 2080./(totalEvents/thisCrossSection); // dataLumi = 2080/pb

  TH1D* thisHist = (TH1D*)(thisFile->Get(Form("%s", hname.c_str())));  
  thisHist->Scale(thisScale);
  cout << thisScale << endl;  
  return thisHist;

}
