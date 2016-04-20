#include <TSystem.h>
#include "skimTree.C"

void runSkimTree(std::string channel){

  std::fstream samplePath("samplePath.txt");

  if( !samplePath.is_open() ){

    std::cout << "Text file does not exist!" << std::endl; 
    return;

  }

  TFile* f = TFile::Open("scalefactors_v4.root");
  TF1* fewk_z = (TF1*)(f->Get("z_ewkcorr/z_ewkcorr_func"));

  std::string thisPath, keyWord;

  while( samplePath >> thisPath ){

    if     ( channel == "muon"     ) keyWord = "SingleEle";
    else if( channel == "electron" ) keyWord = "SingleMuon";

    if( thisPath.find(keyWord.data()) == std::string::npos ){

      std::cout << "Now skim sample: " << thisPath << std::endl;

      skimTree skimthis(thisPath.data());
      skimthis.Loop(channel.data(),fewk_z);

    }

  } // end of while

  std::cout << "All jobs are done!" << std::endl;

}
