#include <TSystem.h>
#include "skimTree.C"

void runSkimTree(std::string channel){

  if( channel != "muon" && channel != "electron" ){

    std::cout << "Wrong channel!" << endl;
    return;

  }

  std::fstream samplePath("samplePath.txt");

  if( !samplePath.is_open() ){

    std::cout << "Text file does not exist!" << std::endl; 
    return;

  }

  std::string thisPath, keyWord;

  while( samplePath >> thisPath ){

    if( channel == "muon" )          keyWord = "SingleEle";
    else if( channel == "electron" ) keyWord = "SingleMuon";

    if( thisPath.find(keyWord.data()) == std::string::npos ){

      std::cout << "Now skim sample: " << thisPath << std::endl;

      skimTree skimthis(thisPath.data());
      skimthis.Loop(channel.data());

    }

  } // end of while

  std::cout << "All jobs are done!" << std::endl;

}
