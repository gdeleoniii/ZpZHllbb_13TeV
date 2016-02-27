#include "skimTree.C"

void runSkimTree(){

  std::fstream samplePath("samplePath.txt");

  if( !samplePath.is_open() ){

    std::cout << "Text file does not exist!" << std::endl; 
    return;

  }

  while( !samplePath.eof() ){

    std::string thisPath;
    samplePath >> thisPath;

    std::cout << "Now skim sample: " << thisPath << std::endl;

    skimTree skimthis(thisPath.data());
    skimthis.Loop();

  }

  std::cout << "All jobs are done!" << std::endl;

}
