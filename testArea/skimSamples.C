#include <string>
#include <TFile.h>
#include <TTree.h>

void skimSamples(std::string inputName){

  TFile* inputSample  = new TFile(inputName.data(),"READ");
  TFile* outputSample = new TFile(Form("%s_skim.root",inputName.c_str()),"RECREATE");

  std::string cut = "";

  TTree* fatTree  = (TTree*) inputSample->Get("tree");
  TTree* skimTree = fatTree->CopyTree(cut.data(), "", 1e10, 0);

  skimTree->Write();
  outputSample->Write();

}
