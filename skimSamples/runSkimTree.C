#include "skimTree.C"

void runSkimTree(string channel, string thisPath, int puScale){

  TFile* f = TFile::Open("scalefactors_v4.root");
  TF1* fewk_z = (TF1*)(f->Get("z_ewkcorr/z_ewkcorr_func"));
  string keyWord = (channel == "muon") ? "SingleElectron" : "SingleMuon";

  if( thisPath.find(keyWord.data()) == string::npos ){

    cout << "Now skim sample: " << thisPath << endl;
    skimTree skimthis(thisPath.data());
    skimthis.Loop(channel.data(),fewk_z,puScale);

  }

  cout << "Done!" << endl;

}
