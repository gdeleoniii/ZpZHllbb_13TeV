#define skimTree_cxx
#include "skimTree.h"

void skimTree::Loop(string channel, TF1* fewk_z){

  if( fChain == 0 ) return;
  if( channel != "muon" && channel != "electron" ) return;

  // Set the name of output root file

  string tmpStr = inputFile_.erase(inputFile_.find_last_not_of("/")+1);
  string infix  = tmpStr.substr(tmpStr.find_last_of("/")+1); 
  string outputFile = ((channel == "muon") ? "skim_mu_" : "skim_ele_") + infix + ".root";

  // Now open new root file

  TFile* newfile_data = new TFile(outputFile.data(), "recreate");

  // Histogram to store total events 

  TH1D* h_totalEv = new TH1D("h_totalEv", "totalEvents", 2, -1, 1);

  // Clone tree and add a new branch

  Float_t ev_weight;

  fChain->LoadTree(0);

  TTree* newtree = fChain->GetTree()->CloneTree(0);
  newtree->Branch("ev_weight", &ev_weight, "ev_weight/F");

  cout << "Saving tree in " << outputFile << endl;

  Long64_t nentries = fChain->GetEntries();
  Long64_t nPassEv  = 0;

  for( Long64_t jentry = nentries-1; jentry >= 0; --jentry ){
        
    if( LoadTree(jentry) < 0 ) break;    

    fChain->GetEntry(jentry);

    if( (unsigned)jentry % 10000 == 0 ) fprintf(stdout, "Still left events %lli of %lli\n", jentry, nentries);
    
    // Apply MC weight, pile-up weight (Correct the pile-up shape of MC), and k factor weight for MC
    
    ev_weight = !isData ? ((mcWeight > 0 ? 1 : -1)*(puWeight((Int_t)pu_nTrueInt))*(infix.find("DYJets") != string::npos ? kWeight(fewk_z) : 1)) : 1;

    h_totalEv->Fill(0., ev_weight);

    // Remove event which is no hard interaction (noise)
    
    if( nVtx <= 0 ) continue;

    // Data filter (to filter non-collision bkg (ECAL/HCAL noise)) and trigger cut
      
    if( channel == "muon"     && !TriggerStatus("HLT_Mu45"  )) continue;
    if( channel == "electron" && !TriggerStatus("HLT_Ele105")) continue;
    if( isData && (
		   !FilterStatus("Flag_CSCTightHaloFilter") || 
		   !FilterStatus("Flag_eeBadScFilter"     ) || 
		   !FilterStatus("Flag_HBHENoiseFilter"   ) ||
		   !FilterStatus("Flag_HBHENoiseIsoFilter") )) continue;
    
    ++nPassEv;

    newtree->Fill();

  }

  h_totalEv->Write();
  newtree->AutoSave();

  delete newfile_data;

  gSystem->Exec("rm -f inputdir.txt");

  cout << "nentries = " << nentries << endl;
  cout << "Number of passed events = " << nPassEv << endl;
  cout << "Reduction rate = " << (float)nPassEv/(float)nentries << endl;

}
