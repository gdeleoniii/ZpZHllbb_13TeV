#define skimTree_cxx
#include "skimTree.h"
#include <TH1.h>

void skimTree::Loop(std::string channel){

  if( fChain == 0 ) return;
  if( channel != "muon" && channel != "electron" ) return;

  // Set the name of output root file

  std::string tmpStr = inputFile_.erase(inputFile_.find_last_not_of("/")+1);
  std::string infix  = tmpStr.substr(tmpStr.find_last_of("/")+1); 
  std::string prefix;

  if( channel == "muon" ) prefix = "skim_mu_";
  else if( channel == "electron" ) prefix = "skim_ele_";

  std::string suffix = ".root";
  std::string outputFile = prefix + infix + suffix;

  // Now open new root file

  TFile* newfile_data = new TFile(outputFile.data(), "recreate");

  // Histogram to store total events 

  TH1D* h_totalEv = new TH1D("h_totalEv", "totalEvents", 2, -1, 1);

  // Clone tree and add a new branch

  Float_t  ev_weight;

  fChain->LoadTree(0);

  TTree* newtree = fChain->GetTree()->CloneTree(0);
  newtree->Branch("ev_weight", &ev_weight, "ev_weight/F");

  std::cout << "Saving tree in " << outputFile << std::endl;

  Long64_t nentries  = fChain->GetEntries();
  Long64_t nPassEv   = 0;

  for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
        
    if( LoadTree(jentry) < 0 ) break;    

    fChain->GetEntry(jentry);

    if( jentry % 100000 == 0 )
      fprintf(stderr, "Processing event %lli of %lli\n", jentry + 1, nentries);
    
    // Corrections for MC 

    if( !isData ){
      
      Float_t mc_weight = mcWeight > 0 ? +1 : -1;
      Float_t pu_weight = pileupWeight((Int_t)pu_nTrueInt); // Correct the pile-up shape of MC

      ev_weight = mc_weight * pu_weight;

    }

    else ev_weight = 1;

    h_totalEv->Fill(0.0, ev_weight);

    // Remove event which is no hard interaction (noise)
    
    if( nVtx < 1 ) continue;

    // Data filter (to filter non-collision bkg (ECAL/HCAL noise)) and trigger cut
      
    bool muTrigger  = TriggerStatus("HLT_Mu45");
    bool eleTrigger = TriggerStatus("HLT_Ele105");
    bool CSCT       = FilterStatus ("Flag_CSCTightHaloFilter");
    bool eeBadSc    = FilterStatus ("Flag_eeBadScFilter");
    bool Noise      = FilterStatus ("Flag_HBHENoiseFilter");
    bool NoiseIso   = FilterStatus ("Flag_HBHENoiseIsoFilter");

    if( channel == "muon"     && !muTrigger  ) continue;
    if( channel == "electron" && !eleTrigger ) continue;
    if( isData && !CSCT     ) continue;
    if( isData && !eeBadSc  ) continue;
    if( isData && !Noise    ) continue;
    if( isData && !NoiseIso ) continue;
    
    nPassEv++;

    newtree->Fill();

  }

  h_totalEv->Write();
  newtree->AutoSave();

  delete newfile_data;
  gSystem->Exec("rm -f inputdir.txt");

  std::cout << "nentries = " << nentries << std::endl;
  std::cout << "Number of passed events = " << nPassEv << std::endl;
  std::cout << "Reduction rate = " << (double)nPassEv/(double)nentries << std::endl;

}
