#define skimTree_cxx
#include "skimTree.h"
#include "../dataFilter.h"
#include "../pileupMCweight.h"

void skimTree::Loop(){

  if( fChain == 0 ) return;

  // Set the name of output root file
  std::string tmpStr = inputFile_.erase(inputFile_.find_last_not_of("/")+1);
  std::string infix  = tmpStr.substr(tmpStr.find_last_of("/")+1); 
  std::string prefix = "skim_";
  std::string suffix = ".root";
  std::string outputFile = prefix + infix + suffix;

  // Now open new root file
  TFile* newfile_data = new TFile(outputFile.data(),"recreate");

  // Clone tree
  TTree* newtree = fChain->CloneTree(0);
  newtree->SetMaxTreeSize(5e9);
  std::cout << "Saving tree in " << outputFile << std::endl;

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nPassEvt = 0;
  //Long64_t nbytes = 0;
  //Long64_t nb = 0;

  Float_t eventWeight;
  //TBranch *branch_evWeight = newtree->Branch("eventWeight",&eventWeight,"eventWeight/F");

  for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
    
    Long64_t ientry = LoadTree(jentry);
    
    if( ientry < 0 ) break;
    //nb = fChain->GetEntry(jentry);
    //nbytes += nb;

    if( jentry % 10000 == 0 )
      fprintf(stderr, "Processing event %lli of %lli\n", jentry + 1, nentries);
    
    // Remove event which is no hard interaction (noise)

    //if( nVtx < 1 ) continue;

    // Correct the pile-up shape of MC

    //eventWeight = pileupWeight(isData, (Int_t)pu_nTrueInt);    
    //branch_evWeight->Fill();

    /*    
    // data filter (to filter non-collision bkg (ECAL/HCAL noise)) and trigger cut
      
    bool muTrigger  = TriggerStatus(data, "HLT_Mu45");
    bool eleTrigger = TriggerStatus(data, "HLT_Ele105");
    bool CSCT       = FilterStatus (data, "Flag_CSCTightHaloFilter");
    bool eeBadSc    = FilterStatus (data, "Flag_eeBadScFilter");
    bool Noise      = FilterStatus (data, "Flag_HBHENoiseFilter");
    bool NoiseIso   = FilterStatus (data, "Flag_HBHENoiseIsoFilter");

    if( !muTrigger || !eleTrigger ) continue;
    if( isData && !CSCT ) continue;
    if( isData && !eeBadSc ) continue;
    if( isData && !Noise ) continue;
    if( isData && !NoiseIso ) continue;

    */
    newtree->Fill();
    nPassEvt++;

  }

  newtree->AutoSave();
  delete newfile_data;

  std::cout << "nentries = " << nentries << std::endl;
  std::cout << "Number of passed events = " << nPassEvt << std::endl;
  std::cout << "Reduction rate = " << (double)nPassEvt/(double)nentries << std::endl;

}
