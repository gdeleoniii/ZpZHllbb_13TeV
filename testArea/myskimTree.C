#define skimTree_cxx
#include "skimTree.h"
#include "../dataFilter.h"
#include "../pileupMCweight.h"

void skimTree::Loop(){

  if(fChain == 0) return;

  std::string remword = ".root";
  size_t pos = inputFile_.find(remword);
  std::string forOutput = inputFile_;  

  if( pos != std::string::npos )
    forOutput.swap(forOutput.erase(pos,remword.length()));   

  std::string endfix = "_filtered.root";
  std::string outputFile = forOutput + endfix;

  // now open new root file
  TFile* newfile_data = new TFile(outputFile.data(),"recreate");

  // clone tree
  TTree* newtree = fChain->CloneTree(0);
  newtree->SetMaxTreeSize(5e9);
  std::cout << "Saving " << outputFile << " tree" << std::endl;

  std::ofstream fout;
  fout.open("wrong.dat");
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nPassEvt = 0;
  Long64_t nbytes = 0;
  Long64_t nb = 0;

  TBranch *b_evWeight = newtree->Branch("eventWeight",&eventWeight,"eventWeight/D");
  //newtree->SetBranchAddress("eventWeight",&eventWeight);

  for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
    
    Long64_t ientry = LoadTree(jentry);
    
    if( ientry < 0 ) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    
    // remove event which is no hard interaction (noise)

    if( nVtx < 1 ) continue;

    // Correct the pile-up shape of MC

    Double_t eventWeight = pileupWeight(isData, (Int_t)pu_nTrueInt);
    
    //h_eventWeight->Fill(0.,eventWeight);

    newtree->GetEntry(jentry);
    b_evWeight->Fill();
    newtree->Write();
    
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


    newtree->Fill();
    nPassEvt++;

  }

  newtree->AutoSave();
  delete newfile_data;
  fout.close();

  std::cout << "nentries = " << nentries << std::endl;
  std::cout << "Number of passed events = " << nPassEvt << std::endl;
  std::cout << "Reduction rate = " << (double)nPassEvt/(double)nentries << std::endl;

}
