#define skimTree_cxx
#include "skimTree.h"

void skimTree::Loop(string channel, TF1* fewk_z, Int_t puScale){

  if( fChain == 0 ) return;
  if( channel != "muon" && channel != "electron" ) return;

  // Set the name of output root file

  string tmpStr = inputFile_.erase(inputFile_.find_last_not_of("/")+1);
  string infix  = tmpStr.substr(tmpStr.find_last_of("/")+1); 
  string outputFile = ((channel == "muon") ? "skim_mu_" : "skim_ele_") + infix + ".root";

  // Now open new root file

  TFile* newfile_data = new TFile(outputFile.data(), "recreate");

  // Histogram to store total events 

  TH1F h_totalEv("h_totalEv", "totalEvents", 1, -1, 1);
  TH1F h_genLepEv("h_genLepEv", "generatorEventsWtLepFlavor", 1, -1, 1);

  // Clone tree and add a new branch

  Float_t ev_weight;
  Bool_t isTrueGenFlavor;

  fChain->LoadTree(0);

  TTree* newtree = fChain->GetTree()->CloneTree(0);
  newtree->Branch("ev_weight", &ev_weight, "ev_weight/F");
  newtree->Branch("isTrueGenFlavor", &isTrueGenFlavor, "isTrueGenFlavor/O");

  Long64_t nPassEv = 0;

  fprintf(stdout, "Saving tree in %s\n", outputFile.c_str());
  fprintf(stdout, "Total events %lli\n", fChain->GetEntries());

  for( Long64_t jentry = fChain->GetEntries()-1; jentry >= 0; --jentry ){
        
    if( LoadTree(jentry) < 0 ) break;    

    fChain->GetEntry(jentry);

    // if( (unsigned)jentry % 100000 == 0 ) fprintf(stdout, "Still left events %lli\n", jentry);
    
    // Apply MC weight, pile-up weight (Correct the pile-up shape of MC), and k factor weight for MC
    
    ev_weight = !isData ? ((mcWeight > 0 ? 1 : -1)*(puWeight((Int_t)pu_nTrueInt,puScale))*(infix.find("DYJets") != string::npos ? kWeight(fewk_z) : 1)) : 1;

    h_totalEv.Fill(0., (!isData ? ((mcWeight > 0 ? 1 : -1)*(puWeight((Int_t)pu_nTrueInt))) : 1));

    // for signal MC only: save generated events according lepton flavor

    if( outputFile.find("ZprimeToZhToZlephbb") != string::npos ){

      Int_t lepId = 0;

      if( channel == "electron" )  lepId = 11;
      else if( channel == "muon" ) lepId = 13;

      isTrueGenFlavor = false;
 
      for( Int_t nGen = 0; nGen < nGenPar; ++nGen ){
      
	if( (*genParId)[nGen] == lepId && ( (*genMomParId)[nGen] == 23 || (*genMomParId)[nGen] == lepId ) && (*genParSt)[nGen] == 1 ){

	  h_genLepEv.Fill(0., (mcWeight > 0 ? 1 : -1)*(puWeight((Int_t)pu_nTrueInt)));
	  isTrueGenFlavor = true;
	  break;

	} // end of if statement
      
      } // end of generator events loop

    } // end of only signal MC
    
    // Remove event which is no hard interaction (noise)
    
    if( nVtx <= 0 ) continue;

    // Trigger cut
      
    if( channel == "muon"     && !TriggerStatus("HLT_Mu45"  ) ) continue;
    if( channel == "electron" && !TriggerStatus("HLT_Ele105") ) continue;

    // Data filter (to filter non-collision bkg (ECAL/HCAL noise))

    if( isData && !FilterStatus("Flag_CSCTightHaloFilter") ) continue;
    if( isData && !FilterStatus("Flag_eeBadScFilter"     ) ) continue;
    if( isData && !FilterStatus("Flag_HBHENoiseFilter"   ) ) continue;
    if( isData && !FilterStatus("Flag_HBHENoiseIsoFilter") ) continue;

    ++nPassEv;

    newtree->Fill();

  }

  h_genLepEv.Write();
  h_totalEv.Write();
  newtree->AutoSave();

  delete newfile_data;

  gSystem->Exec("rm -f inputdir.txt");

  fprintf(stdout, "Number of passed events %lli\n", nPassEv);
  fprintf(stdout, "Reduction rate = %f\n", (float)nPassEv/(float)fChain->GetEntries());

}
