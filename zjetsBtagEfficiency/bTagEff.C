#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/untuplizer.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZee.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZmumu.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassJet.h"

void bTagEff(string inputFile, string outputFile, string channel){
  
  // read the ntuples (in pcncu)
  
  TreeReader data(inputFile.data());

  // Declare the histogram

  TFile f(inputFile.data());
  TH1D* h_totalEvents = (TH1D*)f.Get("h_totalEv");

  float varBinslight[] = {30,50,70,100,140,200,300,670,1000,2000};
  float varBins[] = {30,50,70,100,140,200,300,670,2000};

  int nvarBinslight = sizeof(varBinslight)/sizeof(varBinslight[0])-1;
  int nvarBins = sizeof(varBins)/sizeof(varBins[0])-1;

  TH1F* h_lJetPtnoCSV = new TH1F("h_lJetPtnoCSV", "lJetPtnoCSV", nvarBinslight, varBinslight);
  TH1F* h_lJetPtwtCSV = new TH1F("h_lJetPtwtCSV", "lJetPtwtCSV", nvarBinslight, varBinslight);
  TH1F* h_cJetPtnoCSV = new TH1F("h_cJetPtnoCSV", "cJetPtnoCSV", nvarBins, varBins);
  TH1F* h_cJetPtwtCSV = new TH1F("h_cJetPtwtCSV", "cJetPtwtCSV", nvarBins, varBins);
  TH1F* h_bJetPtnoCSV = new TH1F("h_bJetPtnoCSV", "bJetPtnoCSV", nvarBins, varBins);
  TH1F* h_bJetPtwtCSV = new TH1F("h_bJetPtwtCSV", "bJetPtwtCSV", nvarBins, varBins);

  // begin of event loop

  fprintf(stdout, "Total events %lli\n", data.GetEntriesFast());

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    if( (unsigned)ev % 100000 == 0 )
      fprintf(stdout, "Still left events %lli\n", ev);

    data.GetEntry(ev);

    Float_t        eventWeight     = data.GetFloat("ev_weight");
    TClonesArray*  muP4            = (TClonesArray*) data.GetPtrTObject("muP4");
    TClonesArray*  eleP4           = (TClonesArray*) data.GetPtrTObject("eleP4");
    TClonesArray*  FATjetP4        = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    Int_t          FATnJet         = data.GetInt("FATnJet");    
    Int_t*         FATnSubSDJet    = data.GetPtrInt("FATnSubSDJet");
    vector<float>* FATsubjetSDCSV  = data.GetPtrVectorFloat("FATsubjetSDCSV", FATnJet);
    vector<float>* FATsubjetSDPx   = data.GetPtrVectorFloat("FATsubjetSDPx", FATnJet);
    vector<float>* FATsubjetSDPy   = data.GetPtrVectorFloat("FATsubjetSDPy", FATnJet);
    vector<float>* FATsubjetSDPz   = data.GetPtrVectorFloat("FATsubjetSDPz", FATnJet);
    vector<float>* FATsubjetSDE    = data.GetPtrVectorFloat("FATsubjetSDE", FATnJet);
    vector<int>*   FATsubjetFlavor = data.GetPtrVectorInt("FATsubjetSDHadronFlavor", FATnJet);

    // select good reco level events     
    // select good leptons
      
    vector<int> goodLepID;

    if( channel == "ele" && !isPassZee(data,goodLepID)   ) continue;
    if( channel == "mu"  && !isPassZmumu(data,goodLepID) ) continue;

    TLorentzVector* thisLep = (channel=="ele") ? (TLorentzVector*)eleP4->At(goodLepID[0]) : (TLorentzVector*)muP4->At(goodLepID[0]);
    TLorentzVector* thatLep = (channel=="ele") ? (TLorentzVector*)eleP4->At(goodLepID[1]) : (TLorentzVector*)muP4->At(goodLepID[1]);

    // select good FATjet

    int goodFATJetID = -1;

    if( !isPassJet(data, &goodFATJetID, thisLep, thatLep) ) continue;

    TLorentzVector* thisJet = (TLorentzVector*)FATjetP4->At(goodFATJetID);

    float mllbb;

    if( !noiseCleaning(thisLep, thatLep, thisJet, &mllbb) ) continue;

    // b-tag efficiency part

    for( int is = 0; is < FATnSubSDJet[goodFATJetID]; ++is ){
      
      TLorentzVector thisSubJet;

      thisSubJet.SetPxPyPzE(FATsubjetSDPx[goodFATJetID][is],
			    FATsubjetSDPy[goodFATJetID][is],
			    FATsubjetSDPz[goodFATJetID][is],
			    FATsubjetSDE[goodFATJetID][is]);

      if( FATsubjetFlavor[goodFATJetID][is] != 4 && FATsubjetFlavor[goodFATJetID][is] != 5 ){

	h_lJetPtnoCSV->Fill(thisSubJet.Pt(),eventWeight);
	
	if( FATsubjetSDCSV[goodFATJetID][is] > 0.605 )
	  h_lJetPtwtCSV->Fill(thisSubJet.Pt(),eventWeight);

      }

      else if( FATsubjetFlavor[goodFATJetID][is] == 4 ){

	h_cJetPtnoCSV->Fill(thisSubJet.Pt(),eventWeight);

	if( FATsubjetSDCSV[goodFATJetID][is] > 0.605 )
	  h_cJetPtwtCSV->Fill(thisSubJet.Pt(),eventWeight);

      }
 
      else if( FATsubjetFlavor[goodFATJetID][is] == 5 ){

	h_bJetPtnoCSV->Fill(thisSubJet.Pt(),eventWeight);

	if( FATsubjetSDCSV[goodFATJetID][is] > 0.605 )
	  h_bJetPtwtCSV->Fill(thisSubJet.Pt(),eventWeight);

      } // end of if-else subjet flavor
 
    } // end of subjet loop                           

  } // end of event loop
  
  fprintf(stdout, "Processed all events\n");

  TFile* outFile = new TFile(Form("%s_%sMCbtagEff.root",outputFile.data(),channel.data()), "recreate");

  h_lJetPtnoCSV->Write("lJetPtnoCSV");
  h_lJetPtwtCSV->Write("lJetPtwtCSV");
  h_cJetPtnoCSV->Write("cJetPtnoCSV");
  h_cJetPtwtCSV->Write("cJetPtwtCSV");
  h_bJetPtnoCSV->Write("bJetPtnoCSV");
  h_bJetPtwtCSV->Write("bJetPtwtCSV");
  h_totalEvents->Write("totalEvents");
  
  outFile->Write();

  delete outFile;

}
