R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagCalhead/BTagCalibrationStandalone_cpp.so)
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/untuplizer.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZee.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZmumu.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagCalhead/BTagCalibrationStandalone.h"

void mZHLimit(string inputFile, string outputFile, string channel, int cat, int mzh=0){

  // setup calibration and reader

  BTagCalibration calib("csvv1", "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/CSVV1.csv");

  BTagCalibrationReader reader_udsg(BTagEntry::OP_LOOSE, "central");
  BTagCalibrationReader reader_c(BTagEntry::OP_LOOSE, "central");
  BTagCalibrationReader reader_b(BTagEntry::OP_LOOSE, "central");

  reader_udsg.load(calib, BTagEntry::FLAV_UDSG, "mujets");
  reader_c.load(calib, BTagEntry::FLAV_C, "mujets");
  reader_b.load(calib, BTagEntry::FLAV_B, "mujets");

  // to read b-tag effinciency 

  TFile* f_l = TFile::Open(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/uncertainties/btagging.zjets/bTagEffroot/%s_udsgflavor_zjetsBtagEff.root", channel.data()));
  TFile* f_c = TFile::Open(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/uncertainties/btagging.zjets/bTagEffroot/%s_cflavor_zjetsBtagEff.root", channel.data()));
  TFile* f_b = mzh==0 ? 
    TFile::Open(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/uncertainties/btagging.zjets/bTagEffroot/%s_bflavor_zjetsBtagEff.root", channel.data())) :
    TFile::Open(Form("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/uncertainties/btagging.signal/bTagEffroot/%s_bflavor_signalBtagEff.root", channel.data()));
  
  TGraphAsymmErrors* g_l = (TGraphAsymmErrors*)(f_l->Get(Form("%s_udsgflavor", channel.data())));
  TGraphAsymmErrors* g_c = (TGraphAsymmErrors*)(f_c->Get(Form("%s_cflavor", channel.data())));
  TGraphAsymmErrors* g_b = mzh==0 ?
    (TGraphAsymmErrors*)(f_b->Get(Form("%s_bflavor", channel.data()))) :
    (TGraphAsymmErrors*)(f_b->Get(Form("%s_bflavor_m%i", channel.data(), mzh)));

  // read the ntuples (in pcncu)
  
  TreeReader data(inputFile.data());
  
  TFile f(inputFile.data());
  TH1F* h_totalEvents = (TH1F*)f.Get("h_totalEv");

  // Declare the histogram
     
  TH1F* h_mZprime = new TH1F("h_mZprime", "mZprime", 100, 400, 5000);

  h_mZprime->Sumw2();
  h_mZprime->GetXaxis()->SetTitle("mZprime");

  // begin of event loop

  fprintf(stdout, "Total events %lli\n", data.GetEntriesFast());

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    if( (unsigned)ev % 100000 == 0 )
      fprintf(stdout, "Still left events %lli\n", ev);

    data.GetEntry(ev);

    Bool_t         isData            = data.GetBool("isData");
    Float_t        eventWeight       = data.GetFloat("ev_weight");
    TClonesArray*  muP4              = (TClonesArray*) data.GetPtrTObject("muP4");
    TClonesArray*  eleP4             = (TClonesArray*) data.GetPtrTObject("eleP4");
    Int_t          FATnJet           = data.GetInt("FATnJet");    
    Int_t*         FATnSubSDJet      = data.GetPtrInt("FATnSubSDJet");
    Float_t*       FATjetPRmassCorr  = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    TClonesArray*  FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    vector<bool>&  FATjetPassIDLoose = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    vector<float>* FATsubjetSDCSV    = data.GetPtrVectorFloat("FATsubjetSDCSV", FATnJet);
    vector<float>* FATsubjetSDPx     = data.GetPtrVectorFloat("FATsubjetSDPx", FATnJet);
    vector<float>* FATsubjetSDPy     = data.GetPtrVectorFloat("FATsubjetSDPy", FATnJet);
    vector<float>* FATsubjetSDPz     = data.GetPtrVectorFloat("FATsubjetSDPz", FATnJet);
    vector<float>* FATsubjetSDE      = data.GetPtrVectorFloat("FATsubjetSDE", FATnJet);
    vector<int>*   FATsubjetFlavor   = data.GetPtrVectorInt("FATsubjetSDHadronFlavor", FATnJet);

    // select good leptons
      
    vector<int> goodLepID;

    if( channel == "ele" && !isPassZee(data,goodLepID)   ) continue;
    if( channel == "mu"  && !isPassZmumu(data,goodLepID) ) continue;

    TLorentzVector* thisLep = (channel=="ele") ? (TLorentzVector*)eleP4->At(goodLepID[0]) : (TLorentzVector*)muP4->At(goodLepID[0]);
    TLorentzVector* thatLep = (channel=="ele") ? (TLorentzVector*)eleP4->At(goodLepID[1]) : (TLorentzVector*)muP4->At(goodLepID[1]);

    // select good FATjet

    int goodFATJetID = -1;
    TLorentzVector thisJet(0,0,0,0);

    for( int ij = 0; ij < FATnJet; ++ij ){

      TLorentzVector* myJet = (TLorentzVector*)FATjetP4->At(ij);

      if( myJet->Pt() < 200 ) continue;
      if( fabs(myJet->Eta()) > 2.4 ) continue;
      if( !FATjetPassIDLoose[ij] ) continue;
      if( myJet->DeltaR(*thisLep) < 0.8 || myJet->DeltaR(*thatLep) < 0.8 ) continue;
      if( FATjetPRmassCorr[ij] < 105 || FATjetPRmassCorr[ij] > 135 ) continue;
      
      goodFATJetID = ij;
      thisJet = *myJet;

      break;
 
    } // end of FatnJet loop
 
    if( goodFATJetID < 0 ) continue;

    if( (*thisLep+*thatLep+thisJet).M() < 750 ) continue;

    int nsubBjet = 0;
    float pMC = 1., pData = 1.;

    for( int is = 0; is < FATnSubSDJet[goodFATJetID]; ++is ){
            
      TLorentzVector thisSubJet;
      
      thisSubJet.SetPxPyPzE(FATsubjetSDPx[goodFATJetID][is],
			    FATsubjetSDPy[goodFATJetID][is],
			    FATsubjetSDPz[goodFATJetID][is],
			    FATsubjetSDE[goodFATJetID][is]);

      if( !isData ){

	float btagEff, scaleFactor;

	if( FATsubjetFlavor[goodFATJetID][is] != 4 && FATsubjetFlavor[goodFATJetID][is] != 5 ){

	  btagEff = g_l->Eval(thisSubJet.Pt());
	  scaleFactor = reader_udsg.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, thisSubJet.Eta(), thisSubJet.Pt());

	}

	else if( FATsubjetFlavor[goodFATJetID][is] == 4 ){

	  btagEff = g_c->Eval(thisSubJet.Pt());
	  scaleFactor = reader_c.eval_auto_bounds("central", BTagEntry::FLAV_C, thisSubJet.Eta(), thisSubJet.Pt());

	}

	else if( FATsubjetFlavor[goodFATJetID][is] == 5 ){

	  btagEff = g_b->Eval(thisSubJet.Pt());
	  scaleFactor = reader_b.eval_auto_bounds("central", BTagEntry::FLAV_B, thisSubJet.Eta(), thisSubJet.Pt());

	}

	if( FATsubjetSDCSV[goodFATJetID][is] > 0.605 ){

	  pMC   *= btagEff;
	  pData *= scaleFactor * btagEff;
	  ++nsubBjet;
	
	}
      
	else{
	
	  pMC   *= (1 - btagEff);
	  pData *= (1 - scaleFactor * btagEff);
	
	}

      } // only for MC
      
    } // end of subjet for loop

    // b-tag cut

    if( cat == 1 && nsubBjet != 1 ) continue;
    if( cat == 2 && nsubBjet != 2 ) continue;
    
    h_mZprime->Fill((*thisLep+*thatLep+thisJet).M(), eventWeight*(pData/pMC));

  } // end of event loop

  fprintf(stderr, "Processed all events\n");

  TFile* outFile = new TFile(Form("%s_mZH%sLimit.root", outputFile.data(), channel.data()), "recreate");

  h_mZprime    ->Write("mZprime");
  h_totalEvents->Write("totalEvents");

  outFile->Write();

  delete outFile;
  
}