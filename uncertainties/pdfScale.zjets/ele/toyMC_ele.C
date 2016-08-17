R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagCalhead/BTagCalibrationStandalone_cpp.so)
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <TH1D.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/untuplizer.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/isPassZee.h"
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagCalhead/BTagCalibrationStandalone.h"

float crossSection(string thisPath){

  ifstream textFile("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/xSec.txt");
  string token;
  float crosssection = 0., thisNum = 0.;

  while( textFile >> token >> thisNum ){

    if( thisPath.find(token) != string::npos )
      crosssection = thisNum;

  }

  return crosssection;

}

void toyMC_ele(string inputFile, string outputFile){

  // setup calibration and reader

  BTagCalibration calib("csvv1", "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/CSVV1.csv");

  BTagCalibrationReader reader_udsg(BTagEntry::OP_LOOSE, "central");
  BTagCalibrationReader reader_c(BTagEntry::OP_LOOSE, "central");
  BTagCalibrationReader reader_b(BTagEntry::OP_LOOSE, "central");

  reader_udsg.load(calib, BTagEntry::FLAV_UDSG, "mujets");
  reader_c.load(calib, BTagEntry::FLAV_C, "mujets");
  reader_b.load(calib, BTagEntry::FLAV_B, "mujets");

  // to read b-tag effinciency 

  TFile* f_l = TFile::Open("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagEffroot/ele_udsgflavor_zjetsBtagEff.root");
  TFile* f_c = TFile::Open("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagEffroot/ele_cflavor_zjetsBtagEff.root");
  TFile* f_b = TFile::Open("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/bTagEffroot/ele_bflavor_zjetsBtagEff.root");
  
  TGraphAsymmErrors* g_l = (TGraphAsymmErrors*)(f_l->Get("ele_udsgflavor"));
  TGraphAsymmErrors* g_c = (TGraphAsymmErrors*)(f_c->Get("ele_cflavor"));
  TGraphAsymmErrors* g_b = (TGraphAsymmErrors*)(f_b->Get("ele_bflavor"));

  // read the ntuples (in pcncu)

  TreeReader data(inputFile.data());

  TFile* f = new TFile(inputFile.data());
  TH1D* h_totalEvents = (TH1D*)f->Get("h_totalEv");

  // Create a tree to store variables

  TFile* outFile = new TFile(Form("%s_toyMC.root",outputFile.c_str()), "recreate");
  TTree* tree = new TTree("tree", "TreeForRooFit");

  Int_t   cat;
  Float_t mllbb, prmass;
  Float_t evweight[110] = {0};

  tree->Branch("cat",      &cat,     "cat/I");
  tree->Branch("mllbb",    &mllbb,   "mllbb/F");
  tree->Branch("prmass",   &prmass,  "prmass/F");

  for(int i = 0; i < 110; ++i)
    tree->Branch(Form("evweight%02i",i), &(evweight[i]), Form("evweight%02i/F",i));

  // Calculate the scale correspond to inputFile

  Float_t scale = 2512.*crossSection(outputFile.data())/h_totalEvents->Integral();

  // begin of event loop

  fprintf(stdout, "Total events %lli\n", data.GetEntriesFast());

  for( Long64_t ev = data.GetEntriesFast()-1; ev >= 0; --ev ){

    if( (unsigned)ev % 100000 == 0 )
      fprintf(stdout, "Still left events %lli\n", ev);

    data.GetEntry(ev);

    Float_t*       pdfscaleSysWeight = data.GetPtrFloat("pdfscaleSysWeights"); 
    Float_t        eventWeight       = data.GetFloat("ev_weight"); 
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
    if( !isPassZee(data,goodLepID) ) continue;

    TLorentzVector* thisLep = (TLorentzVector*)eleP4->At(goodLepID[0]);   
    TLorentzVector* thatLep = (TLorentzVector*)eleP4->At(goodLepID[1]);   

    // select good FATjet

    int goodFATJetID = -1;
    TLorentzVector thisJet(0,0,0,0);

    for( int ij = 0; ij < FATnJet; ++ij ){

      TLorentzVector* myJet = (TLorentzVector*)FATjetP4->At(ij);

      if( myJet->Pt() < 200 ) continue;
      if( fabs(myJet->Eta()) > 2.4 ) continue;
      if( !FATjetPassIDLoose[ij] ) continue;
      if( myJet->DeltaR(*thisLep) < 0.8 || myJet->DeltaR(*thatLep) < 0.8 ) continue;
      
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

      float thisSubJetPt = thisSubJet.Pt() < 420 ? thisSubJet.Pt() : 419.999;

      float btagEff, scaleFactor;

      if( FATsubjetFlavor[goodFATJetID][is] != 4 && FATsubjetFlavor[goodFATJetID][is] != 5 ){

	btagEff = g_l->Eval(thisSubJetPt);
	scaleFactor = reader_udsg.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, thisSubJet.Eta(), thisSubJetPt);

      }

      else if( FATsubjetFlavor[goodFATJetID][is] == 4 ){

	btagEff = g_c->Eval(thisSubJetPt);
	scaleFactor = reader_c.eval_auto_bounds("central", BTagEntry::FLAV_C, thisSubJet.Eta(), thisSubJetPt);

      }
      else if( FATsubjetFlavor[goodFATJetID][is] == 5 ){

	btagEff = g_b->Eval(thisSubJetPt);
	scaleFactor = reader_b.eval_auto_bounds("central", BTagEntry::FLAV_B, thisSubJet.Eta(), thisSubJetPt);

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
      
    } // end of subjet loop
    
    // b-tag cut
    
    if     ( nsubBjet == 1 ) cat = 1;
    else if( nsubBjet == 2 ) cat = 2;      
    else                     cat = 0;
    
    mllbb    = (*thisLep+*thatLep+thisJet).M(); 
    prmass   = FATjetPRmassCorr[goodFATJetID];

    for( int i = 0; i < 110; ++i ){
      evweight[i] = (i<9) ?
	eventWeight * scale * (pData/pMC) * pdfscaleSysWeight[i] :
	eventWeight * scale * (pData/pMC) * (pdfscaleSysWeight[i]/pdfscaleSysWeight[9]);
    }

    tree->Fill();

  } // end of event loop

  fprintf(stdout, "Processed all events\n");

  tree->Write();  
  outFile->Write();

  delete f;
  delete outFile;

}
