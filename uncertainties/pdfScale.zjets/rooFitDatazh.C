R__LOAD_LIBRARY(PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(PDFs/PdfDiagonalizer_cc.so)
using namespace RooFit;

void rooFitDatazh(string channel, string catcut, string type, int first, int last, int iter){

  // Suppress all the INFO message

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooMsgService::instance().setSilentMode(true);

  // Input files and sum all backgrounds

  TChain* treeZjets = new TChain("tree");

  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-100to200_13TeV_toyMC.root", channel.data()));
  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-200to400_13TeV_toyMC.root", channel.data()));
  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-400to600_13TeV_toyMC.root", channel.data()));
  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-600toInf_13TeV_toyMC.root", channel.data()));

  RooRealVar mZH("mllbb", "", 900., 3000.);
  
  RooPlot* alphaFrame = mZH.frame();
 
  const int NN = 1+(last-first)/iter;

  TH1 *hsb[NN];
  TH1 *hsg[NN];

  float binContent[21][NN];

  for(int nw = last; nw >= first; nw -= iter){

    // fprintf(stdout, ">>>> Weight %i <<<<\n", nw);

    // Define all the variables from the trees 

    RooRealVar cat("cat", "", 0, 2);
    RooRealVar mJet("prmass", "",  30.,  300.);
    RooRealVar evWeight(Form("evweight%02i",nw), "", -1.e10, 1.e10);
  
    // Set the range in zh mass 

    mZH.setRange("fullRange", 900., 3000.);

    RooArgSet variables(cat, mJet, mZH, evWeight);

    TCut catCut = Form("cat==%s", catcut.c_str());
    TCut sbCut  = "prmass>30 && !(prmass>65 && prmass<135) && prmass<300";
    TCut sigCut = "prmass>105 && prmass<135";

    // Create a dataset from a tree -> to process unbinned likelihood fitting

    RooDataSet dataSetZjetsSB("dataSetZjetsSB", "dataSetZjetsSB", variables, Cut(catCut && sbCut),  WeightVar(evWeight), Import(*treeZjets));
    RooDataSet dataSetZjetsSG("dataSetZjetsSG", "dataSetZjetsSG", variables, Cut(catCut && sigCut), WeightVar(evWeight), Import(*treeZjets));
  
    RooBinning b(21, 900., 3000.);

    //dataSetZjetsSG.plotOn(alphaFrame, Binning(b), MarkerColor((nw==first)?kBlue:kCyan), LineColor((nw==first)?kBlue:kCyan));

    hsb[nw] = dataSetZjetsSB.createHistogram(Form("hsb%02i",nw), mZH, Cut(catCut && sbCut),  WeightVar(evWeight), Import(*treeZjets), Binning(b));
    hsg[nw] = dataSetZjetsSG.createHistogram(Form("hsg%02i",nw), mZH, Cut(catCut && sigCut), WeightVar(evWeight), Import(*treeZjets), Binning(b));

    TH1F *hal = new TH1F("hal","", 21, 900, 3000);
    hal->Divide(hsg[nw],hsb[nw],1,1);
  
    RooDataHist dh("dh","dh", mZH, Import(*hal)) ;
    dh.plotOn(alphaFrame, Binning(b), MarkerColor((nw==first)?kBlue:kCyan), LineColor((nw==first)?kBlue:kCyan));

    for(int binx = 1; binx <= 21; ++binx){
      binContent[binx-1][nw] = hal->GetBinContent(binx);
    }

    hal->Delete();

  } // end of weight for loop

  gStyle->SetOptStat(0);

  float rms[21] = {0};

  TH1F* h_rms = new TH1F("h_rms","",21,900,3000);

  for(int binx = 0; binx < 21; ++binx){

    rms[binx] = TMath::RMS(NN, binContent[binx+1]);
    if(rms[binx] > 99) rms[binx] = -1;
    h_rms->SetBinContent(binx+1, rms[binx]);
    //fprintf(stdout, "bin %i : rms = %f\n", binx+1, rms[binx]);

  }

  h_rms->GetXaxis()->SetLabelSize(0.12);
  h_rms->GetXaxis()->SetTitle("m_{ZH} (GeV)");

  TLegend* leg = new TLegend(0.15,0.15,0.30,0.25);

  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(alphaFrame->findObject(alphaFrame->nameOf(0)), "central values", "l");
  leg->Draw();

  alphaFrame->addObject(leg);
  alphaFrame->SetTitle("");
  alphaFrame->GetXaxis()->SetLabelOffset(999);
  alphaFrame->GetYaxis()->SetTitle("#alpha Ratio");
  alphaFrame->GetYaxis()->SetTitleOffset(1.3);
  alphaFrame->SetMinimum(0.00001);

  TLatex* lar = new TLatex();

  lar->SetTextSize(0.035);
  lar->SetLineWidth(5);
  
  float up_height     = 0.8;
  float dw_correction = 1.375;
  float dw_height     = (1-up_height)*dw_correction;

  TCanvas c("c","",0,0,1000,900);
  c.Divide(1,2);

  TPad* c_up = (TPad*) c.GetListOfPrimitives()->FindObject("c_1");
  TPad* c_dw = (TPad*) c.GetListOfPrimitives()->FindObject("c_2"); 

  c_up->SetPad(0,1-up_height,1,1);
  c_dw->SetPad(0,0,1,dw_height);
  c_dw->SetBottomMargin(0.25);

  c_up->cd()->SetLogy(1);

  alphaFrame->Draw();
  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.75, 0.80, Form("%s  %s btag", channel.data(), catcut.data()));
  lar->DrawLatexNDC(0.75, 0.75, Form("%s", (first==0)?"mur = 1":type.data()));

  c_up->RedrawAxis();
  c_dw->cd();

  h_rms->Draw();
  lar->DrawLatexNDC(0.75, 0.80, "RMS of each bin");

  c.Draw();
  c.Print(Form("alpha_%sScale_%s_cat%s.pdf", type.data(), channel.data(), catcut.data()));

}
