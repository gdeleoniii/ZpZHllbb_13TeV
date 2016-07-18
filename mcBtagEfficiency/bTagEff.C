

void bTagEff(string channel, int catcut, int flavorcut){

  // Input files and sum all backgrounds

  TChain* treeZjets = new TChain("tree");

  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-100to200_13TeV_MCbtagEff.root", channel.data()));
  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-200to400_13TeV_MCbtagEff.root", channel.data()));
  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-400to600_13TeV_MCbtagEff.root", channel.data()));
  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-600toInf_13TeV_MCbtagEff.root", channel.data()));

  int cat, jetFlavor;
  float jetPt, evWeight;

  treeZjets->SetBranchAddress("cat", &cat);
  treeZjets->SetBranchAddress("jetFlavor", &jetFlavor);
  treeZjets->SetBranchAddress("jetPt", &jetPt);
  treeZjets->SetBranchAddress("evweight", &evWeight);

  TH1F* h_jetPtnoBtag = new TH1F("h_jetPtnoBtag", "", 50, 0, 3000);
  TH1F* h_jetPtwtBtag = new TH1F("h_jetPtwtBtag", "", 50, 0, 3000);

  for(int ev = treeZjets->GetEntries()-1; ev >= 0; --ev){

    treeZjets->GetEntry(ev);

    //if( jetFlavor != flavorcut ) continue;

    if( cat == catcut )
      h_jetPtwtBtag->Fill(jetPt,evWeight);

    h_jetPtnoBtag->Fill(jetPt,evWeight);

  }

  TGraphAsymmErrors *g_bTagEff = new TGraphAsymmErrors();

  g_bTagEff->BayesDivide(h_jetPtwtBtag, h_jetPtnoBtag, "B");
  g_bTagEff->SetMarkerStyle(8);
  g_bTagEff->SetMaximum(1.3);
  g_bTagEff->GetYaxis()->SetTitle("Efficiency");  
  g_bTagEff->GetXaxis()->SetTitle("p_{T jet} [GeV]");

  string sflavor;

  if     ( flavorcut == 1 ) sflavor = "udsg";
  else if( flavorcut == 4 ) sflavor = "c";
  else if( flavorcut == 5 ) sflavor = "b";
  else                      sflavor = "0";

  TLegend leg(0.60, 0.70, 0.90, 0.87);

  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.04);
  leg.AddEntry(g_bTagEff, Form("%s %i b-tag %s flavor", channel.c_str(), catcut, sflavor.c_str()), "lp");
  
  TLatex lar;

  lar.SetNDC(kTRUE);
  lar.SetTextSize(0.04);
  lar.SetLineWidth(5);

  TCanvas c("c", "", 0, 0, 800, 600);

  c.cd();
  g_bTagEff->Draw("ap");
  leg.Draw();
  lar.DrawLatex(0.15, 0.83, "CMS");
  lar.DrawLatex(0.15, 0.79, "#it{#bf{Simulation}}");
  c.Print(Form("%s_MCbTagEff_%ibtag_%s.pdf", channel.data(), catcut, sflavor.c_str()));
  
}
