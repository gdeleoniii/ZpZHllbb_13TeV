R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/PdfDiagonalizer_cc.so)
using namespace RooFit;

void rooFitAlpha(string channel, string catcut){

  // Suppress all the INFO message

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooMsgService::instance().setSilentMode(true);

  // Define all the variables from the trees

  RooRealVar cat ("cat", "", 0, 2);
  RooRealVar mJet("prmass", "", 30., 300.);
  RooRealVar evWeight("evweight", "", 0., 1.e3);
  RooRealVar mZH("mllbb", "M_{ZH}", 800., 4000., "GeV");
 
  mZH.setRange("fullRange", 800., 4000.);
 
  RooBinning mZHbin(21, 800., 4000.);
  RooArgSet variables(cat, mJet, mZH, evWeight);
 
  TCut catCut = Form("cat==%s", catcut.c_str());
  TCut sbCut  = "prmass>30 && !(prmass>65 && prmass<135) && prmass<300";
  TCut sigCut = "prmass>105 && prmass<135";

  string region[3] = {"central","up","down"};
  float alpha[11][3];

  TF1* f_alpha = new TF1("f_alpha", "[0]*TMath::Exp(-x/([1]+[2]*x))/TMath::Exp(-x/([3]+[4]*x))", 800, 4000);

  for(int nw = 2; nw >= 0; --nw){

    // Input files and sum all backgrounds

    TChain* treeZjets = new TChain("tree");
    
    treeZjets->Add(Form("Zjets/DYJetsToLL_M-50_HT-100to200_13TeV_%s_%sMiniTree.root", region[nw].data(), channel.data()));
    treeZjets->Add(Form("Zjets/DYJetsToLL_M-50_HT-200to400_13TeV_%s_%sMiniTree.root", region[nw].data(), channel.data()));
    treeZjets->Add(Form("Zjets/DYJetsToLL_M-50_HT-400to600_13TeV_%s_%sMiniTree.root", region[nw].data(), channel.data()));
    treeZjets->Add(Form("Zjets/DYJetsToLL_M-50_HT-600toInf_13TeV_%s_%sMiniTree.root", region[nw].data(), channel.data()));
    
    // Create a dataset from a tree -> to process an unbinned likelihood fitting

    RooDataSet dataSetZjetsSB("dataSetZjetsSB", "dataSetZjetsSB", variables, Cut(catCut && sbCut),  WeightVar(evWeight), Import(*treeZjets));  
    RooDataSet dataSetZjetsSG("dataSetZjetsSG", "dataSetZjetsSG", variables, Cut(catCut && sigCut), WeightVar(evWeight), Import(*treeZjets));

    // Total events number

    RooRealVar nSBMcEvents("nSBMcEvents", "nSBMcEvents", 0., 1.e10);
    RooRealVar nSGMcEvents("nSGMcEvents", "nSGMcEvents", 0., 1.e10);

    nSBMcEvents.setVal(dataSetZjetsSB.sumEntries());
    nSBMcEvents.setConstant(true);
  
    nSGMcEvents.setVal(dataSetZjetsSG.sumEntries());
    nSGMcEvents.setConstant(true);
  
    // Alpha ratio part

    // set fit parameters // [a,b][min,max]

    float sgVaraMin, sgVaraMax;

    if( channel == "ele" ){
      sgVaraMin = 10; sgVaraMax = 50;
    }
    else{ 
      sgVaraMin = 300; sgVaraMax = 350; 
    }

    RooRealVar sbVara("sbVara", "sbVara", 225., 150., 300.);
    RooRealVar sbVarb("sbVarb", "sbVarb", 0.025, 0.01, 0.10);
    RooRealVar sgVara("sgVara", "sgVara", 0.5*(sgVaraMin+sgVaraMax), sgVaraMin, sgVaraMax);
    RooRealVar sgVarb("sgVarb", "sgVarb", 0.05, 0.0001, 0.1);

    // Fit ZH mass in side band

    RooGenericPdf model_ZHSB("model_ZHSB", "model_ZHSB", "TMath::Exp(-@0/(@1+@2*@0))", RooArgSet(mZH,sbVara,sbVarb));
    RooExtendPdf ext_model_ZHSB("ext_model_ZHSB", "ext_model_ZHSB", model_ZHSB, nSBMcEvents);
    
    ext_model_ZHSB.fitTo(dataSetZjetsSB, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));

    // Fit ZH mass in signal region

    RooGenericPdf model_ZHSG("model_ZHSG", "model_ZHSG", "TMath::Exp(-@0/(@1+@2*@0))", RooArgSet(mZH,sgVara,sgVarb));
    RooExtendPdf ext_model_ZHSG("ext_model_ZHSG", "ext_model_ZHSG", model_ZHSG, nSGMcEvents);

    ext_model_ZHSG.fitTo(dataSetZjetsSG, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));

    // Set the model of alpha ratio

    float normConst = ((TF1*)ext_model_ZHSB.asTF(mZH, RooArgList(sbVara, sbVarb)))->Integral(800,4000)/((TF1*)ext_model_ZHSG.asTF(mZH, RooArgList(sgVara, sgVarb)))->Integral(800,4000);

    f_alpha->SetParameters(normConst, sgVara.getVal(), sgVarb.getVal(), sbVara.getVal(), sbVarb.getVal());

    int mzh = 800;
    for( int im = 0; im < 11; ++im ){
      alpha[im][nw] = f_alpha->Eval(mzh);
      mzh += (mzh<2000) ? 200 : 500;
    }

    fprintf(stdout, "sbVara=%f\tsbVarb=%f\tsgVara=%f\tsgVarb=%f\n", sbVara.getVal(), sbVarb.getVal(), sgVara.getVal(), sgVarb.getVal());

    delete treeZjets;

  } // end of weight loop
   
  // Calculate uncertainty of each mass bin

  float Alpha[11], Unc[11], Mzh[11] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};
  float relativeUnc[11];

  for( int im = 0; im < 11; ++im ){

    Alpha[im] = alpha[im][0];
    Unc[im] = (fabs(alpha[im][1]-alpha[im][0])>fabs(alpha[im][2]-alpha[im][0])) ? fabs(alpha[im][1]-alpha[im][0]) : fabs(alpha[im][2]-alpha[im][0]);
    relativeUnc[im] = Unc[im]/Alpha[im];

  } // end of mass points
  
  TGraphErrors *g_alpha = new TGraphErrors(11, Mzh, Alpha, 0, Unc);

  g_alpha->SetTitle("");
  g_alpha->GetXaxis()->SetTitle("");
  g_alpha->GetXaxis()->SetLabelOffset(999);
  g_alpha->GetXaxis()->SetLabelSize(0);
  g_alpha->GetXaxis()->SetLimits(800,4000);
  g_alpha->GetYaxis()->SetTitle("#alpha Ratio");  
  g_alpha->GetYaxis()->SetTitleOffset(1.3);
  g_alpha->SetMinimum(0.05);
  g_alpha->SetMaximum(50);
  g_alpha->SetLineWidth(2);
  g_alpha->SetLineColor(kBlue);
  g_alpha->SetMarkerStyle(8);
  g_alpha->SetMarkerColor(kBlue);
  g_alpha->SetFillStyle(3002);
  
  TGraph* g_unc = new TGraph(11, Mzh, relativeUnc);
  
  g_unc->SetTitle("");
  g_unc->GetXaxis()->SetTitle("m_{ZH} (GeV)");
  g_unc->GetXaxis()->SetLabelSize(0.1);
  g_unc->GetXaxis()->SetLabelOffset(0.005);
  g_unc->GetXaxis()->SetTitleSize(0.125);
  g_unc->GetXaxis()->SetTitleOffset(0.8);
  g_unc->GetXaxis()->SetLimits(800,4000);
  g_unc->GetYaxis()->SetTitle("Relative unc.");
  g_unc->GetYaxis()->SetTitleOffset(0.5);
  g_unc->GetYaxis()->SetLabelSize(0.1);
  g_unc->GetYaxis()->SetTitleSize(0.1);
  g_unc->GetYaxis()->SetNdivisions(505);
  g_unc->SetMinimum(1e-3);
  g_unc->SetMaximum(0.45);
  g_unc->SetLineWidth(2);
  g_unc->SetMarkerStyle(8);
  g_unc->SetMarkerColor(kBlack);

  TLatex lar;

  lar.SetTextSize(0.03);
  lar.SetLineWidth(5);

  float up_height = 0.8;
  float dw_height = (1-up_height)*1.375;

  TCanvas cv("cv","",0,0,1000,900);

  cv.Divide(1,2);

  TPad* cv_up = (TPad*)cv.GetListOfPrimitives()->FindObject("cv_1");
  TPad* cv_dw = (TPad*)cv.GetListOfPrimitives()->FindObject("cv_2"); 

  cv_up->SetPad(0,1-up_height,1,1);
  cv_dw->SetPad(0,0,1,dw_height);
  cv_dw->SetBottomMargin(0.25);

  cv_up->cd()->SetLogy();

  g_alpha->Draw("apz");
  g_alpha->Draw("3same");
  g_alpha->Draw("cxsame");

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.60, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s b-tag", channel.data(), catcut.data()));
  lar.DrawLatexNDC(0.15, 0.82, "b-tagging scale factor");

  cv_up->RedrawAxis();
  cv_dw->cd()->SetLogy(1);

  g_unc->Draw();

  cv.Draw();
  cv.Print(Form("alpha_bTagScale_%s_cat%s.pdf", channel.data(), catcut.data()));

}
