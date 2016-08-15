R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/PdfDiagonalizer_cc.so)
using namespace RooFit;

void rooFitData(string channel, string catcut, string type, int first, int last){

  // Suppress all the INFO message

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooMsgService::instance().setSilentMode(true);

  // Input files and sum all backgrounds

  TChain* treeZjets = new TChain("tree");

  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-100to200_13TeV_toyMC.root", channel.data()));
  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-200to400_13TeV_toyMC.root", channel.data()));
  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-400to600_13TeV_toyMC.root", channel.data()));
  treeZjets->Add(Form("%s/Zjets/DYJetsToLL_M-50_HT-600toInf_13TeV_toyMC.root", channel.data()));

  // Define all the variables from the trees 

  RooRealVar cat("cat", "", 0, 2);
  RooRealVar mJet("prmass", "",  30.,  300.);
  RooRealVar mZH("mllbb", "M_{ZH}", 800., 4000., "GeV");
  RooBinning mZHbin(32, 800., 4000.);

  mZH.setRange("fullRange", 800., 4000.);

  TCut catCut = Form("cat==%s", catcut.c_str());
  TCut sbCut  = "prmass>30 && !(prmass>65 && prmass<135) && prmass<300";
  TCut sigCut = "prmass>105 && prmass<135";

  int N = last-first;
  int iw = N-1;
  float alphaScale[11][N], alphaCentral[11];

  TF1* f_alpha = new TF1("f_alpha", "[0]*TMath::Exp([1]*x+[2]/x)/TMath::Exp([3]*x+[4]/x)", 800, 4000);
     
  for( int nw = last; nw >= first; --nw ){
    
    RooRealVar evWeight(Form("evweight%02i",nw), "", -1.e10, 1.e10);
    RooArgSet variables(cat, mJet, mZH, evWeight);

    // Create a dataset from a tree -> to process unbinned likelihood fitting

    RooDataSet dataSetZjetsSB("dataSetZjetsSB", "dataSetZjetsSB", variables, Cut(catCut && sbCut),  WeightVar(evWeight), Import(*treeZjets));  
    RooDataSet dataSetZjetsSG("dataSetZjetsSG", "dataSetZjetsSG", variables, Cut(catCut && sigCut), WeightVar(evWeight), Import(*treeZjets));
  
    // Total event numbers

    RooRealVar nSBMcEvents("nSBMcEvents", "nSBMcEvents", 0., 1.e10);
    RooRealVar nSGMcEvents("nSGMcEvents", "nSGMcEvents", 0., 1.e10);

    nSBMcEvents.setVal(dataSetZjetsSB.sumEntries());
    nSBMcEvents.setConstant(true);
  
    nSGMcEvents.setVal(dataSetZjetsSG.sumEntries());
    nSGMcEvents.setConstant(true);
  
    // Alpha ratio part

    // Fit ZH mass in side band 

    float bmin, bmax;

    if( channel == "ele" ){
      bmin = (catcut=="1") ?  600. : 1400.;
      bmax = (catcut=="1") ? 1500. : 2600.;
    }

    else if( channel == "mu" ){
      bmin = (catcut=="1") ? 100. : 2100.; 
      bmax = (catcut=="1") ? 900. : 2900.;
    }
    
    RooRealVar a("a", "a", -0.002, -0.005, 0.);
    RooRealVar b("b", "b", (bmin+bmax)*0.5, bmin, bmax);

    RooGenericPdf model_ZHSB("model_ZHSB", "model_ZHSB", "TMath::Exp(@1*@0+@2/@0)", RooArgSet(mZH,a,b));
    RooExtendPdf ext_model_ZHSB("ext_model_ZHSB", "ext_model_ZHSB", model_ZHSB, nSBMcEvents);

    ext_model_ZHSB.fitTo(dataSetZjetsSB, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));

    float p0 = a.getVal();
    float p1 = b.getVal();

    // Fit ZH mass in signal region

    float dmin, dmax;
    
    if( channel == "ele" ){
      dmin = (catcut=="1") ? 3100. : 0.;
      dmax = (catcut=="1") ? 3900. : 1.;
    }
    
    else if( channel == "mu" ){
      dmin = (catcut=="1") ? 0. : 8.;
      dmax = (catcut=="1") ? 1. : 18.;
    }
    
    RooRealVar c("c", "c", -0.002, -0.005, 0.);
    RooRealVar d("d", "d", (dmin+dmax)*0.5, dmin, dmax);

    RooGenericPdf model_ZHSG("model_ZHSG", "model_ZHSG", "TMath::Exp(@1*@0+@2/@0)", RooArgSet(mZH,c,d));
    RooExtendPdf ext_model_ZHSG("ext_model_ZHSG", "ext_model_ZHSG", model_ZHSG, nSGMcEvents);
    
    ext_model_ZHSG.fitTo(dataSetZjetsSG, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));

    float p2 = c.getVal();
    float p3 = d.getVal();
 
    // Set the model of alpha ratio

    float normConst = ((TF1*)ext_model_ZHSB.asTF(mZH,RooArgList(a,b)))->Integral(800,4000) / ((TF1*)ext_model_ZHSG.asTF(mZH,RooArgList(c,d)))->Integral(800,4000);

    f_alpha->SetParameters(normConst,p2,p3,p0,p1);

    int mzh = 800;
    for( int im = 0; im < 11; ++im ){

      if( nw != first )
	alphaScale[im][iw] = f_alpha->Eval(mzh);
      else
	alphaCentral[im] = f_alpha->Eval(mzh);

      mzh += (mzh<2000) ? 200 : 500;

    }

    --iw;

    fprintf(stdout, "weight=%i\t(sb)p0=%f\tp1=%f\t(sg)p2=%f\tp3=%f\n", nw, p0, p1, p2, p3);

  } // end of weight loop

  // Calculate RMS value of each mass bin

  float Unc[11], Mzh[11] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000};
  float relativeUnc[11];

  for( int im = 0; im < 11; ++im ){

    if( type == "mur1" )
      Unc[im] = ( fabs(alphaScale[im][1]-alphaCentral[im]) > fabs(alphaScale[im][0]-alphaCentral[im]) ) ?
	fabs(alphaScale[im][1]-alphaCentral[im]) : fabs(alphaScale[im][0]-alphaCentral[im]);

    else 
      Unc[im] = TMath::RMS(N, alphaScale[im]);

    relativeUnc[im] = Unc[im]/alphaCentral[im];

  } // end of mass points

  TGraphErrors *g_alpha = new TGraphErrors(11, Mzh, alphaCentral, 0, Unc);

  g_alpha->SetTitle("");
  g_alpha->GetXaxis()->SetTitle("");
  g_alpha->GetXaxis()->SetLabelOffset(999);
  g_alpha->GetXaxis()->SetLabelSize(0);
  g_alpha->GetXaxis()->SetLimits(800,4000);
  g_alpha->GetYaxis()->SetTitle("#alpha Ratio");  
  g_alpha->GetYaxis()->SetTitleOffset(1.3);
  g_alpha->SetMinimum(0);
  g_alpha->SetMaximum(1.8);
  g_alpha->SetLineWidth(2);
  g_alpha->SetLineColor(kBlue);
  g_alpha->SetMarkerStyle(8);
  g_alpha->SetMarkerColor(kBlue);
  g_alpha->SetFillStyle(1001);
  g_alpha->SetFillColor(kYellow);  
  
  TGraph* g_unc = new TGraph(11, Mzh, relativeUnc);
  
  g_unc->SetTitle("");
  g_unc->GetXaxis()->SetTitle("m_{ZH} (GeV)");
  g_unc->GetXaxis()->SetLabelSize(0.1);
  g_unc->GetXaxis()->SetLabelOffset(0.005);
  g_unc->GetXaxis()->SetTitleSize(0.125);
  g_unc->GetXaxis()->SetTitleOffset(0.8);
  g_unc->GetXaxis()->SetLimits(800,4000);
  g_unc->GetYaxis()->SetTitle("Relative unc.");
  g_unc->GetYaxis()->SetTitleOffset(0.45);
  g_unc->GetYaxis()->SetLabelSize(0.1);
  g_unc->GetYaxis()->SetTitleSize(0.1);
  g_unc->GetYaxis()->SetNdivisions(505);
  g_unc->SetMinimum(0);
  g_unc->SetMaximum(0.45);
  g_unc->SetLineWidth(2);
  g_unc->SetMarkerStyle(8);
  g_unc->SetMarkerColor(kBlack);

  TLegend* leg = new TLegend(0.15,0.15,0.35,0.25);

  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);

  TLatex* lar = new TLatex();
  
  lar->SetTextSize(0.035);
  lar->SetLineWidth(5);

  float up_height     = 0.8;
  float dw_correction = 1.375;
  float dw_height     = (1-up_height)*dw_correction;

  TCanvas cv("cv","",0,0,1000,900);
  cv.Divide(1,2);

  TPad* cv_up = (TPad*)cv.GetListOfPrimitives()->FindObject("cv_1");
  TPad* cv_dw = (TPad*)cv.GetListOfPrimitives()->FindObject("cv_2"); 

  cv_up->SetPad(0,1-up_height,1,1);
  cv_dw->SetPad(0,0,1,dw_height);
  cv_dw->SetBottomMargin(0.25);

  cv_up->cd();

  g_alpha->Draw("apz");
  g_alpha->Draw("3same");
  g_alpha->Draw("cxsame");
  leg->AddEntry(g_alpha, "alpha ratio with uncertainties", "lf");
  leg->Draw();
  lar->DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{2015}}");
  lar->DrawLatexNDC(0.55, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
  lar->DrawLatexNDC(0.72, 0.80, Form("%s  %s btag", channel.data(), catcut.data()));

  cv_up->RedrawAxis();
  cv_dw->cd();

  g_unc->Draw();

  cv.Draw();
  cv.Print(Form("alpha_%sScale_%s_cat%s.pdf", type.data(), channel.data(), catcut.data()));  

}
