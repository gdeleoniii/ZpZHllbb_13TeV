R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/PdfDiagonalizer_cc.so)
using namespace RooFit;

void signalShape(string channel, string catcut){

  // Suppress all the INFO message

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooMsgService::instance().setSilentMode(true);

  int mass[10] = {800,1000,1200,1400,1600,1800,2000,2500,3000,3500};

  for( int i = 0; i < 10; ++i ){

    // Input files

    TChain* treeSignal = new TChain("tree");

    treeSignal->Add(Form("signal/ZprimeToZhToZlephbb_M-%i_13TeV_%sMiniTree.root", mass[i], channel.data()));

    // Define all the variables from the trees

    RooRealVar cat ("cat", "", 0, 2);
    RooRealVar mZH ("mllbb", "M_{ZH}", 600., 4500., "GeV");
    RooRealVar evWeight("evweight", "", 0., 1.e3);

    mZH.setRange("fullRange", 600., 4500.);

    RooBinning binsmZH(78, 600, 4500);

    RooArgSet variables(cat, mZH, evWeight);

    TCut catCut = Form("cat==%s", catcut.c_str());

    // Create a dataset from a tree -> to process an unbinned likelihood fitting

    RooDataSet dataSetSignal("dataSetSignal", "dataSetSignal", variables, Cut(catCut), WeightVar(evWeight), Import(*treeSignal));

    RooRealVar nEvents("nEvents", "nEvents", 0., 1.e10);

    nEvents.setVal(dataSetSignal.sumEntries());
    nEvents.setConstant(true);

    // Fit the signal shape -> using Crystal Ball function

    RooRealVar m("m", "mean", 0., 1.e5);
    RooRealVar s("s", "sigma", 25., 80.);
    RooRealVar a("a", "alpha", 1., 0.1, 3.);
    RooRealVar n("n", "n", 80., 0., 160.);

    m.setVal(mass[i]);
    m.setConstant(true);

    RooCBShape model("model", "Cystal Ball Function", mZH, m, s, a, n);
    RooExtendPdf ext_model("ext_model", "Extended Cystal Ball Function", model, nEvents);

    ext_model.fitTo(dataSetSignal, SumW2Error(true), Extended(true), Range("fullRange"), Strategy(2), Minimizer("Minuit2"), Save(1));

    fprintf(stdout, "mass=%i\tmean=%f\tsigma=%f\talpha=%f\tn=%f\n", mass[i], m.getVal(), s.getVal(), a.getVal(), n.getVal());

    // Plot the results on frame

    RooPlot* signalFrame = mZH.frame();
  
    dataSetSignal.plotOn(signalFrame, Binning(binsmZH));
    ext_model.plotOn(signalFrame, LineColor(kBlue));

    RooPlot* signalPullFrame = mZH.frame();

    // Output the results

    TLatex lar;
    lar.SetTextSize(0.03);
    lar.SetLineWidth(5);

    float up_height = 0.82;
    float dw_height = (1-up_height)*1.445;

    TCanvas c("c","",0,0,1000,800);
  
    c.Divide(1,2);

    TPad* c_up = (TPad*)c.GetListOfPrimitives()->FindObject("c_1");
    TPad* c_dw = (TPad*)c.GetListOfPrimitives()->FindObject("c_2"); 

    c_up->SetPad(0,1-up_height,1,1);
    c_dw->SetPad(0,0,1,dw_height);
    c_dw->SetBottomMargin(0.25);
    c_up->cd();

    signalFrame->SetTitle("");
    signalFrame->SetMinimum(0);
    signalFrame->GetXaxis()->SetTitle("");
    signalFrame->GetXaxis()->SetLabelOffset(999);
    signalFrame->Draw();

    lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
    lar.DrawLatexNDC(0.65, 0.92, "L = 2.512 fb^{-1} at #sqrt{s} = 13 TeV");
    lar.DrawLatexNDC(0.15, 0.86, Form("%s, %s b-tag", channel.data(), catcut.data()));
    lar.DrawLatexNDC(0.15, 0.82, "Crystal Ball function");

    c_up->RedrawAxis();

    c_dw->cd();

    signalPullFrame->addObject(signalFrame->pullHist(), "P");
    signalPullFrame->SetTitle("");
    signalPullFrame->GetYaxis()->SetTitle("Pulls");
    signalPullFrame->GetYaxis()->SetTitleOffset(0.25);
    signalPullFrame->GetXaxis()->SetLabelSize(0.125);
    signalPullFrame->GetXaxis()->SetTitleSize(0.125);
    signalPullFrame->GetYaxis()->SetLabelSize(0.125);
    signalPullFrame->GetYaxis()->SetTitleSize(0.125);
    signalPullFrame->GetYaxis()->SetNdivisions(505);
    signalPullFrame->SetMinimum(-4);
    signalPullFrame->SetMaximum(4);
    signalPullFrame->Draw();

    c.Draw();
    c.Print(Form("signalShape_%s_%sbtag.pdf%s", channel.data(), catcut.data(), (i==0?"(":(i==9?")":""))));

    delete treeSignal;

  }

}
