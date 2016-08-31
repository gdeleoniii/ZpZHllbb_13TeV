R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/PdfDiagonalizer_cc.so)
using namespace RooFit;

void test(){

  RooRealVar n("n","n",0,1000);
  n.setVal(800);
  n.setConstant(true);

  RooRealVar x("x", "x", 0, 5);
  RooRealVar lamda("lamda", "lamda", -1, -4, -0.1);
  lamda.setConstant(true);
  RooExponential model("model", "model", x, lamda);
  //  RooExtendPdf extmodel("extmodel", "extmodel", model, n);


  RooGenericPdf modelA("modelA", "modelA", "2*@0+3", x);
  RooProdPdf    modelX("modelX", "model*modelA", RooArgList(model, modelA));
  // RooExtendPdf  extmodelX("extmodelX", "extmodelX", modelX, n);

  RooFormulaVar modelB("modelB","2*x+3",x);
  RooEffProd modelY("modelY", "model*modelB", model, modelB);


  RooAbsReal* intModel = model.createIntegral(x);
  cout << model.getVal(x) << "\t" <<  model.getVal() << "\t" << intModel->getVal() << "\t" << model.getVal()/intModel->getVal() << endl;

  RooAbsReal* intModelX = modelY.createIntegral(x);
  cout << intModelX->getVal() << endl;
  /*
  TH1* htest = ext_model_sigData.createHistogram("htest", mZH, Binning(binsmZH), Extended(true));
  htest->SetMinimum(1e-4);
  htest->SetMaximum(10);
  */
  //  cout << normFactor.getVal() << "\t" << htest->Integral() <<  "\t" << (model_sigData.createIntegral(RooArgSet(mZH), Range("fullRange")))->getVal() << "\t" << ext_model_sigData.createIntegral(RooArgSet(mZH), Range("fullRange"))->getVal() << endl;

  // Plot the results on frame 

  RooPlot* fframe = x.frame();

  model.plotOn(fframe, LineColor(kGreen+1));
  //  modelA.plotOn(fframe, Normalization(1.0,RooAbsReal::NumEvent), LineColor(kRed+1));
  // modelX.plotOn(fframe,LineColor(kOrange+1) );

  //  modelB.plotOn(fframe, LineStyle(2), LineColor(kRed+1));
  //  modelY.plotOn(fframe,LineStyle(2),LineColor(kOrange+1) );

  TCanvas cv("cv","",0,0,1000,800);

  cv.cd();
  fframe->Draw();
  cv.Print("test.pdf");


}
