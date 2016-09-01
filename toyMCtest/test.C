R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/PdfDiagonalizer_cc.so)
using namespace RooFit;

void test(){

  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  RooRealVar n("n","n",0,1000);
  n.setVal(100);
  n.setConstant(true);

  RooBinning binsX(20, 0, 5);

  RooRealVar x("x", "x", 0, 5);
  RooRealVar lamda("lamda", "lamda", -1, -4, -0.1);

  lamda.setConstant(true);

  RooExponential model("model", "model", x, lamda);

  // Extended pdf seems only useful when perform fitting. (still need to be proved)
  // As its behaviour is same as normal pdf ( It doesnt matter what event number you put )
  // The integral of extpdf is still the same with its original pdf. This makes the plot didnt change (Pls see comments on plot part)

  RooExtendPdf extmodel("extmodel", "extmodel", model, n);

  cout << "extmodel.createIntegral(x): " << extmodel.createIntegral(x)->getVal() << endl;

  // What is exactly RooProdPdf? when multiply to functions with same observable x, what will be the output? still 2D ?

  RooGenericPdf modelA("modelA", "modelA", "2*@0+3", x);
  RooProdPdf    modelX("modelX", "model*modelA", RooArgList(model, modelA));

  // Using RooEffProd with RooFormulaVar maybe is a correct way
  // What is the difference between RooGenericPdf and RooFormulaVar? Will they have same behavior?

  RooFormulaVar modelB("modelB","2*x+3",x);
  RooEffProd modelY("modelY", "model*modelB", model, modelB);

  // modelY.Print("t");

  RooExtendPdf extmodelY("extmodelY", "extmodelY", modelY, n);

  cout << "extmodelY.createIntegral(x): " << extmodelY.createIntegral(x)->getVal() << endl;

  RooAbsReal* intModel = model.createIntegral(x);

  cout << "model.getVal(): " <<  model.getVal() << endl;  // OK so this is actually the return of exp(-2.5) = 0.082085 
  cout << "model.getVal(x): " << model.getVal(x) << endl; // this is equal to model.getVal()/model.createIntegral(x) but what is the physical meaning?
  cout << "model.createIntegral(x): " << intModel->getVal() << endl;
  cout << "model.getVal()/model.createIntegral(x): " << model.getVal()/intModel->getVal() << endl;

  // modelY = model*modelB = [ exp(-x)*(2x+3) ]
  // the modelB (RooFormulaVar) display correctly on frame.
  // output of integral of model with model.createIntegral(x) (defalut in range [0,5]) is correct: 0.993
  // output of integral of modelB with modelB.createIntegral(x) (defalut in range [0,5]) is correct: 40
  // already proved by hand

  cout << "model.createIntegral(x): " << model.createIntegral(x)->getVal() << endl;
  cout << "modelB.createIntegral(x): " << modelB.createIntegral(x)->getVal() << endl;

  // if modelY is absoutely correct, the line below is to convert the modelY to histogram (see manual: page 126)

  // seems the argument ConditionalObservables(x) must be added, BUT this will give an error (is it harmless?):
  // ERROR:InputArguments -- RooArgSet::checkForDup: ERROR argument with name modelY is already in this set

  // now the histogram is correct: if you integral exp(-x)*(2x+3) in range [0,0.5] by hand, the answer is 1.36. 
  // back to the histogram, you can see that the sum of first two bins is 1.36 

  TH1* htest = modelY.createHistogram("htest", x, Binning(binsX), ConditionalObservables(x));

  // the area of modelY is similar to the area of histogram converted from modelY
  // already proved by hand

  cout << "Integral of modelY: " << modelY.createIntegral(x)->getVal() << endl;
  cout << "Integral of hist: " << htest->Integral() << endl;

  RooPlot* fframe = x.frame();

  // If plotting the pdf on empty frame, the normalization is always 0.05 (why?)
  // But you can set normalization by using the integration of pdf (if true).
  // The additional argument RooAbsReal::Raw is used. (RooAbsReal::NumEvent failed. why?)

  model.plotOn(fframe, Normalization(model.createIntegral(x)->getVal(), RooAbsReal::Raw), LineColor(kGreen+1));
  // extmodel.plotOn(fframe, Normalization(extmodel.createIntegral(x)->getVal(), RooAbsReal::Raw), LineColor(kBlue+1));
  // modelA.plotOn(fframe, LineColor(kRed+1));
  // modelX.plotOn(fframe,LineColor(kOrange+1) );
  // modelB.plotOn(fframe, LineStyle(2), LineColor(kRed+1));
  modelY.plotOn(fframe, Normalization(modelY.createIntegral(x)->getVal(), RooAbsReal::Raw), LineStyle(2),LineColor(kOrange+1) );

  //fframe->Print("v");
  fframe->getNormVars()->Print("1");

  TCanvas cv("cv","",0,0,1000,800);

  cv.cd();
  fframe->Draw();
  cv.Print("test.pdf(");

  cv.Clear();
  cv.cd();
  htest->Draw();
  cv.Print("test.pdf)");

}
