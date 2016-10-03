// Manually calculate propagated fit uncertainties along mZH by self-difinition pdf and RooFitResults
// This is not a general function. The RooFitResult should be consistent with TF1.
// myFormula = Form("%f*((%f*[0]*exp(-x/([1]+[2]*x))-%f*[3]*exp(-x/([4]+[5]*x))-%f*[6]*exp(-x/([7]+[8]*x)))*([12]/[9])*(%f*[9]*exp(-x/([10]+[11]*x)))/(%f*[12]*exp(-x/([13]+[14]*x)))+%f*[15]*exp(-x/([16]+[17]*x))+%f*[18]*exp(-x/([19]+[20]*x)))", normFactorVal, corr_sbData, corr_sbSub1, corr_sbSub2, corr_sgDom, corr_sbDom, corr_sgSub1, corr_sgSub2);

void getFitErrors(TF1 f, const RooFitResult& fitRes, const RooBinning myBins){

  const RooArgList& myList = fitRes.floatParsFinal();
  const TString myFormula = f.GetExpFormula();
  const vector<string> myOrder{
    "nEv_sbData",
      "a_dataSb",
      "b_dataSb",
      "nEv_sbSub1",
      "a_sub1Sb",
      "b_sub1Sb",
      "nEv_sbSub2",
      "a_sub2Sb",
      "b_sub2Sb",
      "nEv_sgDom",
      "a_domSg",
      "b_domSg",
      "nEv_sbDom",
      "a_domSb",
      "b_domSb",
      "nEv_sgSub1",
      "a_sub1Sg",
      "b_sub1Sg",
      "nEv_sgSub2",
      "a_sub2Sg",
      "b_sub2Sg"};
  
  int intIdx[myList.getSize()]; // to store RooFit internal Index, used in M & Mt below
  double cenVal[myList.getSize()], errVal[myList.getSize()];

  for( unsigned int ipar = 0; ipar < myList.getSize(); ++ipar ){

    int idx = myList.index(myOrder[ipar].data());

    intIdx[ipar] = idx;
    cenVal[ipar] = ((RooRealVar&)myList[idx]).getVal();
    errVal[ipar] = ((RooRealVar&)myList[idx]).getError();

  }

  f.SetParameters(cenVal);

  //fitRes.Print();
  //fitRes.correlationMatrix().Print("t");
  //fitRes.covarianceMatrix().Print();
  
  int numBins = myBins.numBins();
  double binWidth = myBins.binWidth(1);
  double x = myBins.lowBound();
  vector<double> sigma;

  for( int nb = 0; nb <= numBins; ++nb ){
    
    TMatrixD M (myList.getSize(),1);
    TMatrixD Mt(1,myList.getSize());

    for( unsigned int ipar = 0; ipar < myList.getSize(); ++ipar ){

      double cenTempUp[sizeof(cenVal)], cenTempDw[sizeof(cenVal)];

      memcpy(cenTempUp, cenVal, sizeof(cenVal));
      memcpy(cenTempDw, cenVal, sizeof(cenVal));      

      cenTempUp[ipar] += errVal[ipar];
      cenTempDw[ipar] -= errVal[ipar];
      
      TF1 f_tempUp("f_tempUp", myFormula.Data(), myBins.lowBound(), myBins.highBound());
      TF1 f_tempDw("f_tempDw", myFormula.Data(), myBins.lowBound(), myBins.highBound());

      f_tempUp.SetParameters(cenTempUp);
      f_tempDw.SetParameters(cenTempDw);

      M(intIdx[ipar],0) = fabs(f.Eval(x)-f_tempUp.Eval(x));
      //M(intIdx[ipar],0) = (fabs(f.Eval(x)-f_tempUp.Eval(x)) > fabs(f.Eval(x)-f_tempDw.Eval(x))) ? fabs(f.Eval(x)-f_tempUp.Eval(x)) : fabs(f.Eval(x)-f_tempDw.Eval(x));
      Mt(0,intIdx[ipar]) = M(intIdx[ipar],0);

    }

    // Very dangereous: are the elements in correlation matrix match to M?
    TMatrixD sigmaSquare = Mt*(fitRes.correlationMatrix()*M);

    sigma.push_back(TMath::Sqrt(sigmaSquare(0,0)));
    
    fprintf(stdout, "mZH=%i\tsigma^2=%.3f\tsigma=%.3f\n", (int)x, sigmaSquare(0,0), sigma[nb]); 

    x += binWidth;

  }

  double a = myBins.lowBound();
  double b = a + binWidth;

  TH1D* h_shape[3];

  for( int i = 0; i < 3; ++i ){
  
    h_shape[i] = new TH1D(Form("h_shape%i",i), "", numBins, myBins.lowBound(), myBins.highBound());

  }
  
  for( int n = 1; n <= numBins; ++n ){

    h_shape[0]->SetBinContent(n, f.Integral(a,b)/binWidth);
    h_shape[1]->SetBinContent(n, f.Integral(a,b)/binWidth + sigma[n]);
    h_shape[2]->SetBinContent(n, f.Integral(a,b)/binWidth - sigma[n]);

    a += binWidth;
    b += binWidth;
    
  }

  TFile f_shape("histo_mZH_fitParamUnc.root", "recreate");

  h_shape[0]->Write("h_mZH_fitParam_central");
  h_shape[1]->Write("h_mZH_fitParam_up");
  h_shape[2]->Write("h_mZH_fitParam_down");
  
}
