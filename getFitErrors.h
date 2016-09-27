// Manually calculate propagated fit uncertainties along mZH by self-difinition pdf and RooFitResults
// This is not a general function. The RooFitResult should be consistent with TF1.

vector<double> getFitErrors(TF1 f, const RooFitResult& fitRes, const RooBinning myBins){

  const RooArgList& myList = fitRes.floatParsFinal();
  const TString myFormula = f.GetExpFormula();
  const vector<string> myOrder{"nEv_sbData","a_dataSb","b_dataSb","nEv_sbSub1","a_sub1Sb","nEv_sbSub2","a_sub2Sb","nEv_sgDom","a_domSg","b_domSg","nEv_sbDom","a_domSb","b_domSb","nEv_sgSub1","a_sub1Sg","nEv_sgSub2","a_sub2Sg"};

  double cenVal[myList.getSize()], errVal[myList.getSize()];

  for( unsigned int ipar = 0; ipar < myList.getSize(); ++ipar ){

    int idx = myList.index(myOrder[ipar].data());

    cenVal[ipar] = ((RooRealVar&)myList[idx]).getVal();
    errVal[ipar] = ((RooRealVar&)myList[idx]).getError();

  }

  f.SetParameters(cenVal);

  double x = myBins.lowBound();
  vector<double> sigma;

  for( int nb = 0; nb <= myBins.numBins(); ++nb ){
    
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

      M(ipar,0) = (fabs(f.Eval(x)-f_tempUp.Eval(x)) > fabs(f.Eval(x)-f_tempDw.Eval(x))) ? fabs(f.Eval(x)-f_tempUp.Eval(x)) : fabs(f.Eval(x)-f_tempDw.Eval(x));
      Mt(0,ipar) = M(ipar,0);
      
    }

    // Very dangereous: are the elements in correlation matrix match to M?
    TMatrixD sigmaSquare = Mt*(fitRes.correlationMatrix()*M);

    sigma.push_back(TMath::Sqrt(sigmaSquare(0,0)));

    fprintf(stdout, "mZH=%i\tsigma=%.3f\n", (int)x, sigma[nb]); 

    x += myBins.binWidth(1);

  }

  return sigma;

}
