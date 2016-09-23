TH1D* createHistogram(TF1* f, double xMin, double xMax, int nBins){

  double binWidth = (xMax-xMin)/nBins;
  
  TH1D* h = new TH1D("h", "", nBins, xMin, xMax);

  double a = xMin;
  double b = xMin+binWidth;
  
  for( int n = 1; n <= nBins; ++n ){

    h->SetBinContent(n, f->Integral(a, b)/binWidth);

    a += binWidth;
    b += binWidth;
    
  }

  return h;

}
