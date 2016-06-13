/**
   \class    standalone_LumiReWeighting standalone_LumiReWeighting.h "PhysicsTools/Utilities/interface/standalone_LumiReWeighting.h"
   \brief    Class to provide lumi weighting for analyzers to weight "flat-to-N" MC samples to data

   This class will trivially take two histograms:
   1. The generated "flat-to-N" distributions from a given processing (or any other generated input)
   2. A histogram generated from the "estimatePileup" macro here:

   https://twiki.cern.ch/twiki/bin/view/CMS/LumiCalc#How_to_use_script_estimatePileup

   and produce weights to convert the input distribution (1) to the latter (2).

   \author Shin-Shan Eiko Yu, Salvatore Rappoccio, modified by Mike Hildreth
  
*/
#ifndef standalone_LumiReWeighting_cxx
#define standalone_LumiReWeighting_cxx
#include "TH1.h"
#include "TFile.h"
#include <string>
#include "standalone_LumiReWeighting.h"

//=======================================================
// For 2015 Data and MC
//=======================================================

double MC2015[52]={
  4.8551E-07,
  1.74806E-06,
  3.30868E-06,
  1.62972E-05,
  4.95667E-05,
  0.000606966,
  0.003307249,
  0.010340741,
  0.022852296,
  0.041948781,
  0.058609363,
  0.067475755,
  0.072817826,
  0.075931405,
  0.076782504,
  0.076202319,
  0.074502547,
  0.072355135,
  0.069642102,
  0.064920999,
  0.05725576,
  0.047289348,
  0.036528446,
  0.026376131,
  0.017806872,
  0.011249422,
  0.006643385,
  0.003662904,
  0.001899681,
  0.00095614,
  0.00050028,
  0.000297353,
  0.000208717,
  0.000165856,
  0.000139974,
  0.000120481,
  0.000103826,
  8.88868E-05,
  7.53323E-05,
  6.30863E-05,
  5.21356E-05,
  4.24754E-05,
  3.40876E-05,
  2.69282E-05,
  2.09267E-05,
  1.5989E-05,
  4.8551E-06,
  2.42755E-06,
  4.8551E-07,
  2.42755E-07,
  1.21378E-07,
  4.8551E-08
};

double Data2015[52]={
  147883,
  666514,
  886596,
  1.31858e+06,
  2.24533e+06,
  5.14783e+06,
  1.6764e+07, 
  6.81266e+07,
  2.00024e+08,
  3.57504e+08,
  4.51868e+08,
  4.6053e+08,
  3.93032e+08,
  2.75298e+08,
  1.56242e+08,
  7.26553e+07,
  2.93184e+07,
  1.19483e+07,
  5.79335e+06,
  3.04766e+06,
  1.40506e+06,
  512252, 
  146792, 
  35303.5,
  8271.34,
  2235.61,
  721.381,
  258.85, 
  97.2709,
  36.8716,
  13.7275,
  4.93171,
  1.69241,
  0.551894,
  0.170601, 
  0.0499358,
  0.0138339,
  0.00362662, 
  0.000899606,
  0.000211148,
  4.68928e-05,
  9.85392e-06,
  1.95929e-06,
  3.6862e-07, 
  6.56248e-08,
  1.10534e-08,
  1.76248e-09,
  2.61497e-10,
  4.768e-11,
  0, 
  0, 
  0
};

double Data2015Up[52]={
  124589,
  611524,
  813990,
  1.17492e+06,
  1.86543e+06,
  3.78652e+06,
  1.06091e+07,
  3.92702e+07,
  1.31164e+08,
  2.76072e+08,
  3.9622e+08, 
  4.43974e+08,
  4.20871e+08,
  3.3768e+08, 
  2.24975e+08,
  1.23823e+08,
  5.7309e+07, 
  2.38058e+07,
  1.0299e+07, 
  5.26091e+06,
  2.85236e+06,
  1.36371e+06,
  526643,
  163061,
  42592.3,
  10584.7,
  2923.75,
  955.685,
  351.3, 
  136.772, 
  54.3136, 
  21.4253, 
  8.23874, 
  3.05183, 
  1.08184, 
  0.36577, 
  0.117761,
  0.0360763, 
  0.0105133, 
  0.00291401, 
  0.000768174,
  0.000192592,
  4.59223e-05,
  1.04141e-05,
  2.2461e-06, 
  4.60737e-07,
  8.98872e-08,
  1.66763e-08,
  2.94239e-09,
  4.96039e-10,
  9.24658e-11,
  1.00203e-12
};

double Data2015Down[52]={
  173216,
  727332,
  978926,
  1.50593e+06,
  2.79307e+06,
  7.41815e+06,
  2.86975e+07,
  1.17312e+08,
  2.87819e+08,
  4.37121e+08,
  4.90852e+08,
  4.50031e+08,
  3.35001e+08,
  1.98139e+08,
  9.33848e+07,
  3.68086e+07,
  1.40969e+07,
  6.43364e+06,
  3.27039e+06,
  1.44846e+06,
  495066, 
  130081, 
  28609.8,
  6333.37,
  1680.79,
  533.046,
  185.383,
  66.6646,
  23.8662,
  8.28654,
  12.74673,
  0.862062,
  0.255191,
  0.0711318,
  0.0186567,
  0.00460324,
  0.00106833,
  0.000233212,
  4.78844e-05,
  9.24779e-06,
  1.67991e-06,
  2.8704e-07,
  4.6131e-08,
  6.97509e-09,
  9.94095e-10,
  1.33153e-10,
  1.59604e-11,
  0,
  0,
  0,
  0,
  0
};

standalone_LumiReWeighting::standalone_LumiReWeighting(int mode) {

  std::cout << "=======================================================================" << std::endl;
  
  std::vector<double> MC_distr;
  std::vector<double> Lumi_distr;

  MC_distr.clear();
  Lumi_distr.clear();
  switch (mode)
    {
    case 0:
      std::cout << "Using central value " << std::endl;
      break;
    case 1:
      std::cout << "Using +1 sigma 5% value " << std::endl;
      break;
    case -1:
      std::cout << "Using -1 sigma 5% value " << std::endl;
      break;
    default:
      std::cout << "Using central value " << std::endl;
      break;
    } // end of switch

  Int_t NBins = 52;
  
  for( int i=0; i< NBins; ++i) {
    switch (mode){
    case 0:
      Lumi_distr.push_back(Data2015[i]);
      break;
    case 1:
      Lumi_distr.push_back(Data2015Up[i]);
      break;
    case -1:
      Lumi_distr.push_back(Data2015Down[i]);
      break;
    default:
      Lumi_distr.push_back(Data2015[i]);
      break;
    } // end of switch

    MC_distr.push_back(MC2015[i]);
  } // end of loop over bins

  // no histograms for input: use vectors  
  // now, make histograms out of them:
  // first, check they are the same size...

  if( MC_distr.size() != Lumi_distr.size() ){   
    std::cout << "MC_distr.size() = " << MC_distr.size() << std::endl;
    std::cout << "Lumi_distr.size() = " << Lumi_distr.size() << std::endl;
    std::cerr <<"ERROR: standalone_LumiReWeighting: input vectors have different sizes. Quitting... \n";

  }

  weights_ = new TH1D(Form("luminumer_%d",mode),
 		      Form("luminumer_%d",mode),
 		      NBins,0.0, double(NBins));

  TH1D* den = new TH1D(Form("lumidenom_%d",mode),
 		       Form("lumidenom_%d",mode),
 		       NBins,0.0, double(NBins));
  
  for(int ibin = 1; ibin<NBins+1; ++ibin ) {
    weights_->SetBinContent(ibin, Lumi_distr[ibin-1]);
    den->SetBinContent(ibin,MC_distr[ibin-1]);
  }

  std::cout << "Data Input " << std::endl;
  for(int ibin = 1; ibin<NBins+1; ++ibin){
    std::cout << "   " << ibin-1 << " " << weights_->GetBinContent(ibin) << std::endl;
  }
  std::cout << "MC Input " << std::endl;
  for(int ibin = 1; ibin<NBins+1; ++ibin){
    std::cout << "   " << ibin-1 << " " << den->GetBinContent(ibin) << std::endl;
  }

  // check integrals, make sure things are normalized

  double deltaH = weights_->Integral();
  if(fabs(1.0 - deltaH) > 0.02 ) { //*OOPS*...
    weights_->Scale( 1.0/ weights_->Integral() );
  }

  double deltaMC = den->Integral();
  if(fabs(1.0 - deltaMC) > 0.02 ) {
    den->Scale(1.0/ den->Integral());
  }

  weights_->Divide( den );  // so now the average weight should be 1.0    

  std::cout << "Reweighting: Computed Weights per In-Time Nint " << std::endl;

  for(int ibin = 1; ibin<NBins+1; ++ibin){
    std::cout << "   " << ibin-1 << " " << weights_->GetBinContent(ibin) << std::endl;
  }

  std::cout << "=======================================================================" << std::endl;

}

standalone_LumiReWeighting::~standalone_LumiReWeighting(){}

double standalone_LumiReWeighting::weight( double npv ){
  int bin = weights_->GetXaxis()->FindBin( npv );
  return weights_->GetBinContent( bin );
}

#endif
