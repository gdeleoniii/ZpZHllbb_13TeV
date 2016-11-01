void myEvNum(string channel, int btag){

  double mzh[7] = {800,1000,1200,1400,1600,1800,2000};

  TFile* f_data = TFile::Open(Form("sourceRoot/mZH%sbtag%i.root",channel.data(),btag));
  TFile* f_bkgs = TFile::Open(Form("sourceRoot/background_FitDev_%s_cat%i.root",channel.data(),btag));

  TH1D* h_data = (TH1D*)(f_data->Get("data_obs"));
  TH1D* h_bkgs = (TH1D*)(f_bkgs->Get("background_FitDev"));

  FILE* f_out = fopen(Form("%s_%ibtag_nEvents.txt",channel.data(),btag), "w");

  fprintf(f_out, "mass(>)\tN_Data\tDataStatUnc\tN_Zjets\n");

  for( int i = 0; i < 7; ++i ){

    float ev_data = h_data->Integral(h_data->FindBin(mzh[i]), h_data->GetNbinsX());
    float ev_bkgs = h_bkgs->Integral(h_bkgs->FindBin(mzh[i]), h_bkgs->GetNbinsX());

    fprintf(f_out, "%i\t%.3f\t%.3f\t%.3f\n", (int)mzh[i], ev_data, sqrt(ev_data), ev_bkgs);
      
  }

  fclose(f_out);

}

void eventNumbers(){

  myEvNum("ele",1);
  myEvNum("ele",2);
  myEvNum("mu",1);
  myEvNum("mu",2);

}
