float leptonWeight(TH2F* h2_lep, TLorentzVector* thisLep=NULL, int lepScale=0){

  float thisWeight = 1;
  float thisLepPt;

  float minPt = h2_lep->GetYaxis()->GetBinLowEdge(h2_lep->GetYaxis()->GetFirst());
  float maxPt = h2_lep->GetYaxis()->GetBinUpEdge(h2_lep->GetYaxis()->GetLast());

  if( thisLep->Pt() < minPt ) 
    thisLepPt = minPt + 0.001;

  else if( thisLep->Pt() > maxPt )
    thisLepPt = maxPt - 0.001;

  else
    thisLepPt = thisLep->Pt();

  int thisEtaBin = h2_lep->GetXaxis()->FindBin(thisLep->Eta());
  int thisPtBin  = h2_lep->GetYaxis()->FindBin(thisLepPt);

  thisWeight  = h2_lep->GetBinContent(thisEtaBin, thisPtBin);
  thisWeight += (lepScale>0) ? lepScale*h2_lep->GetBinErrorUp(thisEtaBin, thisPtBin) : lepScale*h2_lep->GetBinErrorLow(thisEtaBin, thisPtBin);
     
  return thisWeight;

}
