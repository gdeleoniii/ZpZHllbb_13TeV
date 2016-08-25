float leptonWeight(TH2F* h2_lep, TLorentzVector* thisLep=NULL, int lepScale=0){

  float thisLepPt;

  float minPt = h2_lep->GetYaxis()->GetBinLowEdge(h2_lep->GetYaxis()->GetFirst());
  float maxPt = h2_lep->GetYaxis()->GetBinUpEdge(h2_lep->GetYaxis()->GetLast());

  if( thisLep->Pt() < minPt ) 
    thisLepPt = minPt + 0.001;

  else if( thisLep->Pt() > maxPt )
    thisLepPt = maxPt - 0.001;

  else
    thisLepPt = thisLep->Pt();

  int thisPtBin  = h2_lep->GetYaxis()->FindBin(thisLepPt);
  int thisEtaBin = h2_lep->GetXaxis()->FindBin(thisLep->Eta());

  return h2_lep->GetBinContent(thisEtaBin, thisPtBin) + ( lepScale * h2_lep->GetBinError(thisEtaBin, thisPtBin));
  
}
