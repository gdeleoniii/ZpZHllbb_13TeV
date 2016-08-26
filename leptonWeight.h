// This function calculate and return the efficiency/trigger scale factor of each lepton.
// Input arguments: 
//    h2_lep   -> the 2d scale factor histogram provided by POG. 
//    thisLep  -> 4-vector of this lepton.
//    absEta   -> take the absolute of eta in default (false for electron)
//    lepScale -> scale the weight to [central=0(default), up=1, down=-1]

float leptonWeight(TH2F* h2_lep, TLorentzVector* thisLep=NULL, bool absEta=true, int lepScale=0){

  float thisLepPt, thisLepEta;

  float minPt = h2_lep->GetYaxis()->GetBinLowEdge(h2_lep->GetYaxis()->GetFirst());
  float maxPt = h2_lep->GetYaxis()->GetBinUpEdge(h2_lep->GetYaxis()->GetLast());

  if( thisLep->Pt() < minPt ) 
    thisLepPt = minPt + 0.001;

  else if( thisLep->Pt() > maxPt )
    thisLepPt = maxPt - 0.001;

  else
    thisLepPt = thisLep->Pt();

  thisLepEta = (absEta) ? fabs(thisLep->Eta()) : thisLep->Eta();

  int thisPtBin  = h2_lep->GetYaxis()->FindBin(thisLepPt);
  int thisEtaBin = h2_lep->GetXaxis()->FindBin(thisLepEta);

  float thisWeight = h2_lep->GetBinContent(thisEtaBin, thisPtBin) + ( lepScale * h2_lep->GetBinError(thisEtaBin, thisPtBin));

  return thisWeight;
  
}
