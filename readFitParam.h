// This small class is used to read fit parameters from text file.
// The keyWord is restricted by user. This function is only for special purpose.

class param{

 public:

  param(string, string);
  float value(string);

 private:

  string flavor_;
  string btag_;

};

param::param(string flavor, string btag){

  flavor_ = flavor;
  btag_   = btag;  

}

float param::value(string keyWord){

  ifstream textFile("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/mZhFitParam.txt");

  string textFlavor, textBtag;
  float  a_domSb, b_domSb, a_domSg, b_domSg, a_subSb, b_subSb, a_subSg, b_subSg, a_datSb, b_datSb, lamda;
  float  thisParam;

  // ignore first line of text file

  textFile.ignore(1000,'\n');

  while( textFile >> textFlavor >> textBtag >> a_domSb >> b_domSb >> a_domSg >> b_domSg >> a_subSb >> b_subSb >> a_subSg >> b_subSg >> a_datSb >> b_datSb >> lamda ){
  
    if( flavor_ == textFlavor && btag_ == textBtag ){
    
      if     ( keyWord == "a_domSb"    ) thisParam = a_domSb; 
      else if( keyWord == "a_domSbMin" ) thisParam = a_domSb*0.5;
      else if( keyWord == "a_domSbMax" ) thisParam = a_domSb*1.5;

      else if( keyWord == "b_domSb"    ) thisParam = b_domSb;
      else if( keyWord == "b_domSbMin" ) thisParam = b_domSb*0.5;
      else if( keyWord == "b_domSbMax" ) thisParam = b_domSb*1.5;

      else if( keyWord == "a_domSg"    ) thisParam = a_domSg;
      else if( keyWord == "a_domSgMin" ) thisParam = a_domSg*0.5;
      else if( keyWord == "a_domSgMax" ) thisParam = a_domSg*1.5;

      else if( keyWord == "b_domSg"    ) thisParam = b_domSg;
      else if( keyWord == "b_domSgMin" ) thisParam = b_domSg*0.5;
      else if( keyWord == "b_domSgMax" ) thisParam = b_domSg*1.5;

      else if( keyWord == "a_subSb"    ) thisParam = a_subSb;
      else if( keyWord == "a_subSbMin" ) thisParam = a_subSb*0.5;
      else if( keyWord == "a_subSbMax" ) thisParam = a_subSb*1.5;

      else if( keyWord == "b_subSb"    ) thisParam = b_subSb;
      else if( keyWord == "b_subSbMin" ) thisParam = b_subSb*0.5;
      else if( keyWord == "b_subSbMax" ) thisParam = b_subSb*1.5;

      else if( keyWord == "a_subSg"    ) thisParam = a_subSg;
      else if( keyWord == "a_subSgMin" ) thisParam = a_subSg*0.5;
      else if( keyWord == "a_subSgMax" ) thisParam = a_subSg*1.5;

      else if( keyWord == "b_subSg"    ) thisParam = b_subSg;
      else if( keyWord == "b_subSgMin" ) thisParam = b_subSg*0.5;
      else if( keyWord == "b_subSgMax" ) thisParam = b_subSg*1.5;

      else if( keyWord == "a_datSb"    ) thisParam = a_datSb;
      else if( keyWord == "a_datSbMin" ) thisParam = a_datSb*0.5;
      else if( keyWord == "a_datSbMax" ) thisParam = a_datSb*1.5;

      else if( keyWord == "b_datSb"    ) thisParam = b_datSb;
      else if( keyWord == "b_datSbMin" ) thisParam = b_datSb*0.5;
      else if( keyWord == "b_datSbMax" ) thisParam = b_datSb*1.5;

      else if( keyWord == "lamda"    ) thisParam = lamda;
      else if( keyWord == "lamdaMin" ) thisParam = lamda*1.5;
      else if( keyWord == "lamdaMax" ) thisParam = lamda*0.5;

      else thisParam = -999;
    
    }

  }

  return thisParam;

}
