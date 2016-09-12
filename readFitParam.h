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
  float  sbVara, sbVarb, sgVara, sgVarb, daVara, daVarb;
  float  thisParam;

  // ignore first line of text file

  textFile.ignore(1000,'\n');

  while( textFile >> textFlavor >> textBtag >> sbVara >> sbVarb >> sgVara >> sgVarb >> daVara >> daVarb ){
  
    if( flavor_ == textFlavor && btag_ == textBtag ){
    
      if     ( keyWord == "sbVara"    ) thisParam = sbVara; 
      else if( keyWord == "sbVaraMin" ) thisParam = sbVara*0.5;
      else if( keyWord == "sbVaraMax" ) thisParam = sbVara*1.5;

      else if( keyWord == "sbVarb"    ) thisParam = sbVarb;
      else if( keyWord == "sbVarbMin" ) thisParam = sbVarb*0.5;
      else if( keyWord == "sbVarbMax" ) thisParam = sbVarb*1.5;

      else if( keyWord == "sgVara"    ) thisParam = sgVara;
      else if( keyWord == "sgVaraMin" ) thisParam = sgVara*0.5;
      else if( keyWord == "sgVaraMax" ) thisParam = sgVara*1.5;

      else if( keyWord == "sgVarb"    ) thisParam = sgVarb;
      else if( keyWord == "sgVarbMin" ) thisParam = sgVarb*0.5;
      else if( keyWord == "sgVarbMax" ) thisParam = sgVarb*1.5;

      else if( keyWord == "daVara"    ) thisParam = daVara;
      else if( keyWord == "daVaraMin" ) thisParam = daVara*0.5;
      else if( keyWord == "daVaraMax" ) thisParam = daVara*1.5;

      else if( keyWord == "daVarb"    ) thisParam = daVarb;
      else if( keyWord == "daVarbMin" ) thisParam = daVarb*0.5;
      else if( keyWord == "daVarbMax" ) thisParam = daVarb*1.5;

      else thisParam = -999;
    
    }

  }

  return thisParam;

}
