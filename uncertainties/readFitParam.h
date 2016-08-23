// This small function read fit parameters from text file.
// The keyWord is restricted by user. This function is only for special purpose.

float param(string flavor, string btag, string keyWord){

  ifstream textFile("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/uncertainties/mZhFitParam.txt");

  string textFlavor, textBtag;
  float  sbVaraMin, sbVaraMax, sbVarbMin, sbVarbMax, sgVaraMin, sgVaraMax, sgVarbMin, sgVarbMax, daVaraMin, daVaraMax, daVarbMin, daVarbMax;
  float  thisParam;

  // ignore first line of text file

  textFile.ignore(1000,'\n');

  while( textFile >> textFlavor >> textBtag >> sbVaraMin >> sbVaraMax >> sbVarbMin >> sbVarbMax >> sgVaraMin >> sgVaraMax >> sgVarbMin >> sgVarbMax >> daVaraMin >> daVaraMax >> daVarbMin >> daVarbMax ){
  
    if( flavor == textFlavor && btag == textBtag ){
      
      if     ( keyWord == "sbVaraMin" ) thisParam = sbVaraMin;
      else if( keyWord == "sbVaraMax" ) thisParam = sbVaraMax;
      else if( keyWord == "sbVarbMin" ) thisParam = sbVarbMin;
      else if( keyWord == "sbVarbMax" ) thisParam = sbVarbMax;
      else if( keyWord == "sgVaraMin" ) thisParam = sgVaraMin;
      else if( keyWord == "sgVaraMax" ) thisParam = sgVaraMax;
      else if( keyWord == "sgVarbMin" ) thisParam = sgVarbMin;
      else if( keyWord == "sgVarbMax" ) thisParam = sgVarbMax;
      else if( keyWord == "daVaraMin" ) thisParam = daVaraMin;
      else if( keyWord == "daVaraMax" ) thisParam = daVaraMax;
      else if( keyWord == "daVarbMin" ) thisParam = daVarbMin;
      else if( keyWord == "daVarbMax" ) thisParam = daVarbMax;
      else thisParam = -999;
    
    }

  }

  return thisParam;

}
