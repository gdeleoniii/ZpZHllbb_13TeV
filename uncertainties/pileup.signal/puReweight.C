#include <iostream>
#include <fstream>
#include <standalone_LumiReWeighting.cc>

void puReweight(){

  standalone_LumiReWeighting LumiWeights_central(0);
  standalone_LumiReWeighting LumiWeights_scaleup(1);
  standalone_LumiReWeighting LumiWeights_scaledw(-1);
  
  fstream file;

  file.open("LumiWeights.txt", ios::out);

  for(int ntrue = 0; ntrue < 54; ++ntrue){
 
    file << ntrue << "\t\t" 
	 << LumiWeights_central.weight(ntrue) << "\t\t" 
	 << LumiWeights_scaleup.weight(ntrue) << "\t\t" 
	 << LumiWeights_scaledw.weight(ntrue) << endl;
    
  }

  file.close();

}
