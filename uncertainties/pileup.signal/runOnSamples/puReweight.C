#include <iostream>
#include <cstdio>
#include <standalone_LumiReWeighting.cc>

void puReweight(){

  standalone_LumiReWeighting LumiWeights_central(0);
  standalone_LumiReWeighting LumiWeights_scaleup(1);
  standalone_LumiReWeighting LumiWeights_scaledw(-1);
  
  FILE* fout = fopen("LumiWeights.txt", "w");

  for(int ntrue = 0; ntrue < 52; ++ntrue){
 
    fprintf(fout, "%d\t%g\t%g\t%g\n",  
	    ntrue,
	    LumiWeights_central.weight(ntrue),
	    LumiWeights_scaleup.weight(ntrue),
	    LumiWeights_scaledw.weight(ntrue));
    
  }

  fclose(fout);

}
