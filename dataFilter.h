#include <iostream>
#include <string>
#include <vector>
#include "untuplizer.h"

bool FilterStatus(TreeReader &data, std::string TRIGNAME){

  std::string* trigName    = data.GetPtrString("hlt_trigName");
  vector<bool>& trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));

  bool triggerStat = false;
    
  for(int it = 0; it < data.GetPtrStringSize(); it++){
    
    std::string thisTrig = trigName[it];
    bool results = trigResult[it];
      
    if( thisTrig.find(TRIGNAME) != std::string::npos && results ){

      triggerStat = true;
      break;

    }
      
  }

  return triggerStat;

}
