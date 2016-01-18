#include <iostream>
#include <string>
#include <vector>
#include "untuplizer.h"

bool TriggerStatus(TreeReader &data, std::string TRIGGERNAME){

  const std::string* trigName    = data.GetPtrString("hlt_trigName");
  const vector<bool>& trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));

  bool triggerStat = false;

  for(int it = 0; it < data.GetPtrStringSize(); it++){

    std::string thisTrig = trigName[it];
    bool results = trigResult[it];

    if( thisTrig.find(TRIGGERNAME) != std::string::npos && results ){

      triggerStat = true;
      break;

    }

  }

  return triggerStat;

}

bool FilterStatus(TreeReader &data, std::string FILTERNAME){

  const std::string* filterName    = data.GetPtrString("hlt_filterName");
  const vector<bool>& filterResult = *((vector<bool>*) data.GetPtr("hlt_filterResult"));

  bool filterStat = false;
    
  for(int it = 0; it < data.GetPtrStringSize(); it++){
    
    std::string thisFilter = filterName[it];
    bool results = filterResult[it];
      
    if( thisFilter.find(FILTERNAME) != std::string::npos && results ){

      filterStat = true;
      break;

    }
      
  }

  return filterStat;

}
