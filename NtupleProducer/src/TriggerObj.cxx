#include "include/NtupleProducer.h"
#include <iostream>
#include <iomanip>

ClassImp(TriggerObj)

TriggerObj::TriggerObj()
{
}

TriggerObj::~TriggerObj()
{
}

void TriggerObj::read()
{
   
   // 
   _ID = idx;
   
   _pT  = ntP->triggerobject_pt->at(_ID);
   _eta = ntP->triggerobject_eta->at(_ID); 
   _phi = ntP->triggerobject_phi->at(_ID);
   
   _collection = ntP->triggerobject_collection->at(_ID);
   
   //std::cout <<"pt, eta, phi, collection "<< std::endl;          
   //std::cout <<_pT<<" "<<_eta<<" "<<_phi<<" "<<_collection<< std::endl;
   
   //
   _filterIds_n    =  ntP->triggerobject_filterIds_n->at(_ID);
   _filterLabels_n =  ntP->triggerobject_filterLabels_n->at(_ID);
   _pathNamesAll_n =  ntP->triggerobject_pathNamesAll_n->at(_ID);
   
   //std::cout <<"Number of filters, labels, paths per trigger object "<< std::endl;
   //std::cout <<_filterIds_n  <<" "<< _filterLabels_n <<" "<< _pathNamesAll_n << std::endl;
   
   _pathNamesAll_offset = 0;
   _filterLabels_offset = 0;
   _filterIds_offset    = 0;
  
   for (int i=0; i<_ID; i++)
   {
     _pathNamesAll_offset += ntP->triggerobject_pathNamesAll_n->at(i);
     _filterIds_offset    += ntP->triggerobject_filterIds_n->at(i);
     _filterLabels_offset += ntP->triggerobject_filterLabels_n->at(i);
     }   
     
   //std::cout <<"offsets "<< std::endl;
   //std::cout <<_filterIds_offset  <<" "<< _filterLabels_offset <<" "<< _pathNamesAll_offset << std::endl;
     
   //    
   std::vector<string> pathNamesAll      = *ntP->triggerobject_pathNamesAll;
   std::vector<bool> pathNamesAll_isL3   = *ntP->triggerobject_pathNamesAll_isL3;
   std::vector<bool> pathNamesAll_isLF   = *ntP->triggerobject_pathNamesAll_isLF;
   std::vector<bool> pathNamesAll_isBoth = *ntP->triggerobject_pathNamesAll_isBoth;
   std::vector<bool> pathNamesAll_isNone = *ntP->triggerobject_pathNamesAll_isNone;
   std::vector<int> filterIds            = *ntP->triggerobject_filterIds;
   std::vector<std::string> filterLabels = *ntP->triggerobject_filterLabels;

   
   _pathNamesAll = std::vector<std::string>(pathNamesAll.begin() + _pathNamesAll_offset, pathNamesAll.begin() + _pathNamesAll_offset + _pathNamesAll_n);
   
   _pathNamesAll_isL3   = std::vector<bool>(pathNamesAll_isL3.begin() + _pathNamesAll_offset, pathNamesAll_isL3.begin() + _pathNamesAll_offset + _pathNamesAll_n);
   _pathNamesAll_isLF   = std::vector<bool>(pathNamesAll_isLF.begin() + _pathNamesAll_offset, pathNamesAll_isLF.begin() + _pathNamesAll_offset + _pathNamesAll_n);
   _pathNamesAll_isBoth = std::vector<bool>(pathNamesAll_isBoth.begin() + _pathNamesAll_offset, pathNamesAll_isBoth.begin() + _pathNamesAll_offset + _pathNamesAll_n);
   _pathNamesAll_isNone = std::vector<bool>(pathNamesAll_isNone.begin() + _pathNamesAll_offset, pathNamesAll_isNone.begin() + _pathNamesAll_offset + _pathNamesAll_n);
    
   _filterIds    = std::vector<int>(filterIds.begin()+_filterIds_offset, filterIds.begin()+_filterIds_offset+_filterIds_n);
   _filterLabels = std::vector<std::string>(filterLabels.begin()+_filterLabels_offset, filterLabels.begin()+_filterLabels_offset+_filterLabels_n);

    
   /*
   
   std::cout<<"========= paths size " <<  _pathNamesAll.size() << std::endl;  
     
   for(int j=0;j<_pathNamesAll.size();j++)
   { 
    std::cout<<_pathNamesAll.at(j)<<" "<< _pathNamesAll_isL3.at(j)<<" "<< 
   	       _pathNamesAll_isLF.at(j)<<" "<<_pathNamesAll_isBoth.at(j)<<" "<< _pathNamesAll_isNone.at(j)  <<std::endl;}
  
   std::cout<<"========= filters size"<< std::endl;
   for(int j=0;j<_filterIds.size();j++)
   { 
    std::cout<< _filterIds.at(j)  <<std::endl;}
   
   std::cout<<"========= labels size "<< std::endl;
   for(int j=0;j<_filterLabels.size();j++)
   { 
    std::cout<< _filterLabels.at(j)  <<std::endl;}
   
   */
     
}


void TriggerObj::init()
{  
   
   _pT  = -999.; 
   _eta = -999.; 
   _phi = -999.; 
   _collection = "empty";
   
   _pathNamesAll_n = 0; 
   _filterLabels_n = 0;
   _filterIds_n = 0; 
   
   _pathNamesAll_offset = 0; //to decode flattrees
   _filterLabels_offset = 0; 
   _filterIds_offset = 0; 
   
   _filterIds.clear();
   _filterLabels.clear();
    
   _pathNamesAll.clear();    
   _pathNamesAll_isL3.clear();
   _pathNamesAll_isLF.clear();
   _pathNamesAll_isBoth.clear();
   _pathNamesAll_isNone.clear();
   	
}



bool TriggerObj::sel()
{  
  bool isTrigger = false;
  
 
  //if (std::find(_filterIds.begin(),_filterIds.end(), 82)!= _filterIds.end() || std::find(_filterIds.begin(),_filterIds.end(), 83) != _filterIds.end()) isTrigger = true;
  
  
  for(int j=0;j<_pathNamesAll_n;j++)
  { 
    //SL
    std::size_t ok1 = _pathNamesAll.at(j).find("HLT_IsoMu20");
    std::size_t ok2 = _pathNamesAll.at(j).find("HLT_IsoTkMu20");
    std::size_t ok3 = _pathNamesAll.at(j).find("HLT_Ele23_CaloIdL_TrackIdL_IsoVL");
    std::size_t ok4 = _pathNamesAll.at(j).find("HLT_Ele23_WPLoose_Gsf");
    
    //Tri-lep
    std::size_t ok5 = _pathNamesAll.at(j).find("HLT_DiMu9_Ele9_CaloIdL_TrackIdL");
    std::size_t ok6 = _pathNamesAll.at(j).find("HLT_Mu8_DiEle12_CaloIdL_TrackIdL");
    std::size_t ok7 = _pathNamesAll.at(j).find("HLT_TripleMu_12_10_5");
    std::size_t ok8 = _pathNamesAll.at(j).find("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL");
    
    //Di-lep
    std::size_t ok9  = _pathNamesAll.at(j).find("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL");
    std::size_t ok10 = _pathNamesAll.at(j).find("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL");
    std::size_t ok11 = _pathNamesAll.at(j).find("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL");
    std::size_t ok12 = _pathNamesAll.at(j).find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL");
    std::size_t ok13 = _pathNamesAll.at(j).find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL");
    
    std::size_t ok14 = _pathNamesAll.at(j).find("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter");
    std::size_t ok15 = _pathNamesAll.at(j).find("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter");
    
    //
    if ( ok1!=std::string::npos  || ok2!=std::string::npos  || ok3!=std::string::npos  || ok4!=std::string::npos  ||
         ok5!=std::string::npos  || ok6!=std::string::npos  || ok7!=std::string::npos  || ok8!=std::string::npos  ||
         ok9!=std::string::npos  || ok10!=std::string::npos || ok11!=std::string::npos || ok12!=std::string::npos ||
	 ok13!=std::string::npos || ok14!=std::string::npos || ok15!=std::string::npos) 
    {
      isTrigger = true;
      break;
     }
    }
          
 
  //std::cout <<"isTrigger " << isTrigger << std::endl;
  
  return isTrigger;
}
