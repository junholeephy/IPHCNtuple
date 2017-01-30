#ifndef MVA
#define MVA

#include <memory>
#include <iostream>
#include <sstream>

#include "TRegexp.h"
#include "TString.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TFile.h"
#include "TThreadSlots.h"
#include "TROOT.h"
#include "Compression.h"

// Input variables
float max_Lep_eta;
float mindr_lep1_jet;
float mindr_lep2_jet;
float met;
float avg_dr_jet;
float MT_met_lep1;
float LepGood_conePt0;
float LepGood_conePt1;
float nJet25_Recl;
float mhtJet25_Recl;

// Input variables (spectators)
float iF_Recl0;
float iF_Recl1;
float iF_Recl2;

// Output variables
float signal_2lss_TT_MVA;
float signal_2lss_TTV_MVA;
float signal_3l_TT_MVA;
float signal_3l_TTV_MVA;

// TMVA readers
TMVA::Reader* mva_2lss_tt;
TMVA::Reader* mva_2lss_ttV;
TMVA::Reader* mva_3l_tt;
TMVA::Reader* mva_3l_ttV;

// Function definitions
TMVA::Reader* Book_2LSS_TT_MVAReader(  std::string basePath, std::string weightFileName, std::string type);
TMVA::Reader* Book_2LSS_TTV_MVAReader( std::string basePath, std::string weightFileName, std::string type);
TMVA::Reader* Book_3L_TT_MVAReader(    std::string basePath, std::string weightFileName, std::string type);
TMVA::Reader* Book_3L_TTV_MVAReader(   std::string basePath, std::string weightFileName, std::string type);
void Load_MVA();

#endif
