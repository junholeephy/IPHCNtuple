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

#include "../include/SignalExtractionMVA.h"

//https://root.cern.ch/doc/v606/classTMVA_1_1Reader.html

TMVA::Reader* Book_2LSS_TT_MVAReader(std::string basePath, std::string weightFileName, std::string type)
{

    // weight/2lss_ttbar_BDTG.weights.xml
    TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

    reader->AddVariable("max(abs(LepGood_eta[iF_Recl[0]]),abs(LepGood_eta[iF_Recl[1]]))",  &max_Lep_eta      );
    reader->AddVariable("nJet25_Recl",                                                     &nJet25_Recl      );
    reader->AddVariable("mindr_lep1_jet",                                                  &mindr_lep1_jet   );
    reader->AddVariable("mindr_lep2_jet",                                                  &mindr_lep2_jet   );
    reader->AddVariable("min(met_pt,400)",                                                 &met              );
    reader->AddVariable("avg_dr_jet",                                                      &avg_dr_jet       );
    reader->AddVariable("MT_met_lep1",                                                     &MT_met_lep1      );

    reader->AddSpectator("iF_Recl[0]",                                                     &iF_Recl0         );
    reader->AddSpectator("iF_Recl[1]",                                                     &iF_Recl1         );
    reader->AddSpectator("iF_Recl[2]",                                                     &iF_Recl2         );

    reader->BookMVA("BDTG method", basePath+"/"+weightFileName);

    return reader;
}

TMVA::Reader* Book_2LSS_TTV_MVAReader(std::string basePath, std::string weightFileName, std::string type)
{

    // weight/2lss_ttV_BDTG.weights.xml
    TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

    reader->AddVariable("max(abs(LepGood_eta[iF_Recl[0]]),abs(LepGood_eta[iF_Recl[1]]))", &max_Lep_eta      );
    reader->AddVariable("MT_met_lep1",                                                    &MT_met_lep1      );
    reader->AddVariable("nJet25_Recl",                                                    &nJet25_Recl      );
    reader->AddVariable("mindr_lep1_jet",                                                 &mindr_lep1_jet   );
    reader->AddVariable("mindr_lep2_jet",                                                 &mindr_lep2_jet   );
    reader->AddVariable("LepGood_conePt[iF_Recl[0]]",                                     &LepGood_conePt0  );
    reader->AddVariable("LepGood_conePt[iF_Recl[1]]",                                     &LepGood_conePt1  );

    reader->AddSpectator("iF_Recl[0]",                                                    &iF_Recl0         );
    reader->AddSpectator("iF_Recl[1]",                                                    &iF_Recl1         );
    reader->AddSpectator("iF_Recl[2]",                                                    &iF_Recl2         );

    reader->BookMVA("BDTG method", basePath+"/"+weightFileName);

    return reader;
}

TMVA::Reader* Book_3L_TT_MVAReader(std::string basePath, std::string weightFileName, std::string type)
{

    // weight/3l_ttbar_BDTG.weights.xml
    TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

    reader->AddVariable("max(abs(LepGood_eta[iF_Recl[0]]),abs(LepGood_eta[iF_Recl[1]]))", &max_Lep_eta      );
    reader->AddVariable("MT_met_lep1",                                                    &MT_met_lep1      );
    reader->AddVariable("nJet25_Recl",                                                    &nJet25_Recl      );
    reader->AddVariable("mhtJet25_Recl",                                                  &mhtJet25_Recl    );
    reader->AddVariable("avg_dr_jet",                                                     &avg_dr_jet       );
    reader->AddVariable("mindr_lep1_jet",                                                 &mindr_lep1_jet   );
    reader->AddVariable("mindr_lep2_jet",                                                 &mindr_lep2_jet   );

    reader->AddSpectator("iF_Recl[0]",                                                    &iF_Recl0         );
    reader->AddSpectator("iF_Recl[1]",                                                    &iF_Recl1         );
    reader->AddSpectator("iF_Recl[2]",                                                    &iF_Recl2         );

    reader->BookMVA("BDTG method", basePath+"/"+weightFileName);

    return reader;
}

TMVA::Reader* Book_3L_TTV_MVAReader(std::string basePath, std::string weightFileName, std::string type)
{

    // weight/3l_ttV_BDTG.weights.xml
    TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

    reader->AddVariable("max(abs(LepGood_eta[iF_Recl[0]]),abs(LepGood_eta[iF_Recl[1]]))", &max_Lep_eta      );
    reader->AddVariable("MT_met_lep1",                                                    &MT_met_lep1      );
    reader->AddVariable("nJet25_Recl",                                                    &nJet25_Recl      );
    reader->AddVariable("mindr_lep1_jet",                                                 &mindr_lep1_jet   );
    reader->AddVariable("mindr_lep2_jet",                                                 &mindr_lep2_jet   );
    reader->AddVariable("LepGood_conePt[iF_Recl[0]]",                                     &LepGood_conePt0  );
    reader->AddVariable("LepGood_conePt[iF_Recl[2]]",                                     &LepGood_conePt0  );

    reader->AddSpectator("iF_Recl[0]",                                                    &iF_Recl0         );
    reader->AddSpectator("iF_Recl[1]",                                                    &iF_Recl1         );
    reader->AddSpectator("iF_Recl[2]",                                                    &iF_Recl2         );

    reader->BookMVA("BDTG method", basePath+"/"+weightFileName);

    return reader;
}

void Load_MVA()
{

    //std::cout << "Temporarily redirecting stdout to avoid huge TMVA dump when loading MVA readers..." << std::endl;
    std::stringstream tmpBuffer;
    std::streambuf* oldStdout = std::cout.rdbuf(tmpBuffer.rdbuf());

    std::string NtupleAnalyzerMVAPath = std::string("../src/");
    mva_2lss_tt  = Book_2LSS_TT_MVAReader(  NtupleAnalyzerMVAPath, "weight/2lss_ttbar_BDTG.weights.xml", "test");
    mva_2lss_ttV = Book_2LSS_TTV_MVAReader( NtupleAnalyzerMVAPath, "weight/2lss_ttV_BDTG.weights.xml",   "test");
    mva_3l_tt    = Book_3L_TT_MVAReader(    NtupleAnalyzerMVAPath, "weight/3l_ttbar_BDTG.weights.xml",   "test");
    mva_3l_ttV   = Book_3L_TTV_MVAReader(   NtupleAnalyzerMVAPath, "weight/3l_ttV_BDTG.weights.xml",     "test");
    //mva_2lss_tt  = Book_2LSS_TT_MVAReader(NtupleAnalyzerMVAPath,  "weight/2lss_ttbar_BDTG.weights.xml", "test");
    //mva_2lss_ttV = Book_2LSS_TTV_MVAReader(NtupleAnalyzerMVAPath, "weight/2lss_ttV_BDTG.weights.xml",   "test");
    //mva_3l_tt    = Book_3L_TT_MVAReader(NtupleAnalyzerMVAPath,    "weight/3l_ttbar_BDTG.weights.xml",   "test");
    //mva_3l_ttV   = Book_3L_TTV_MVAReader(NtupleAnalyzerMVAPath,   "weight/3l_ttV_BDTG.weights.xml",     "test");

    std::cout.rdbuf(oldStdout);
    //std::cout << "Stdout now restored." << std::endl;

    return;
}

