#include "../include/TTbarHiggsMultileptonAnalysis.h"
#include "TSystem.h"
#include "SignalExtractionMVA.cxx"
#include "Helper.cxx"
#include "BTagging.cxx"

#define kCat_3l_2b_2j 0
#define kCat_3l_1b_2j 1
#define kCat_3l_2b_1j 2
#define kCat_3l_1b_1j 3
#define kCat_3l_2b_0j 4
#define kCat_4l_2b 5
#define kCat_4l_1b 6
#define kCat_2lss_2b_4j 7
#define kCat_2lss_1b_4j 8
#define kCat_2lss_2b_3j 9
#define kCat_2lss_1b_3j 10
#define kCat_2lss_2b_2j 11

using namespace std;

TTbarHiggsMultileptonAnalysis::TTbarHiggsMultileptonAnalysis() 
{
    _printLHCO_MC = false;
    _processLHCO_MC = -1;

    _printLHCO_RECO = false;
    _processLHCO_RECO = -1;

}

TTbarHiggsMultileptonAnalysis::TTbarHiggsMultileptonAnalysis(TString inputFileName, TChain *tree, TString sampleName, TString treeName, TString outputFileName, bool isdata, float xsec, float lumi, int nowe, int nmax)
{    

    //
    _isdata = isdata;
    _xsec = xsec;
    _lumi = lumi;
    _nowe = nowe;
    _nmax = nmax;
    _outputFileName = outputFileName;
    _sampleName = sampleName;

    //
    _printLHCO_MC = false;
    _processLHCO_MC = -1;

    _printLHCO_RECO = false;
    _processLHCO_RECO = -1;

    //
    //_file_PVreweighting = TFile::Open("/home-pbs/lebihan/someone/ttH_070116/ttH/NtupleAnalyzer/test/PUweight.root");
    //_h_PV = (TH1F*)_file_PVreweighting->Get("PU_reweighting");

    //
    tree = new TChain(treeName.Data());

    std::ifstream infile;
    infile.open(inputFileName.Data());
    std::string ifile = "";
    while( getline(infile, ifile) )
    {
        std::string fnameStr = std::string(ifile);
        tree->Add(fnameStr.c_str());
    }   
    infile.close();
    Init(tree);

    theHistoManager = new HistoManager();

    TString outputfileNameRoot = _outputFileName+".root";
    outputfile = new TFile(outputfileNameRoot.Data(), "recreate");  

}

void TTbarHiggsMultileptonAnalysis::InitLHCO(int process_MC, int process_RECO) 
{  
    if (!_isdata)
    {    
        _printLHCO_MC = true;
        _processLHCO_MC = process_MC;     
        TString fout_MC_path(_outputFileName+"_LHCO_MC.txt"); 
        fout_MC.open(fout_MC_path.Data());}


        _printLHCO_RECO = true;
        _processLHCO_RECO = process_RECO;
        TString fout_RECO_path(_outputFileName+"_LHCO_RECO.txt");
        fout_RECO.open(fout_RECO_path.Data());


        fline00 = "#   typ	  eta	 phi	   pt  jmass  ntrk  btag   had/em  dummy dummy";
        del = "    ";
        trig = "8";   
}

void TTbarHiggsMultileptonAnalysis::createHistograms()
{    
    outputfile->cd();
    initializeOutputTree();

    // General
    theHistoManager->addHisto("CutFlow",                                     "noSel",        "",   "",  10,   0,     10);

    // Preselection variables

    theHistoManager->addHisto("MuonPt",                                      "noSel",        "",   "",   25,   0,   200);
    theHistoManager->addHisto("MuonEta",                                     "noSel",        "",   "",   25,  -3,     3);
    theHistoManager->addHisto("MuonMVA",                                     "noSel",        "",   "",   20,   0,     1);
    theHistoManager->addHisto("ElectronPt",                                  "noSel",        "",   "",   25,   0,   200);
    theHistoManager->addHisto("ElectronEta",                                 "noSel",        "",   "",   25,  -3,     3);
    theHistoManager->addHisto("ElectronMVA",                                 "noSel",        "",   "",   20,   0,     1);
    theHistoManager->addHisto("TauPt",                                       "noSel",        "",   "",   25,   0,   200);
    theHistoManager->addHisto("TauEta",                                      "noSel",        "",   "",   25,  -3,     3);
    theHistoManager->addHisto("JetPt",                                       "noSel",        "",   "",   25,   0,   200);
    theHistoManager->addHisto("JetEta",                                      "noSel",        "",   "",   25,  -3,     3);
    theHistoManager->addHisto("JetCSVv2",                                    "noSel",        "",   "",   20,   0,     1);
    theHistoManager->addHisto("MET",                                         "noSel",        "",   "",   50,   0,   200);

    // Signal Region Two Leptons

    theHistoManager->addHisto2D("SelectedLeptonsVsJets",           "noSel",   "ttH2lss",   "",    8,    0,    7,    8,    0,    7);
    theHistoManager->addHisto2D("SelectedLeptonsVsBJets",          "noSel",   "ttH2lss",   "",    8,    0,    7,    8,    0,    7);

    theHistoManager->addHisto("WeightCSV_min",                                    "",     "ttH2l",   "",    30,    0,    3);
    theHistoManager->addHisto("WeightCSV_max",                                    "",     "ttH2l",   "",    30,    0,    3);

    theHistoManager->addHisto("CutFlow",                          "FullThreeLeptons",   "ttH2lss",   "",    10,   -1,    9);

    theHistoManager->addHisto("CutFlow",                          "FullThreeLeptons",   "ttH2lee",   "",    10,   -1,    9);
    theHistoManager->addHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "ttH2lee",   "",    20,    0,  200);
    theHistoManager->addHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "ttH2lee",   "",    20,    0,  200);
    theHistoManager->addHisto("MET",                              "FullThreeLeptons",   "ttH2lee",   "",    20,    0,  200);
    theHistoManager->addHisto("MHT",                              "FullThreeLeptons",   "ttH2lee",   "",    20,    0,  200);
    theHistoManager->addHisto("MetLD",                            "FullThreeLeptons",   "ttH2lee",   "",    15, -0.2,  1.4);
    theHistoManager->addHisto("TauMultiplicity",                  "FullThreeLeptons",   "ttH2lee",   "",    10,   -1,    9);
    theHistoManager->addHisto("JetMultiplicity",                  "FullThreeLeptons",   "ttH2lee",   "",    10,   -1,    9);
    theHistoManager->addHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "ttH2lee",   "",    10,   -1,    9);
    theHistoManager->addHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "ttH2lee",   "",    10,   -1,    9);
    theHistoManager->addHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "ttH2lee",   "",    10,   -1,    9);

    theHistoManager->addHisto("CutFlow",                          "FullThreeLeptons",   "ttH2lem",   "",    10,   -1,    9);
    theHistoManager->addHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "ttH2lem",   "",    20,    0,  200);
    theHistoManager->addHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "ttH2lem",   "",    20,    0,  200);
    theHistoManager->addHisto("MET",                              "FullThreeLeptons",   "ttH2lem",   "",    20,    0,  200);
    theHistoManager->addHisto("MHT",                              "FullThreeLeptons",   "ttH2lem",   "",    20,    0,  200);
    theHistoManager->addHisto("MetLD",                            "FullThreeLeptons",   "ttH2lem",   "",    15, -0.2,  1.4);
    theHistoManager->addHisto("TauMultiplicity",                  "FullThreeLeptons",   "ttH2lem",   "",    10,   -1,    9);
    theHistoManager->addHisto("JetMultiplicity",                  "FullThreeLeptons",   "ttH2lem",   "",    10,   -1,    9);
    theHistoManager->addHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "ttH2lem",   "",    10,   -1,    9);
    theHistoManager->addHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "ttH2lem",   "",    10,   -1,    9);
    theHistoManager->addHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "ttH2lem",   "",    10,   -1,    9);

    theHistoManager->addHisto("CutFlow",                          "FullThreeLeptons",   "ttH2lmm",   "",    10,   -1,    9);
    theHistoManager->addHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "ttH2lmm",   "",    20,    0,  200);
    theHistoManager->addHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "ttH2lmm",   "",    20,    0,  200);
    theHistoManager->addHisto("MET",                              "FullThreeLeptons",   "ttH2lmm",   "",    20,    0,  200);
    theHistoManager->addHisto("MHT",                              "FullThreeLeptons",   "ttH2lmm",   "",    20,    0,  200);
    theHistoManager->addHisto("MetLD",                            "FullThreeLeptons",   "ttH2lmm",   "",    15, -0.2,  1.4);
    theHistoManager->addHisto("TauMultiplicity",                  "FullThreeLeptons",   "ttH2lmm",   "",    10,   -1,    9);
    theHistoManager->addHisto("JetMultiplicity",                  "FullThreeLeptons",   "ttH2lmm",   "",    10,   -1,    9);
    theHistoManager->addHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "ttH2lmm",   "",    10,   -1,    9);
    theHistoManager->addHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "ttH2lmm",   "",    10,   -1,    9);
    theHistoManager->addHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "ttH2lmm",   "",    10,   -1,    9);

    // TT variables

    theHistoManager->addHisto2D("SelectedLeptonsVsJets",           "noSel",   "TT_2l_CR",   "",    8,    0,    7,    8,    0,    7);
    theHistoManager->addHisto2D("SelectedLeptonsVsBJets",          "noSel",   "TT_2l_CR",   "",    8,    0,    7,    8,    0,    7);

    theHistoManager->addHisto("CutFlow",                          "FullThreeLeptons",   "TT_2l_CR",   "",    10,   -1,    9);
    theHistoManager->addHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "TT_2l_CR",   "",    20,    0,  200);
    theHistoManager->addHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "TT_2l_CR",   "",    20,    0,  200);
    theHistoManager->addHisto("MET",                              "FullThreeLeptons",   "TT_2l_CR",   "",    20,    0,  200);
    theHistoManager->addHisto("MetLD",                            "FullThreeLeptons",   "TT_2l_SB",   "",    15, -0.2,  1.4);
    theHistoManager->addHisto("TauMultiplicity",                  "FullThreeLeptons",   "TT_2l_CR",   "",    10,   -1,    9);
    theHistoManager->addHisto("JetMultiplicity",                  "FullThreeLeptons",   "TT_2l_CR",   "",    10,   -1,    9);
    theHistoManager->addHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "TT_2l_CR",   "",    10,   -1,    9);
    theHistoManager->addHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "TT_2l_CR",   "",    10,   -1,    9);
    theHistoManager->addHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "TT_2l_CR",   "",    10,   -1,    9);

    // Signal Region Three Leptons

    theHistoManager->addHisto("CutFlow",                          "FullThreeLeptons",   "ttH3l",   "",   10,   -1,    9);

    theHistoManager->addHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "ttH3l",   "",   20,    0,  200);
    theHistoManager->addHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "ttH3l",   "",   20,    0,  200);
    theHistoManager->addHisto("ThirdLeptonPt",                    "FullThreeLeptons",   "ttH3l",   "",   20,    0,  200);
    theHistoManager->addHisto("MET",                              "FullThreeLeptons",   "ttH3l",   "",   20,    0,  200);
    theHistoManager->addHisto("MHT",                              "FullThreeLeptons",   "ttH3l",   "",   20,    0,  200);
    theHistoManager->addHisto("MetLD",                            "FullThreeLeptons",   "ttH3l",   "",   15, -0.2,  1.4);
    theHistoManager->addHisto("TauMultiplicity",                  "FullThreeLeptons",   "ttH3l",   "",   10,   -1,    9);
    theHistoManager->addHisto("JetMultiplicity",                  "FullThreeLeptons",   "ttH3l",   "",   10,   -1,    9);
    theHistoManager->addHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "ttH3l",   "",   10,   -1,    9);
    theHistoManager->addHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "ttH3l",   "",   10,   -1,    9);
    theHistoManager->addHisto("SumOfLeptonsCharges",              "FullThreeLeptons",   "ttH3l",   "",   10,   -5,    5);
    theHistoManager->addHisto("SumOfThreeLeptonsCharges",         "FullThreeLeptons",   "ttH3l",   "",   10,   -5,    5);
    theHistoManager->addHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "ttH3l",   "",   10,   -1,    9);

    // cut 1: / cut 2: / cut 3: ...

    // WZ variables

    theHistoManager->addHisto("CutFlow",                                     "noSel", "WZZZ_CR",   "",   10,   -1,    9);
    theHistoManager->addHisto("CutFlow",                          "ThreePreselected", "WZZZ_CR",   "",    3,   -1,    2);
    theHistoManager->addHisto("CutFlow",                           "PassingTightMVA", "WZZZ_CR",   "",    3,   -1,    2);
    theHistoManager->addHisto("CutFlow",                                  "PassingZ", "WZZZ_CR",   "",    3,   -1,    2);
    theHistoManager->addHisto("CutFlow",                              "PassingMETLD", "WZZZ_CR",   "",    3,   -1,    2);
    theHistoManager->addHisto("CutFlow",                               "PassingJets", "WZZZ_CR",   "",    3,   -1,    2);
    theHistoManager->addHisto("CutFlow",                        "PassingMETLDorJets", "WZZZ_CR",   "",    3,   -1,    2);
    theHistoManager->addHisto("CutFlow",                        "PassingMETLDorSFOS", "WZZZ_CR",   "",    3,   -1,    2);
    theHistoManager->addHisto("CutFlow",                          "PassingbJetsVeto", "WZZZ_CR",   "",    3,   -1,    2);

    theHistoManager->addHisto("ZCandidateInvariantMass",                  "PassingZ", "WZZZ_CR",   "",  15,   60,   120);
    theHistoManager->addHisto("ZCandidateInvariantMass",        "PassingMETLDorJets", "WZZZ_CR",   "",  15,   60,   120);
    theHistoManager->addHisto("ZCandidateInvariantMass",        "PassingMETLDorSFOS", "WZZZ_CR",   "",  15,   60,   120);
    theHistoManager->addHisto("ZCandidateInvariantMass",          "PassingbJetsVeto", "WZZZ_CR",   "",  15,   60,   120);

    theHistoManager->addHisto("ZCandidateInvariantMass",                  "FinalCut", "WZZZ_CR",   "",  15,   60,   120);
    theHistoManager->addHisto("ZCandidateTransverseMomentum",             "FinalCut", "WZZZ_CR",   "",  12,    0,   500);
    theHistoManager->addHisto("MET",                                      "FinalCut", "WZZZ_CR",   "",  20,    0,   200);
    theHistoManager->addHisto("MHT",                                      "FinalCut", "WZZZ_CR",   "",  20,    0,   200);
    theHistoManager->addHisto("MetLD",                                    "FinalCut", "WZZZ_CR",   "",  15, -0.2,   1.4);
    theHistoManager->addHisto("TauMultiplicity",                          "FinalCut", "WZZZ_CR",   "",   8,    0,     8);
    theHistoManager->addHisto("JetMultiplicity",                          "FinalCut", "WZZZ_CR",   "",   8,    0,     8);

    theHistoManager->addHisto("MTW",                                      "FinalCut", "WZZZ_CR",   "",  15,    0,   300);
    theHistoManager->addHisto("InvariantMassOfSelectedLeptons",           "FinalCut", "WZZZ_CR",   "",  15,    0,   600);
    theHistoManager->addHisto("SumOfSelectedLeptonsCharges",              "FinalCut", "WZZZ_CR",   "",  10,   -5,     5);
    theHistoManager->addHisto("SumVecPtSelectedLeptons",                  "FinalCut", "WZZZ_CR",   "",  12,    0,   500);

    theHistoManager->addHisto("InvMassRemainingLepton",                   "FinalCut", "WZZZ_CR",   "",  15,    0,   300);
    theHistoManager->addHisto("SumOfRemainingLeptonsCharges",             "FinalCut", "WZZZ_CR",   "",  10,   -5,     5);
    theHistoManager->addHisto("SumVecPtRemainingLeptons",                 "FinalCut", "WZZZ_CR",   "",  25,    0,   500);

    // WZ variables - relaxed CR

    theHistoManager->addHisto("CutFlow",                                     "noSel", "WZrel_CR",   "",   10,   -1,    9);
    theHistoManager->addHisto("CutFlow",                          "ThreePreselected", "WZrel_CR",   "",    3,   -1,    2);
    theHistoManager->addHisto("CutFlow",                           "PassingTightMVA", "WZrel_CR",   "",    3,   -1,    2);
    theHistoManager->addHisto("CutFlow",                                  "PassingZ", "WZrel_CR",   "",    3,   -1,    2);
    theHistoManager->addHisto("CutFlow",                              "PassingMETLD", "WZrel_CR",   "",    3,   -1,    2);
    theHistoManager->addHisto("CutFlow",                               "PassingJets", "WZrel_CR",   "",    3,   -1,    2);
    theHistoManager->addHisto("CutFlow",                        "PassingMETLDorJets", "WZrel_CR",   "",    3,   -1,    2);
    theHistoManager->addHisto("CutFlow",                        "PassingMETLDorSFOS", "WZrel_CR",   "",    3,   -1,    2);
    theHistoManager->addHisto("CutFlow",                          "PassingbJetsVeto", "WZrel_CR",   "",    3,   -1,    2);

    theHistoManager->addHisto("ZCandidateInvariantMass",                  "PassingZ", "WZrel_CR",   "",  15,   60,   120);
    theHistoManager->addHisto("ZCandidateInvariantMass",        "PassingMETLDorJets", "WZrel_CR",   "",  15,   60,   120);
    theHistoManager->addHisto("ZCandidateInvariantMass",        "PassingMETLDorSFOS", "WZrel_CR",   "",  15,   60,   120);
    theHistoManager->addHisto("ZCandidateInvariantMass",          "PassingbJetsVeto", "WZrel_CR",   "",  15,   60,   120);

    theHistoManager->addHisto("ZCandidateInvariantMass",                  "FinalCut", "WZrel_CR",   "",  15,   60,   120);
    theHistoManager->addHisto("ZCandidateTransverseMomentum",             "FinalCut", "WZrel_CR",   "",  12,    0,   500);
    theHistoManager->addHisto("MET",                                      "FinalCut", "WZrel_CR",   "",  20,    0,   200);
    theHistoManager->addHisto("MHT",                                      "FinalCut", "WZrel_CR",   "",  20,    0,   200);
    theHistoManager->addHisto("MetLD",                                    "FinalCut", "WZrel_CR",   "",  15, -0.2,   1.4);
    theHistoManager->addHisto("TauMultiplicity",                          "FinalCut", "WZrel_CR",   "",   8,    0,     8);
    theHistoManager->addHisto("JetMultiplicity",                          "FinalCut", "WZrel_CR",   "",   8,    0,     8);

    theHistoManager->addHisto("InvariantMassOfSelectedLeptons",           "FinalCut", "WZrel_CR",   "",  30,    0,   600);
    theHistoManager->addHisto("SumOfSelectedLeptonsCharges",              "FinalCut", "WZrel_CR",   "",  10,   -5,     5);
    theHistoManager->addHisto("SumVecPtSelectedLeptons",                  "FinalCut", "WZrel_CR",   "",  25,    0,   500);

    theHistoManager->addHisto("InvMassRemainingLepton",                   "FinalCut", "WZrel_CR",   "",  15,    0,   300);
    theHistoManager->addHisto("SumOfRemainingLeptonsCharges",             "FinalCut", "WZrel_CR",   "",  10,   -5,     5);
    theHistoManager->addHisto("SumVecPtRemainingLeptons",                 "FinalCut", "WZrel_CR",   "",  25,    0,   500);

    // TTZ variables - relaxed CR

    theHistoManager->addHisto("CutFlow",                                     "noSel",   "TTZ_CR",   "",   10,   -1,    9);
    theHistoManager->addHisto("CutFlow",                          "ThreePreselected",   "TTZ_CR",   "",    3,   -1,    2);
    theHistoManager->addHisto("CutFlow",                           "PassingTightMVA",   "TTZ_CR",   "",    3,   -1,    2);
    theHistoManager->addHisto("CutFlow",                                  "PassingZ",   "TTZ_CR",   "",    3,   -1,    2);
    theHistoManager->addHisto("CutFlow",                              "PassingMETLD",   "TTZ_CR",   "",    3,   -1,    2);
    theHistoManager->addHisto("CutFlow",                               "PassingJets",   "TTZ_CR",   "",    3,   -1,    2);
    theHistoManager->addHisto("CutFlow",                        "PassingMETLDorJets",   "TTZ_CR",   "",    3,   -1,    2);
    theHistoManager->addHisto("CutFlow",                        "PassingMETLDorSFOS",   "TTZ_CR",   "",    3,   -1,    2);
    theHistoManager->addHisto("CutFlow",                          "PassingbJetsVeto",   "TTZ_CR",   "",    3,   -1,    2);

    theHistoManager->addHisto("ZCandidateInvariantMass",                  "PassingZ",   "TTZ_CR",   "",  15,   60,   120);
    theHistoManager->addHisto("ZCandidateInvariantMass",        "PassingMETLDorJets",   "TTZ_CR",   "",  15,   60,   120);
    theHistoManager->addHisto("ZCandidateInvariantMass",        "PassingMETLDorSFOS",   "TTZ_CR",   "",  15,   60,   120);
    theHistoManager->addHisto("ZCandidateInvariantMass",          "PassingbJetsVeto",   "TTZ_CR",   "",  15,   60,   120);

    theHistoManager->addHisto("LeadingLeptonPt",                          "FinalCut",   "TTZ_CR",   "",  20,    0,   200);
    theHistoManager->addHisto("SubleadingLeptonPt",                       "FinalCut",   "TTZ_CR",   "",  20,    0,   200);

    theHistoManager->addHisto("ZCandidateInvariantMass",                  "FinalCut",   "TTZ_CR",   "",  15,   60,   120);
    theHistoManager->addHisto("ZCandidateTransverseMomentum",             "FinalCut",   "TTZ_CR",   "",  12,    0,   500);
    theHistoManager->addHisto("MET",                                      "FinalCut",   "TTZ_CR",   "",  20,    0,   200);
    theHistoManager->addHisto("MHT",                                      "FinalCut",   "TTZ_CR",   "",  20,    0,   200);
    theHistoManager->addHisto("MetLD",                                    "FinalCut",   "TTZ_CR",   "",  15, -0.2,   1.4);
    theHistoManager->addHisto("TauMultiplicity",                          "FinalCut",   "TTZ_CR",   "",   8,    0,     8);
    theHistoManager->addHisto("JetMultiplicity",                          "FinalCut",   "TTZ_CR",   "",   8,    0,     8);

    theHistoManager->addHisto("InvariantMassOfSelectedLeptons",           "FinalCut",   "TTZ_CR",   "",  30,    0,   600);
    theHistoManager->addHisto("SumOfSelectedLeptonsCharges",              "FinalCut",   "TTZ_CR",   "",  10,   -5,     5);
    theHistoManager->addHisto("SumVecPtSelectedLeptons",                  "FinalCut",   "TTZ_CR",   "",  25,    0,   500);

    theHistoManager->addHisto("InvMassRemainingLepton",                   "FinalCut",   "TTZ_CR",   "",  15,    0,   300);
    theHistoManager->addHisto("SumOfRemainingLeptonsCharges",             "FinalCut",   "TTZ_CR",   "",  10,   -5,     5);
    theHistoManager->addHisto("SumVecPtRemainingLeptons",                 "FinalCut",   "TTZ_CR",   "",  25,    0,   500);

    theHistoManager->addHisto("LeadingLeptonPt",                          "FinalCut", "TTZ4j_CR",   "",  20,    0,   200);
    theHistoManager->addHisto("SubleadingLeptonPt",                       "FinalCut", "TTZ4j_CR",   "",  20,    0,   200);

    theHistoManager->addHisto("MET",                                      "FinalCut", "TTZ4j_CR",   "",  20,    0,   200);
    theHistoManager->addHisto("JetMultiplicity",                          "FinalCut", "TTZ4j_CR",   "",   8,    0,     8);
    theHistoManager->addHisto("ZCandidateInvariantMass",                  "FinalCut", "TTZ4j_CR",   "",  15,   60,   120);

    // PU reweighting
    theHistoManager->addHisto("NumberOfPrimaryVertex",                       "noSel",        "",   "", 100,    0,    99);

    // 2D histo

    theHistoManager->addHisto2D("SelectedLeptonsVsJets",                     "noSel",   "ttH3l",   "",    8,    0,    7,    8,    0,    7);
    theHistoManager->addHisto2D("SelectedLeptonsVsBJets",                    "noSel",   "ttH3l",   "",    8,    0,    7,    8,    0,    7);
    theHistoManager->addHisto2D("SelectedLeptonsVsJets",                     "noSel", "WZZZ_CR",   "",    8,    0,    7,    8,    0,    7);
    theHistoManager->addHisto2D("SelectedLeptonsVsBJets",                    "noSel", "WZZZ_CR",   "",    8,    0,    7,    8,    0,    7);

    theHistoManager->addHisto2D("InvMassLastLeptonVSZMass",               "FinalCut", "WZZZ_CR",   "",   15,    0,  300,   15,   60,  120);
    theHistoManager->addHisto2D("SumPtLepVSZMass",                        "FinalCut", "WZZZ_CR",   "",   25,    0,  500,   15,   60,  120);
    theHistoManager->addHisto2D("METLDVSZMass",                           "FinalCut", "WZZZ_CR",   "",   15, -0.2,  1.4,   15,   60,  120);

    // from Daniel 

    theHistoManager->addHisto("nLep",         "PreSel", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("nLep loose",   "PreSel", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("lep1Pt",       "PreSel", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep2Pt",       "PreSel", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep3Pt",       "PreSel", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1Eta",      "PreSel", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep2Eta",      "PreSel", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep3Eta",      "PreSel", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lepQ",         "PreSel", "", "", 9, -4.5, 4.5);
    theHistoManager->addHisto("Mll"     ,     "PreSel", "", "", 10, 0., 200.);
    theHistoManager->addHisto("Mll"     ,     "PreSel 3l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("Mll"     ,     "PreSel btag", "", "", 10, 0., 200.);
    theHistoManager->addHisto("Mll"     ,     "PreSel btag 3l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("nJets",        "PreSel", "", "", 8, -0.5, 7.5);
    theHistoManager->addHisto("nLooseB",      "PreSel", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("nMediumB",     "PreSel", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("JetPt",        "PreSel", "", "", 28, 20., 300.);
    theHistoManager->addHisto("JetEta",       "PreSel", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("CSVv2" ,       "PreSel", "", "", 102, -0.01, 1.01);
    theHistoManager->addHisto("METpx",        "PreSel", "", "",14,-140.,140.);
    theHistoManager->addHisto("METpy",        "PreSel", "", "",14,-140.,140.);
    theHistoManager->addHisto("MET"  ,        "PreSel", "", "",10,0.,400.);
    theHistoManager->addHisto("METphi",       "PreSel", "", "",16,-3.2,3.2);
    theHistoManager->addHisto("METsum",       "PreSel", "", "",15,0.,3000.);
    theHistoManager->addHisto("MHT",          "PreSel", "", "", 10, 0., 200.);
    theHistoManager->addHisto("MET LD",       "PreSel", "", "", 10, 0., 2.);

    // TTH3l
    theHistoManager->addHisto("nLep",         "TTH3l", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("nFake",        "TTH3l", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("lep1Pt",       "TTH3l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep2Pt",       "TTH3l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep3Pt",       "TTH3l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep4Pt",       "TTH3l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1Eta",      "TTH3l", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep2Eta",      "TTH3l", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep3Eta",      "TTH3l", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep4Eta",      "TTH3l", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lepQ",         "TTH3l", "", "", 9, -4.5, 4.5);
    theHistoManager->addHisto("nJets",        "TTH3l", "", "", 8, -0.5, 7.5);
    theHistoManager->addHisto("nLooseB",      "TTH3l", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("nMediumB",     "TTH3l", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("JetPt",        "TTH3l", "", "", 28, 20., 300.);
    theHistoManager->addHisto("JetEta",       "TTH3l", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("METpx",        "TTH3l", "", "",14,-140.,140.);
    theHistoManager->addHisto("METpy",        "TTH3l", "", "",14,-140.,140.);
    theHistoManager->addHisto("MET"  ,        "TTH3l", "", "",10,0.,400.);
    theHistoManager->addHisto("METphi",       "TTH3l", "", "",16,-3.2,3.2);
    theHistoManager->addHisto("METsum",       "TTH3l", "", "",15,0.,3000.);
    theHistoManager->addHisto("MHT",          "TTH3l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("MET LD",       "TTH3l", "", "", 10, 0., 2.);
    theHistoManager->addHisto("lep1 Mt",      "TTH3l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1 dRmin",   "TTH3l", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep2 dRmin",   "TTH3l", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep Eta max",  "TTH3l", "", "", 12, 0., 2.4);
    theHistoManager->addHisto("jet dR av",    "TTH3l", "", "", 10, 0., 5.);
    theHistoManager->addHisto("LepId",        "TTH3l", "", "", 4, -0.5, 3.5);
    theHistoManager->addHisto("HT",           "TTH3l", "", "",10, 0., 1000.);
    theHistoManager->addHisto("Mll min",      "TTH3l", "", "",10, 0., 200.);
    theHistoManager->addHisto("best MZ",      "TTH3l", "", "",12,60., 120.);

    // WZ
    theHistoManager->addHisto("nLep",         "WZ", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("nFake",        "WZ", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("lep1Pt",       "WZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep2Pt",       "WZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep3Pt",       "WZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1Eta",      "WZ", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep2Eta",      "WZ", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep3Eta",      "WZ", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lepQ",         "WZ", "", "", 9, -4.5, 4.5);
    theHistoManager->addHisto("nJets",        "WZ", "", "", 8, -0.5, 7.5);
    theHistoManager->addHisto("nLooseB",      "WZ", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("nMediumB",     "WZ", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("JetPt",        "WZ", "", "", 28, 20., 300.);
    theHistoManager->addHisto("JetEta",       "WZ", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("METpx",        "WZ", "", "",14,-140.,140.);
    theHistoManager->addHisto("METpy",        "WZ", "", "",14,-140.,140.);
    theHistoManager->addHisto("MET"  ,        "WZ", "", "",10,0.,400.);
    theHistoManager->addHisto("METphi",       "WZ", "", "",16,-3.2,3.2);
    theHistoManager->addHisto("METsum",       "WZ", "", "",15,0.,3000.);
    theHistoManager->addHisto("MHT",          "WZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("MET LD",       "WZ", "", "", 10, 0., 2.);
    theHistoManager->addHisto("lep1 Mt",      "WZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1 dRmin",   "WZ", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep2 dRmin",   "WZ", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep Eta max",  "WZ", "", "", 12, 0., 2.4);
    theHistoManager->addHisto("jet dR av",    "WZ", "", "", 10, 0., 5.);
    theHistoManager->addHisto("W Mt",         "WZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("Z Pt",         "WZ", "", "", 20,  0., 500.);
    theHistoManager->addHisto("Pt Sum(l)",    "WZ", "", "", 20,  0., 500.);
    theHistoManager->addHisto("inv.mass(l)",  "WZ", "", "", 20,  0., 600.);
    theHistoManager->addHisto("LepId",        "WZ", "", "", 4, -0.5, 3.5);
    theHistoManager->addHisto("HT",           "WZ", "", "",10, 0., 1000.);
    theHistoManager->addHisto("Mll min",      "WZ", "", "",10, 0., 200.);
    theHistoManager->addHisto("best MZ",      "WZ", "", "",12,60., 120.);

    // WZrelaxed
    theHistoManager->addHisto("nLep",         "WZrelaxed", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("nFake",        "WZrelaxed", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("lep1Pt",       "WZrelaxed", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep2Pt",       "WZrelaxed", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep3Pt",       "WZrelaxed", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1Eta",      "WZrelaxed", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep2Eta",      "WZrelaxed", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep3Eta",      "WZrelaxed", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lepQ",         "WZrelaxed", "", "", 9, -4.5, 4.5);
    theHistoManager->addHisto("nJets",        "WZrelaxed", "", "", 8, -0.5, 7.5);
    theHistoManager->addHisto("nLooseB",      "WZrelaxed", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("nMediumB",     "WZrelaxed", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("JetPt",        "WZrelaxed", "", "", 28, 20., 300.);
    theHistoManager->addHisto("JetEta",       "WZrelaxed", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("METpx",        "WZrelaxed", "", "",14,-140.,140.);
    theHistoManager->addHisto("METpy",        "WZrelaxed", "", "",14,-140.,140.);
    theHistoManager->addHisto("MET"  ,        "WZrelaxed", "", "",10,0.,400.);
    theHistoManager->addHisto("METphi",       "WZrelaxed", "", "",16,-3.2,3.2);
    theHistoManager->addHisto("METsum",       "WZrelaxed", "", "",15,0.,3000.);
    theHistoManager->addHisto("MHT",          "WZrelaxed", "", "", 10, 0., 200.);
    theHistoManager->addHisto("MET LD",       "WZrelaxed", "", "", 10, 0., 2.);
    theHistoManager->addHisto("lep1 Mt",      "WZrelaxed", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1 dRmin",   "WZrelaxed", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep2 dRmin",   "WZrelaxed", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep Eta max",  "WZrelaxed", "", "", 12, 0., 2.4);
    theHistoManager->addHisto("jet dR av",    "WZrelaxed", "", "", 10, 0., 5.);
    theHistoManager->addHisto("W Mt",         "WZrelaxed", "", "", 10, 0., 200.);
    theHistoManager->addHisto("Z Pt",         "WZrelaxed", "", "", 20,  0., 500.);
    theHistoManager->addHisto("Pt Sum(l)",    "WZrelaxed", "", "", 20,  0., 500.);
    theHistoManager->addHisto("inv.mass(l)",  "WZrelaxed", "", "", 20,  0., 600.);
    theHistoManager->addHisto("LepId",        "WZrelaxed", "", "", 4, -0.5, 3.5);
    theHistoManager->addHisto("HT",           "WZrelaxed", "", "",10, 0., 1000.);
    theHistoManager->addHisto("Mll min",      "WZrelaxed", "", "",10, 0., 200.);
    theHistoManager->addHisto("best MZ",      "WZrelaxed", "", "",12,60., 120.);

    // TTZ
    theHistoManager->addHisto("nLep",         "TTZ", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("nFake",        "TTZ", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("lep1Pt",       "TTZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep2Pt",       "TTZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep3Pt",       "TTZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep4Pt",       "TTZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1Eta",      "TTZ", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep2Eta",      "TTZ", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep3Eta",      "TTZ", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep4Eta",      "TTZ", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lepQ",         "TTZ", "", "", 9, -4.5, 4.5);
    theHistoManager->addHisto("nJets",        "TTZ", "", "", 8, -0.5, 7.5);
    theHistoManager->addHisto("nLooseB",      "TTZ", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("nMediumB",     "TTZ", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("JetPt",        "TTZ", "", "", 28, 20., 300.);
    theHistoManager->addHisto("JetEta",       "TTZ", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("METpx",        "TTZ", "", "",14,-140.,140.);
    theHistoManager->addHisto("METpy",        "TTZ", "", "",14,-140.,140.);
    theHistoManager->addHisto("MET"  ,        "TTZ", "", "",10,0.,400.);
    theHistoManager->addHisto("METphi",       "TTZ", "", "",16,-3.2,3.2);
    theHistoManager->addHisto("METsum",       "TTZ", "", "",15,0.,3000.);
    theHistoManager->addHisto("MHT",          "TTZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("MET LD",       "TTZ", "", "", 10, 0., 2.);
    theHistoManager->addHisto("lep1 Mt",      "TTZ", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1 dRmin",   "TTZ", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep2 dRmin",   "TTZ", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep Eta max",  "TTZ", "", "", 12, 0., 2.4);
    theHistoManager->addHisto("jet dR av",    "TTZ", "", "", 10, 0., 5.);
    theHistoManager->addHisto("LepId",        "TTZ", "", "", 4, -0.5, 3.5);
    theHistoManager->addHisto("HT",           "TTZ", "", "",10, 0., 1000.);
    theHistoManager->addHisto("Mll min",      "TTZ", "", "",10, 0., 200.);
    theHistoManager->addHisto("best MZ",      "TTZ", "", "",12,60., 120.);

    // Zl
    theHistoManager->addHisto("nLep",         "Zl", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("nFake",        "Zl", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("lep1Pt",       "Zl", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep2Pt",       "Zl", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep3Pt",       "Zl", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1Eta",      "Zl", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep2Eta",      "Zl", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep3Eta",      "Zl", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lepQ",         "Zl", "", "", 9, -4.5, 4.5);
    theHistoManager->addHisto("nJets",        "Zl", "", "", 8, -0.5, 7.5);
    theHistoManager->addHisto("nLooseB",      "Zl", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("nMediumB",     "Zl", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("JetPt",        "Zl", "", "", 28, 20., 300.);
    theHistoManager->addHisto("JetEta",       "Zl", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("METpx",        "Zl", "", "",14,-140.,140.);
    theHistoManager->addHisto("METpy",        "Zl", "", "",14,-140.,140.);
    theHistoManager->addHisto("MET"  ,        "Zl", "", "",10,0.,400.);
    theHistoManager->addHisto("METphi",       "Zl", "", "",16,-3.2,3.2);
    theHistoManager->addHisto("METsum",       "Zl", "", "",15,0.,3000.);
    theHistoManager->addHisto("MHT",          "Zl", "", "", 10, 0., 200.);
    theHistoManager->addHisto("MET LD",       "Zl", "", "", 10, 0., 2.);
    theHistoManager->addHisto("lep1 Mt",      "Zl", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1 dRmin",   "Zl", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep2 dRmin",   "Zl", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep Eta max",  "Zl", "", "", 12, 0., 2.4);
    theHistoManager->addHisto("jet dR av",    "Zl", "", "", 10, 0., 5.);
    theHistoManager->addHisto("LepId",        "Zl", "", "", 4, -0.5, 3.5);
    theHistoManager->addHisto("HT",           "Zl", "", "",10, 0., 1000.);
    theHistoManager->addHisto("Mll min",      "Zl", "", "",10, 0., 200.);
    theHistoManager->addHisto("best MZ",      "Zl", "", "",12,60., 120.);

    // TTH2l
    theHistoManager->addHisto("nLep",         "TTH2l", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("nFake",        "TTH2l", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("lep1Pt",       "TTH2l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep2Pt",       "TTH2l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep3Pt",       "TTH2l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1Eta",      "TTH2l", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep2Eta",      "TTH2l", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep3Eta",      "TTH2l", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lepQ",         "TTH2l", "", "", 9, -4.5, 4.5);
    theHistoManager->addHisto("nJets",        "TTH2l", "", "", 8, -0.5, 7.5);
    theHistoManager->addHisto("nLooseB",      "TTH2l", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("nMediumB",     "TTH2l", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("JetPt",        "TTH2l", "", "", 28, 20., 300.);
    theHistoManager->addHisto("JetEta",       "TTH2l", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("METpx",        "TTH2l", "", "",14,-140.,140.);
    theHistoManager->addHisto("METpy",        "TTH2l", "", "",14,-140.,140.);
    theHistoManager->addHisto("MET"  ,        "TTH2l", "", "",10,0.,400.);
    theHistoManager->addHisto("METphi",       "TTH2l", "", "",16,-3.2,3.2);
    theHistoManager->addHisto("METsum",       "TTH2l", "", "",15,0.,3000.);
    theHistoManager->addHisto("MHT",          "TTH2l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("MET LD",       "TTH2l", "", "", 10, 0., 2.);
    theHistoManager->addHisto("lep1 Mt",      "TTH2l", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1 dRmin",   "TTH2l", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep2 dRmin",   "TTH2l", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep Eta max",  "TTH2l", "", "", 12, 0., 2.4);
    theHistoManager->addHisto("jet dR av",    "TTH2l", "", "", 10, 0., 5.);
    theHistoManager->addHisto("LepId",        "TTH2l", "", "", 4, -0.5, 3.5);
    theHistoManager->addHisto("HT",           "TTH2l", "", "",10, 0., 1000.);
    theHistoManager->addHisto("Mll min",      "TTH2l", "", "",10, 0., 200.);
    theHistoManager->addHisto("best MZ",      "TTH2l", "", "",12,60., 120.);

    // TTdilep
    theHistoManager->addHisto("nLep",         "TTemu", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("nFake",        "TTemu", "", "", 10, -0.5, 9.5);
    theHistoManager->addHisto("lep1Pt",       "TTemu", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep2Pt",       "TTemu", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep3Pt",       "TTemu", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1Eta",      "TTemu", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep2Eta",      "TTemu", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lep3Eta",      "TTemu", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("lepQ",         "TTemu", "", "", 9, -4.5, 4.5);
    theHistoManager->addHisto("nJets",        "TTemu", "", "", 8, -0.5, 7.5);
    theHistoManager->addHisto("nLooseB",      "TTemu", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("nMediumB",     "TTemu", "", "", 6, -0.5, 5.5);
    theHistoManager->addHisto("JetPt",        "TTemu", "", "", 28, 20., 300.);
    theHistoManager->addHisto("JetEta",       "TTemu", "", "", 12, -2.4, 2.4);
    theHistoManager->addHisto("METpx",        "TTemu", "", "",14,-140.,140.);
    theHistoManager->addHisto("METpy",        "TTemu", "", "",14,-140.,140.);
    theHistoManager->addHisto("MET"  ,        "TTemu", "", "",10,0.,400.);
    theHistoManager->addHisto("METphi",       "TTemu", "", "",16,-3.2,3.2);
    theHistoManager->addHisto("METsum",       "TTemu", "", "",15,0.,3000.);
    theHistoManager->addHisto("MHT",          "TTemu", "", "", 10, 0., 200.);
    theHistoManager->addHisto("MET LD",       "TTemu", "", "", 10, 0., 2.);
    theHistoManager->addHisto("lep1 Mt",      "TTemu", "", "", 10, 0., 200.);
    theHistoManager->addHisto("lep1 dRmin",   "TTemu", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep2 dRmin",   "TTemu", "", "", 10, 0., 4.);
    theHistoManager->addHisto("lep Eta max",  "TTemu", "", "", 12, 0., 2.4);
    theHistoManager->addHisto("jet dR av",    "TTemu", "", "", 10, 0., 5.);
    theHistoManager->addHisto("LepId",        "TTemu", "", "", 4, -0.5, 3.5);
    theHistoManager->addHisto("HT",           "TTemu", "", "",10, 0., 1000.);
    theHistoManager->addHisto("Mll min",      "TTemu", "", "",10, 0., 200.);
    theHistoManager->addHisto("best MZ",      "TTemu", "", "",12,60., 120.);

    // MVA
    theHistoManager->addHisto("Signal_2lss_TT_MVA",                       "FinalCut", "ttH2lss",   "",  20,   -1,     1);
    theHistoManager->addHisto("Signal_2lss_TTV_MVA",                      "FinalCut", "ttH2lss",   "",  20,   -1,     1);
    theHistoManager->addHisto("Signal_3l_TT_MVA",                         "FinalCut",   "ttH3l",   "",  20,   -1,     1);
    theHistoManager->addHisto("Signal_3l_TTV_MVA",                        "FinalCut",   "ttH3l",   "",  20,   -1,     1);

    std::string inputFileHF = "/home-pbs/xcoubez/ttHAnalysis_Git/ttHAnalysis_76X_Moriond_ALLCORRECTIONS/IPHCNtuple/NtupleAnalyzer/src/weight/csv_rwt_fit_hf_76x_2016_02_08.root";
    std::string inputFileLF = "/home-pbs/xcoubez/ttHAnalysis_Git/ttHAnalysis_76X_Moriond_ALLCORRECTIONS/IPHCNtuple/NtupleAnalyzer/src/weight/csv_rwt_fit_lf_76x_2016_02_08.root";

    TFile* f_CSVwgt_HF = new TFile ((inputFileHF).c_str());
    TFile* f_CSVwgt_LF = new TFile ((inputFileLF).c_str());

    fillCSVhistos(f_CSVwgt_HF, f_CSVwgt_LF);
}


void TTbarHiggsMultileptonAnalysis::writeHistograms()
{  
    outputfile->cd();

    std::vector<TH1F*> the1DHisto =  theHistoManager->getHisto1D_list();
    std::vector<TH2F*> the2DHisto =  theHistoManager->getHisto2D_list();

    for(unsigned int i=0; i<the1DHisto.size(); i++)  the1DHisto[i]->Write();
    for(unsigned int i=0; i<the2DHisto.size(); i++)  the2DHisto[i]->Write();

    tOutput->Write();
    outputfile->Close();
}


void TTbarHiggsMultileptonAnalysis::Init(TChain *tree)
{
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;

    vEvent    = new std::vector<Event>();
    vElectron = new std::vector<Electron>();
    vMuon     = new std::vector<Muon>();
    vTau      = new std::vector<Tau>();
    vJet      = new std::vector<Jet>();
    vTruth    = new std::vector<Truth>();

    fChain->SetBranchAddress("Event",    &vEvent   );
    fChain->SetBranchAddress("Electron", &vElectron);
    fChain->SetBranchAddress("Muon",     &vMuon    );
    fChain->SetBranchAddress("Tau",      &vTau     );
    fChain->SetBranchAddress("Jet",      &vJet     );
    fChain->SetBranchAddress("Truth",    &vTruth   );

    Load_MVA();
}


void TTbarHiggsMultileptonAnalysis::Loop()
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntries();
    int nentries_max = nentries;
    if ( _nmax != -1 && _nmax < nentries ) nentries_max = _nmax;

    std::cout << "Number of input events = " << nentries << std::endl;
    std::cout << "Number of processed events = " << nentries_max << std::endl;

    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<nentries_max;jentry++) 
    {

        Long64_t ientry = fChain->LoadTree(jentry);
        if (ientry < 0) break;

        if(jentry%10000 == 0) std::cout << "number of processed events " << jentry << std::endl;
        //std::cout << "number of processed events " << jentry << std::endl;

        //if(jentry > 100000) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        //
        int pvn = vEvent->at(0).pv_n();
        theHistoManager->fillHisto("NumberOfPrimaryVertex", "noSel", "", "",  pvn, 1);

	mc_event = vEvent->at(0).id();

        if ( !_isdata )
        {
            weight = _lumi*_xsec/_nowe;
            mc_weight = vEvent->at(0).mc_weight();
            //weight_PV = _h_PV->GetBinContent(pvn);
            weight = weight * mc_weight; //*weight_PV;

	    weight_scale_muF0p5 = vEvent->at(0).weight_scale_muF0p5();
	    weight_scale_muF2 = vEvent->at(0).weight_scale_muF2();
	    weight_scale_muR0p5 = vEvent->at(0).weight_scale_muR0p5();
	    weight_scale_muR2 = vEvent->at(0).weight_scale_muR2();
        }
        else 
        {
            weight    = 1.;
            mc_weight = 1.;
            weight_PV = 1.; 

            /////////////////////////////////////
            // remove double counting in data

            // MuonEG
            // HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v
            // HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v
            // HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v
            // HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v
            // 
            // DoubleMuon
            // HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v
            // HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v
            // HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v
            // HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v
            // HLT_TripleMu_12_10_5_v
            // 
            // DoubleEG
            // HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
            // HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v
            // HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v
            // 
            // SingleMu
            // HLT_IsoMu20_v
            // HLT_IsoTkMu20_v
            // 
            // SingleEG
            // HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v
            // HLT_Ele23_WPLoose_Gsf_v

            // détricotage
            bool TRIGm = false, TRIGe = false, TRIGeData = false, TRIGmTk = false; 
            bool TRIGee = false, TRIGmm = false, TRIGme = false, TRIGem = false, TRIGmmTk = false;
            bool TRIGeee = false, TRIGmme = false, TRIGeem = false, TRIGmmm = false;

            int a = ( vEvent->at(0).ev_trigger_pass_byname_1() )%10;
            int b = ((vEvent->at(0).ev_trigger_pass_byname_1() -a)/10)%10;
            int c = ((vEvent->at(0).ev_trigger_pass_byname_1() -a-10*b)/100)%10;
            int d = ((vEvent->at(0).ev_trigger_pass_byname_1() -a-10*b-100*c)/1000)%10;
            int e =   vEvent->at(0).ev_trigger_pass_byname_1();

            if (a==1 || a==3 || a==6 || a==8) TRIGeee  = true; // HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*
            if (a==2 || a==3 || a==7 || a==8) TRIGme   = true; // HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v* 
            if (a>=5 )                        TRIGem   = true; // HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v* 
            if (b==1 || b==3 || b==6 || b==8) TRIGee   = true; // HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v* 
            if (b==2 || b==3 || b==7 || b==8) TRIGmm   = true; // HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v* 
            if (b>=5 )                        TRIGmmTk = true; // HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v* 
            if (c==1 || c==3 || c==6 || c==8) TRIGm    = true; // HLT_IsoMu20_v* 
            if (c==2 || c==3 || c==7 || c==8) TRIGmTk  = true; // HLT_IsoTkMu20_v* 
            if (c>=5 )                        TRIGe    = true; // HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v* 
            if (d==1 || d==3 || d==6 || d==8) TRIGeData= true; // HLT_Ele23_WPLoose_Gsf_v* 
            if (d==2 || d==3 || d==7 || d==8) TRIGmme  = true; // HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v*  but doesn't exist ?
            if (d>=5 )                        TRIGeem  = true; // HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v* but doesn't exist ?
            if (e>=10000)                     TRIGmmm  = true; // HLT_TripleMu_12_10_5_v*             but doesn't exist ?

            bool E = false, M = false, EE = false, MM = false, EM = false;
            if ( TRIGme || TRIGem || TRIGeem || TRIGmme ) EM = true;
            if ( TRIGmm || TRIGmmTk || TRIGmmm )          MM = true;
            if ( TRIGee || TRIGeee )	              EE = true;
            if ( TRIGm  || TRIGmTk )                      M  = true;
            if ( TRIGe  || TRIGeData )                    E  = true;

            // new code from Xavier (with next Ntuple production)
            //         EM = ( vEvent->at(0).ev_pass_eem() || vEvent->at(0).ev_pass_em()     || vEvent->at(0).ev_pass_mme()    || vEvent->at(0).ev_pass_me() );
            //         MM = ( vEvent->at(0).ev_pass_mm()  || vEvent->at(0).ev_pass_mmTk()   || vEvent->at(0).ev_pass_mmnoDz() || vEvent->at(0).ev_pass_mmTknoDz() || vEvent->at(0).ev_pass_mmm() ) ;
            //         EE = ( vEvent->at(0).ev_pass_ee()  || vEvent->at(0).ev_pass_eenoDz() || vEvent->at(0).ev_pass_eee() );
            //         M  = ( vEvent->at(0).ev_pass_m()   || vEvent->at(0).ev_pass_mTk() );
            //         E  = ( vEvent->at(0).ev_pass_e()   || vEvent->at(0).ev_pass_eData() );

            // old code from Xavier
            //         int tab[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            //         int a = vEvent->at(0).ev_trigger_pass_byname_1();
            //         int n = 0, size = 0;
            // 
            //         do{
            //             n = a%10;
            //             a = a/10;
            //             tab[size] = n;
            //             size = size + 1;
            //         }while(a!=0);
            // 
            //         bool E = false, M = false, EE = false, MM = false, EM = false;
            //         int result_trigger = 0;
            //         if (size > 2) {if ( tab[3] == 1                ) E  = true;}
            //         if (size > 1) {if ( tab[2] == 1 || tab[2] == 2 ) M  = true;}
            //         if (size > 0) {if ( tab[1] == 1                ) EE = true;}
            //         if (size > 1) {if ( tab[1] == 2 || tab[1] == 5 ) MM = true;}
            //         if ( tab[0] == 2 || tab[0] == 5 )                EM = true;

            bool emdataset = _sampleName.Contains("MuonEG");
            bool mmdataset = _sampleName.Contains("DoubleM");
            bool eedataset = _sampleName.Contains("DoubleE");
            bool mdataset  = _sampleName.Contains("SingleM");
            bool edataset  = _sampleName.Contains("SingleE");

            int result_trigger = 0;
            if ( EM  &&                               (emdataset) ) result_trigger = 1;
            if ( !EM && MM  &&                        (mmdataset) ) result_trigger = 1;
            if ( !EM && !MM && EE  &&                 (eedataset) ) result_trigger = 1;
            if ( !EM && !MM && !EE && M  &&           (mdataset ) ) result_trigger = 1;
            if ( !EM && !MM && !EE && !M && E &&      (edataset ) ) result_trigger = 1;

            if(result_trigger == 1)
            {
                weight = 1;
                //std::cout << " EM: " << EM 
                //          << " MM: " << MM 
                //          << " EE: " << EE
                //          << " M:  " << M
                //          << " E:  " << E                    << std::endl;
                //std::cout << " sampleName: " << _sampleName   << std::endl;
                //std::cout << " weight:     " << weight       << std::endl; 
            }
            else
            {
                weight = 0;
                //std::cout << " EM: " << EM 
                //          << " MM: " << MM  
                //          << " EE: " << EE
                //          << " M:  " << M
                //          << " E:  " << E                    << std::endl;
                //std::cout << " sampleName: " << _sampleName   << std::endl;
                //std::cout << " weight:     " << weight       << std::endl;  
            }
        }


        // ###########################################################
        // #  _       _ _   _       _ _           _   _              #
        // # (_)_ __ (_) |_(_) __ _| (_)___  __ _| |_(_) ___  _ __   #
        // # | | '_ \| | __| |/ _` | | / __|/ _` | __| |/ _ \| '_ \  #
        // # | | | | | | |_| | (_| | | \__ \ (_| | |_| | (_) | | | | #
        // # |_|_| |_|_|\__|_|\__,_|_|_|___/\__,_|\__|_|\___/|_| |_| #
        // #                                                         #
        // ###########################################################

        vLeptons.clear();
        vSelectedMuons.clear();
        vSelectedElectrons.clear();
        vSelectedLeptons.clear();
        vFakeMuons.clear();
        vFakeElectrons.clear();
        vFakeLeptons.clear();	
        vSelectedNonBTagJets.clear();
        vSelectedBTagJets.clear();
        vSelectedMediumBTagJets.clear();
        vSelectedJets.clear();

        is_2lss_TTH_SR    = false;
        is_2lss_JM_SB     = false;
        is_2lss_LepMVA_SB = false;
        is_emu_TT_CR      = false;

        is_3l_TTH_SR      = false;
        is_3l_WZ_CR       = false; 
        is_3l_WZrel_CR    = false;
        is_3l_TTZ_CR      = false;
        is_Zl_CR          = false;

        // ######################################
        // #  _        _                        #
        // # | |_ _ __(_) __ _  __ _  ___ _ __  #
        // # | __| '__| |/ _` |/ _` |/ _ \ '__| #
        // # | |_| |  | | (_| | (_| |  __/ |    #
        // #  \__|_|  |_|\__, |\__, |\___|_|    #
        // #             |___/ |___/            #
        // #                                    #
        // ######################################

        is_trigger = false;
        if ( vEvent->at(0).ev_trigger_pass_byname_1() >= 1 ) is_trigger = true;

        // #####################################
        // #  _ __ ___  _   _  ___  _ __  ___  #
        // # | '_ ` _ \| | | |/ _ \| '_ \/ __| #
        // # | | | | | | |_| | (_) | | | \__ \ #
        // # |_| |_| |_|\__,_|\___/|_| |_|___/ #
        // #                                   #
        // #####################################

        for(unsigned int imuon=0; imuon < vMuon->size() ; imuon++)
        {   
            Lepton l; l.setLepton(&vMuon->at(imuon),imuon,0);

            if ( vMuon->at(imuon).isTightTTH() )
            {
                vSelectedMuons.push_back(vMuon->at(imuon));
                vSelectedLeptons.push_back(l);
            }
            else if ( vMuon->at(imuon).isFakeableTTH() )
            {
                vFakeMuons.push_back(vMuon->at(imuon));
                vFakeLeptons.push_back(l);
            }

            vLeptons.push_back(l);

            theHistoManager->fillHisto("MuonPt",                            "noSel",        "",   "",  vMuon->at(imuon).pt(),             weight);
            theHistoManager->fillHisto("MuonEta",                           "noSel",        "",   "",  vMuon->at(imuon).eta(),            weight);
            theHistoManager->fillHisto("MuonMVA",                           "noSel",        "",   "",  vMuon->at(imuon).lepMVA(),         weight);
        }     

        // ##############################################
        // #       _           _                        #
        // #   ___| | ___  ___| |_ _ __ ___  _ __  ___  #
        // #  / _ \ |/ _ \/ __| __| '__/ _ \| '_ \/ __| #
        // # |  __/ |  __/ (__| |_| | | (_) | | | \__ \ #
        // #  \___|_|\___|\___|\__|_|  \___/|_| |_|___/ #
        // #                                            #
        // ##############################################

        for(unsigned int ielectron=0; ielectron < vElectron->size() ; ielectron++)
        {   
            Lepton l; l.setLepton(&vElectron->at(ielectron),ielectron,1);

            if ( vElectron->at(ielectron).isTightTTH() )
            {
                vSelectedElectrons.push_back(vElectron->at(ielectron));	     
                vSelectedLeptons.push_back(l);
            }
            else if ( vElectron->at(ielectron).isFakeableTTH() )
            {
                vFakeElectrons.push_back(vElectron->at(ielectron));	     
                vFakeLeptons.push_back(l);
            }

            vLeptons.push_back(l);

            theHistoManager->fillHisto("ElectronPt",                        "noSel",        "",   "",  vElectron->at(ielectron).pt(),     weight);
            theHistoManager->fillHisto("ElectronEta",                       "noSel",        "",   "",  vElectron->at(ielectron).eta(),    weight);
            theHistoManager->fillHisto("ElectronMVA",                       "noSel",        "",   "",  vElectron->at(ielectron).lepMVA(), weight);
        }  

        // ########################
        // #  _                   #
        // # | |_ __ _ _   _ ___  #
        // # | __/ _` | | | / __| #
        // # | || (_| | |_| \__ \ #
        // #  \__\__,_|\__,_|___/ #
        // #                      #
        // ########################            

        for(unsigned int itau=0; itau < vTau->size() ; itau++)
        {
            Lepton l; l.setLepton(&vTau->at(itau),itau,1);

            vSelectedTaus.push_back(vTau->at(itau));
            //vSelectedLeptons.push_back(l);

            //vLeptons.push_back(l);

            theHistoManager->fillHisto("TauPt",                             "noSel",        "",   "",  vTau->at(itau).pt(),               weight);
            theHistoManager->fillHisto("TauEta",                            "noSel",        "",   "",  vTau->at(itau).eta(),              weight);
        }

        // #############################################
        // #                _           _              #
        // #   ___  _ __ __| | ___ _ __(_)_ __   __ _  #
        // #  / _ \| '__/ _` |/ _ \ '__| | '_ \ / _` | #
        // # | (_) | | | (_| |  __/ |  | | | | | (_| | #
        // #  \___/|_|  \__,_|\___|_|  |_|_| |_|\__, | #
        // #                                    |___/  #
        // #                                           #
        // #############################################

        std::sort(vLeptons.begin(), vLeptons.end(), SortingLeptonPt);
        std::sort(vSelectedLeptons.begin(), vSelectedLeptons.end(), SortingLeptonPt);
        std::sort(vFakeLeptons.begin(), vFakeLeptons.end(), SortingLeptonPt);

        if ( vLeptons.size() >= 2) {
            for (unsigned int ilep = 0; ilep<vLeptons.size()-1; ilep++) {
                if ( vLeptons.at(ilep).pt() < vLeptons.at(ilep+1).pt() ) {
                    std::cout << "Run Event " << vEvent->at(0).run() <<  " " << vEvent->at(0).id() 
                        << " nlep/tight/fake " << vLeptons.size() <<  "/" << vSelectedLeptons.size() <<  "/" <<  vFakeLeptons.size() << std::endl; 
                    std::cout << "  all pt[" << ilep << "]: " << vLeptons.at(ilep).pt() 
                        << " pt[" << ilep+1 << "]: " << vLeptons.at(ilep+1).pt() << std::endl;
                }
            }
        }

        if ( vSelectedLeptons.size() >= 2) {
            for (unsigned int ilep = 0; ilep<vSelectedLeptons.size()-1; ilep++) {
                if ( vSelectedLeptons.at(ilep).pt() < vSelectedLeptons.at(ilep+1).pt() ) {
                    std::cout << "Run Event " << vEvent->at(0).run() <<  " " << vEvent->at(0).id() 
                        << " nlep/tight/fake " << vLeptons.size() <<  "/" << vSelectedLeptons.size() <<  "/" <<  vFakeLeptons.size() << std::endl; 
                    std::cout << "tight pt[" << ilep << "]: " << vSelectedLeptons.at(ilep).pt() 
                        << " pt[" << ilep+1 << "]: " << vSelectedLeptons.at(ilep+1).pt() << std::endl;
                }
            }
        }

        if ( vFakeLeptons.size() >= 2) {
            for (unsigned int ilep = 0; ilep<vFakeLeptons.size()-1; ilep++) {
                if ( vFakeLeptons.at(ilep).pt() < vFakeLeptons.at(ilep+1).pt() ) {
                    std::cout << "Run Event " << vEvent->at(0).run() <<  " " << vEvent->at(0).id() 
                        << " nlep/tight/fake " << vLeptons.size() <<  "/" << vSelectedLeptons.size() <<  "/" <<  vFakeLeptons.size() << std::endl; 
                    std::cout << " fake pt[" << ilep << "]: " << vFakeLeptons.at(ilep).pt() 
                        << " pt[" << ilep+1 << "]: " << vFakeLeptons.at(ilep+1).pt() << std::endl;
                }
            }
        }

        // ################################
        // #                              #
        // #  _           _      _        #
        // # | |__       (_) ___| |_ ___  #
        // # | '_ \ _____| |/ _ \ __/ __| #
        // # | |_) |_____| |  __/ |_\__ \ #
        // # |_.__/     _/ |\___|\__|___/ #
        // #           |__/               #
        // #                              #
        // ################################

        nLooseBJets  = 0;
        nMediumBJets = 0;

        for(unsigned int ijet=0; ijet < vJet->size() ; ijet++)
        {

            // version for 74x as in the AN...:
            if( vJet->at(ijet).CSVv2() > 0.423 ) nLooseBJets++;
            if( vJet->at(ijet).CSVv2() > 0.814 ) nMediumBJets++;

            if(vJet->at(ijet).CSVv2() >= 0.423 ) vSelectedBTagJets.push_back(vJet->at(ijet));
            else                                 vSelectedNonBTagJets.push_back(vJet->at(ijet));
            if(vJet->at(ijet).CSVv2() >= 0.814 ) vSelectedMediumBTagJets.push_back(vJet->at(ijet));

            // to be updated for 76x ???
            //             if( vJet->at(ijet).CSVv2() > 0.460 ) nLooseBJets++;
            //             if( vJet->at(ijet).CSVv2() > 0.800 ) nMediumBJets++;
            // 
            //             if(vJet->at(ijet).CSVv2() >= 0.460 ) vSelectedBTagJets.push_back(vJet->at(ijet));
            //             else                                 vSelectedNonBTagJets.push_back(vJet->at(ijet));
            //             if(vJet->at(ijet).CSVv2() >= 0.800 ) vSelectedMediumBTagJets.push_back(vJet->at(ijet));

            vSelectedJets.push_back(vJet->at(ijet));

            theHistoManager->fillHisto("JetPt",                             "noSel",        "",   "",  vJet->at(ijet).pt(),               weight);
            theHistoManager->fillHisto("JetEta",                            "noSel",        "",   "",  vJet->at(ijet).eta(),              weight);
            theHistoManager->fillHisto("JetCSVv2",                          "noSel",        "",   "",  vJet->at(ijet).CSVv2(),            weight);
        }

        theHistoManager->fillHisto("MET",                               "noSel",        "",   "",  vEvent->at(0).metpt(),            weight );

        theHistoManager->fillHisto("CutFlow",                        "noSel", "", "", 1, 1);

        if( vMuon->size()+vElectron->size() == 3 )     
        {
            theHistoManager->fillHisto("CutFlow",             "ThreePreselected", "", "", 1, 1);
        }

        // ################################################################################
        // #  ____  ____    ____  ____ _____                   _       _     _            #
        // # |___ \|  _ \  | __ )|  _ \_   _| __   ____ _ _ __(_) __ _| |__ | | ___  ___  #
        // #   __) | | | | |  _ \| | | || |   \ \ / / _` | '__| |/ _` | '_ \| |/ _ \/ __| #
        // #  / __/| |_| | | |_) | |_| || |    \ V / (_| | |  | | (_| | |_) | |  __/\__ \ #
        // # |_____|____/  |____/|____/ |_|     \_/ \__,_|_|  |_|\__,_|_.__/|_|\___||___/ #
        // #                                                                              #
        // ################################################################################

        max_Lep_eta     = 0. ;
        //numJets_float   = 0. ;
        nJet25_Recl     = 0 ;
        mindr_lep1_jet  = 0. ;
        mindr_lep2_jet  = 0. ;
        met             = 0. ;
        avg_dr_jet      = 0. ;
        MT_met_lep1     = 0. ;
        LepGood_conePt0 = 0. ;
        LepGood_conePt1 = 0. ;
        mhtJet25_Recl   = 0. ;  

        // ############################################
        // #           _           _   _              #
        // #  ___  ___| | ___  ___| |_(_) ___  _ __   #
        // # / __|/ _ \ |/ _ \/ __| __| |/ _ \| '_ \  #
        // # \__ \  __/ |  __/ (__| |_| | (_) | | | | #
        // # |___/\___|_|\___|\___|\__|_|\___/|_| |_| #
        // #                                          #
        // ############################################

        TwoLeptonsSameSignSelection_TTH2l(jentry);
        //TwoLeptonsSameSignSelection_LepMVA_sideband(jentry);
        //TwoLeptonsSameSignSelection_JetMultiplicity_sideband(jentry);
        DiLeptonSelection_TT_CR(jentry);

        ThreeLeptonSelection_TTH3l(jentry);
        ThreeLeptonSelection_CR_WZ(jentry);
        ThreeLeptonSelection_CR_WZrelaxed(jentry);
        ThreeLeptonSelection_CR_Zl(jentry);
        ThreeLeptonSelection_TTZ(jentry);

        //std::cout <<is_CR_TTl<<" "<< is_Zl_CR <<" " << is_CR_WZ<<" " << is_TTH3l<< std::endl;
        //if (is_TTH3l==true ) std::cout <<"is_TTH3l" << std::endl;
        if ( is_2lss_TTH_SR || is_3l_TTH_SR ) fillOutputTree();

        // ########################################################################################################
        // #  __  __           _              _       _     _     _     _   _  ____ ___        _          __  __  #
        // # |  \/  | __ _  __| |_      _____(_) __ _| |__ | |_  | |   | | | |/ ___/ _ \   ___| |_ _   _ / _|/ _| #
        // # | |\/| |/ _` |/ _` \ \ /\ / / _ \ |/ _` | '_ \| __| | |   | |_| | |  | | | | / __| __| | | | |_| |_  #
        // # | |  | | (_| | (_| |\ V  V /  __/ | (_| | | | | |_  | |___|  _  | |__| |_| | \__ \ |_| |_| |  _|  _| #
        // # |_|  |_|\__,_|\__,_| \_/\_/ \___|_|\__, |_| |_|\__| |_____|_| |_|\____\___/  |___/\__|\__,_|_| |_|   #
        // #                                    |___/                                                             #
        // #                                                                                                      #
        // ########################################################################################################

        if ( !_isdata && _printLHCO_MC && ThreeLeptonSelection_TTH3l_MC()) PrintLHCOforMadweight_MC(jentry);

        // Common Selection:
        if ( !(vLeptons.size() >= 2
                    && vSelectedJets.size() >= 2) ) continue;

        float MET = vEvent->at(0).metpt();
        float METphi = vEvent->at(0).metphi();
        float METx = MET * TMath::Cos(METphi);
        float METy = MET * TMath::Sin(METphi);
        float METsum = vEvent->at(0).metsumet();
        int nlepsel = 0;
        float jet_px = 0, jet_py = 0, lep_px = 0, lep_py = 0, MHT = 0, met_ld = 0, jetht = 0;
        for (int i=0; i<vSelectedLeptons.size(); i++) {
            lep_px += vSelectedLeptons.at(i).p4().Px();
            lep_py += vSelectedLeptons.at(i).p4().Py();
            if ( vSelectedLeptons.at(i).pt() > 10. ) nlepsel++;
        }
        if ( nlepsel >= 2 && vSelectedLeptons.at(0).pt() < 20. ) nlepsel = -1.;
        TLorentzVector jetp4;
        for (int ijet=0; ijet < vSelectedJets.size(); ijet++) {
            jetp4.SetPtEtaPhiE(vSelectedJets.at(ijet).pt(), vSelectedJets.at(ijet).eta(), vSelectedJets.at(ijet).phi(), vSelectedJets.at(ijet).E());
            jet_px += jetp4.Px();
            jet_py += jetp4.Py();
            jetht += vSelectedJets.at(ijet).pt();
            theHistoManager->fillHisto("JetPt",  "Trig", "", "", vJet->at(ijet).pt(), weight);
        }
        MHT = sqrt( (jet_px+lep_px)*(jet_px+lep_px) + (jet_py+lep_py)*(jet_py+lep_py) );
        met_ld = 0.00397 * MET + 0.00265 * MHT;

        float Mllmin = 1000., Mllbest = 1000., Deltabest = 1000.;
        theHistoManager->fillHisto("nLep",   "PreSel", "", "", nlepsel, weight);
        theHistoManager->fillHisto("nLep loose", "PreSel", "", "", vLeptons.size(), weight);

        if ( nlepsel >= 2 ) {
            theHistoManager->fillHisto("lep1Pt", "PreSel", "", "", vSelectedLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep1Eta","PreSel", "", "", vSelectedLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lep2Pt", "PreSel", "", "", vSelectedLeptons.at(1).pt(), weight);
            theHistoManager->fillHisto("lep2Eta","PreSel", "", "", vSelectedLeptons.at(1).eta(), weight);
            if ( nlepsel >= 3 ) theHistoManager->fillHisto("lep3Pt", "PreSel", "", "", vSelectedLeptons.at(2).pt(), weight);
            if ( nlepsel >= 3 ) theHistoManager->fillHisto("lep3Eta","PreSel", "", "", vSelectedLeptons.at(2).eta(), weight);
            float lepq = vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge();
            if ( nlepsel >= 3 ) lepq += vSelectedLeptons.at(2).charge();
            theHistoManager->fillHisto("lepQ", "PreSel", "", "", lepq, weight);
            for (int i=0; i<vSelectedLeptons.size()-1; i++) {
                for (int j=i+1; j<vSelectedLeptons.size(); j++) {
                    float mll = ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M();
                    if ( mll < Mllmin ) Mllmin = mll;
                    if ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id()) {
                        if ( fabs(mll - 91.188) < Deltabest ) {
                            Mllbest = mll;
                            Deltabest = fabs(mll - 91.188);
                        }
                        theHistoManager->fillHisto("Mll", "PreSel", "", "", mll, weight);
                        if ( nlepsel >= 3 ) theHistoManager->fillHisto("Mll", "PreSel 3l", "", "", mll, weight);
                        if ( nLooseBJets>=2 || nMediumBJets>=1 ) theHistoManager->fillHisto("Mll", "PreSel btag", "", "", mll, weight);
                        if ( (nLooseBJets>=2 || nMediumBJets>=1 ) && nlepsel >= 3 ) theHistoManager->fillHisto("Mll", "PreSel btag 3l", "", "", mll, weight);
                    }
                }
            }
        }

        theHistoManager->fillHisto("nJets",    "PreSel", "", "", vSelectedJets.size(), weight);
        theHistoManager->fillHisto("nLooseB",  "PreSel", "", "", nLooseBJets, weight);
        theHistoManager->fillHisto("nMediumB", "PreSel", "", "", nMediumBJets, weight);
        for (int ijet=0; ijet < vSelectedJets.size(); ijet++) {
            theHistoManager->fillHisto("JetPt",  "PreSel", "", "", vJet->at(ijet).pt(), weight);
            theHistoManager->fillHisto("JetEta", "PreSel", "", "", vJet->at(ijet).eta(), weight);
            theHistoManager->fillHisto("CSVv2",  "PreSel", "", "", vJet->at(ijet).CSVv2(), weight);
        }
        theHistoManager->fillHisto("METpx",  "PreSel", "", "", METx, weight);
        theHistoManager->fillHisto("METpy",  "PreSel", "", "", METy, weight);
        theHistoManager->fillHisto("MET"  ,  "PreSel", "", "", MET, weight);
        theHistoManager->fillHisto("METphi", "PreSel", "", "", METphi, weight);
        theHistoManager->fillHisto("METsum", "PreSel", "", "", METsum, weight);
        theHistoManager->fillHisto("MHT",    "PreSel", "", "", MHT, weight);
        theHistoManager->fillHisto("MET LD", "PreSel", "", "", met_ld, weight);

        int lepid2 = -1, lepid3 = -1;
        if ( nlepsel >= 2 ) {
            if ( abs(vSelectedLeptons.at(0).id()) == 11 
                    && abs(vSelectedLeptons.at(1).id()) == 11 ) lepid2 = 0; // ee
            if ( abs(vSelectedLeptons.at(0).id()) == 11 
                    && abs(vSelectedLeptons.at(1).id()) == 13 ) lepid2 = 1; // emu
            if ( abs(vSelectedLeptons.at(0).id()) == 13 
                    && abs(vSelectedLeptons.at(1).id()) == 11 ) lepid2 = 1;
            if ( abs(vSelectedLeptons.at(0).id()) == 13 
                    && abs(vSelectedLeptons.at(1).id()) == 13 ) lepid2 = 2; // mumu
        }
        if ( nlepsel >= 3 ) {
            if ( lepid2 == 0 && abs(vSelectedLeptons.at(2).id()) == 11 ) lepid3 = 0; // eee
            if ( lepid2 == 0 && abs(vSelectedLeptons.at(2).id()) == 13 ) lepid3 = 1; // eemu
            if ( lepid2 == 1 && abs(vSelectedLeptons.at(2).id()) == 11 ) lepid3 = 1;
            if ( lepid2 == 1 && abs(vSelectedLeptons.at(2).id()) == 13 ) lepid3 = 2; // emumu
            if ( lepid2 == 2 && abs(vSelectedLeptons.at(2).id()) == 11 ) lepid3 = 2;
            if ( lepid2 == 2 && abs(vSelectedLeptons.at(2).id()) == 13 ) lepid3 = 3; // mumumu
        }

        // ################################################################################
        // #  ____  ____    ____  ____ _____                   _       _     _            #
        // # |___ \|  _ \  | __ )|  _ \_   _| __   ____ _ _ __(_) __ _| |__ | | ___  ___  #
        // #   __) | | | | |  _ \| | | || |   \ \ / / _` | '__| |/ _` | '_ \| |/ _ \/ __| #
        // #  / __/| |_| | | |_) | |_| || |    \ V / (_| | |  | | (_| | |_) | |  __/\__ \ #
        // # |_____|____/  |____/|____/ |_|     \_/ \__,_|_|  |_|\__,_|_.__/|_|\___||___/ #
        // #                                                                              #
        // ################################################################################

        float dr;
        float lep1_mtw = -1., lep1_dr_min = 100., lep2_dr_min = 100., lep_eta_max = -1.;

        if (   vSelectedLeptons.size()     >= 2
                && vSelectedLeptons.at(0).pt() >  20 
                && vSelectedLeptons.at(1).pt() >  10 ) 
        {
            lep_eta_max = fabs(vSelectedLeptons.at(0).eta());

            if ( fabs(vSelectedLeptons.at(1).eta()) > lep_eta_max ) lep_eta_max = fabs(vSelectedLeptons.at(1).eta());

            for (int ijet=0; ijet < vJet->size() ; ijet++) 
            {
                dr = GetDeltaR( vSelectedLeptons.at(0).eta(), vSelectedLeptons.at(0).phi(), vJet->at(ijet).eta(), vJet->at(ijet).phi() );
                if ( dr < lep1_dr_min ) lep1_dr_min = dr;

                dr = GetDeltaR( vSelectedLeptons.at(1).eta(), vSelectedLeptons.at(1).phi(), vJet->at(ijet).eta(), vJet->at(ijet).phi() );
                if ( dr < lep2_dr_min ) lep2_dr_min = dr;
            }

            lep1_mtw = sqrt( 2 * vSelectedLeptons.at(0).p4().Pt() * MET * (1 - cos( vSelectedLeptons.at(0).phi() - METphi )));
        }	   

        int njj = 0; 
        float jet_dr_av = 0.;
        for (int ijet=0; ijet < vJet->size()-1 ; ijet++) 
        {
            for (int kjet=ijet+1; kjet < vJet->size() ; kjet++) 
            {
                jet_dr_av += GetDeltaR( vJet->at(ijet).eta(), vJet->at(ijet).phi(), vJet->at(kjet).eta(), vJet->at(kjet).phi() );
                njj++;
            }
        }
        if ( njj > 0 ) jet_dr_av = jet_dr_av / njj;

        // ########################################################
        // #  _     _     _                                       #
        // # | |__ (_)___| |_ ___   __ _ _ __ __ _ _ __ ___  ____ #
        // # | '_ \| / __| __/ _ \ / _` | '__/ _` | '_ ` _ \|_  / #
        // # | | | | \__ \ || (_) | (_| | | | (_| | | | | | |/ /  #
        // # |_| |_|_|___/\__\___/ \__, |_|  \__,_|_| |_| |_/___| #
        // #                       |___/                          #
        // #                                                      #
        // ########################################################

        // TTH3l
        if ( is_3l_TTH_SR ) {

            if ( _isdata ) {
                std::cout << "Run Event " << vEvent->at(0).run() <<  " " << vEvent->at(0).id() 
                    << " nlep/tight/fake " << vLeptons.size() <<  "/" << vSelectedLeptons.size() <<  "/" <<  vFakeLeptons.size() 
                    << " njet/loose/medium " << vSelectedJets.size() 
                    <<  "/" << nLooseBJets <<  "/" << nMediumBJets << std::endl;
                std::cout << "lep id " << lepid3 << " pT " << vSelectedLeptons.at(0).pt() <<  " " 
                    << vSelectedLeptons.at(1).pt() <<  " "
                    << vSelectedLeptons.at(2).pt() <<  " " << std::endl;
                std::cout << " " << std::endl;
            }

            theHistoManager->fillHisto("nLep",   "TTH3l", "", "", vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("nFake",  "TTH3l", "", "", vFakeLeptons.size()+vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("lep1Pt", "TTH3l", "", "", vSelectedLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep2Pt", "TTH3l", "", "", vSelectedLeptons.at(1).pt(), weight);
            theHistoManager->fillHisto("lep3Pt", "TTH3l", "", "", vSelectedLeptons.at(2).pt(), weight);
            theHistoManager->fillHisto("lep1Eta", "TTH3l", "", "", vSelectedLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lep2Eta", "TTH3l", "", "", vSelectedLeptons.at(1).eta(), weight);
            theHistoManager->fillHisto("lep3Eta", "TTH3l", "", "", vSelectedLeptons.at(2).eta(), weight);
            if ( vSelectedLeptons.size() >= 4 ) {
                theHistoManager->fillHisto("lep4Pt",  "TTH3l", "", "", vSelectedLeptons.at(3).pt(), weight);
                theHistoManager->fillHisto("lep4Eta", "TTH3l", "", "", vSelectedLeptons.at(3).eta(), weight);
                theHistoManager->fillHisto("lepQ",    "TTH3l", "", "", vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge()+vSelectedLeptons.at(2).charge()+vSelectedLeptons.at(3).charge(), weight);
            }
            else theHistoManager->fillHisto("lepQ", "TTH3l", "", "", vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge()+vSelectedLeptons.at(2).charge(), weight);
            theHistoManager->fillHisto("nJets",     "TTH3l", "", "", vSelectedJets.size(), weight);
            theHistoManager->fillHisto("nLooseB",   "TTH3l", "", "", nLooseBJets, weight);
            theHistoManager->fillHisto("nMediumB",  "TTH3l", "", "", nMediumBJets, weight);
            for (unsigned int ijet=0; ijet < vJet->size() ; ijet++) {
                theHistoManager->fillHisto("JetPt",  "TTH3l", "", "", vJet->at(ijet).pt(), weight);
                theHistoManager->fillHisto("JetEta", "TTH3l", "", "", vJet->at(ijet).eta(), weight);
            }
            theHistoManager->fillHisto("METpx",  "TTH3l", "", "", METx, weight);
            theHistoManager->fillHisto("METpy",  "TTH3l", "", "", METy, weight);
            theHistoManager->fillHisto("MET"  ,  "TTH3l", "", "", MET, weight);
            theHistoManager->fillHisto("METphi", "TTH3l", "", "", METphi, weight);
            theHistoManager->fillHisto("METsum", "TTH3l", "", "", METsum, weight);
            theHistoManager->fillHisto("MHT",    "TTH3l", "", "", MHT, weight);
            if ( vSelectedJets.size() <= 3 ) theHistoManager->fillHisto("MET LD", "TTH3l", "", "", met_ld, weight);
            theHistoManager->fillHisto("lep1 Mt",     "TTH3l", "", "", lep1_mtw, weight);
            theHistoManager->fillHisto("lep1 dRmin",  "TTH3l", "", "", lep1_dr_min, weight);
            theHistoManager->fillHisto("lep2 dRmin",  "TTH3l", "", "", lep2_dr_min, weight);
            theHistoManager->fillHisto("lep Eta max", "TTH3l", "", "", lep_eta_max, weight);
            theHistoManager->fillHisto("jet dR av",   "TTH3l", "", "", jet_dr_av, weight);
            theHistoManager->fillHisto("LepId", "TTH3l", "", "", lepid3, weight);
            theHistoManager->fillHisto("HT",     "TTH3l", "", "", jetht, weight);
            theHistoManager->fillHisto("Mll min","TTH3l", "", "", Mllmin, weight);
            theHistoManager->fillHisto("best MZ","TTH3l", "", "", Mllbest, weight);
        }

        // WZ
        if ( is_3l_WZ_CR ) {
            theHistoManager->fillHisto("nLep",   "WZ", "", "", vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("nFake",  "WZ", "", "", vFakeLeptons.size()+vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("lep1Pt", "WZ", "", "", vSelectedLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep2Pt", "WZ", "", "", vSelectedLeptons.at(1).pt(), weight);
            theHistoManager->fillHisto("lep3Pt", "WZ", "", "", vSelectedLeptons.at(2).pt(), weight);
            theHistoManager->fillHisto("lep1Eta", "WZ", "", "", vSelectedLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lep2Eta", "WZ", "", "", vSelectedLeptons.at(1).eta(), weight);
            theHistoManager->fillHisto("lep3Eta", "WZ", "", "", vSelectedLeptons.at(2).eta(), weight);
            theHistoManager->fillHisto("lepQ", "WZ", "", "", vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge()+vSelectedLeptons.at(2).charge(), weight);
            int lepw = -1;
            float zpt = -1.;
            for (int i=0; i<vSelectedLeptons.size()-1; i++) {
                for (int j=i+1; j<vSelectedLeptons.size(); j++) {
                    if ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id()) {
                        float mz = ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M();
                        if ( fabs(mz - 91.188) < 10. ) {
                            if ( i == 0 && j == 1 ) lepw = 2;
                            if ( i == 0 && j == 2 ) lepw = 1;    
                            if ( i == 1 && j == 2 ) lepw = 0;
                            zpt = ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).Pt();
                        }
                    }
                }
            }
            theHistoManager->fillHisto("nJets",     "WZ", "", "", vSelectedJets.size(), weight);
            theHistoManager->fillHisto("nLooseB",   "WZ", "", "", nLooseBJets, weight);
            theHistoManager->fillHisto("nMediumB",  "WZ", "", "", nMediumBJets, weight);
            for (unsigned int ijet=0; ijet < vJet->size() ; ijet++) {
                theHistoManager->fillHisto("JetPt",  "WZ", "", "", vJet->at(ijet).pt(), weight);
                theHistoManager->fillHisto("JetEta", "WZ", "", "", vJet->at(ijet).eta(), weight);
            }
            theHistoManager->fillHisto("METpx",  "WZ", "", "", METx, weight);
            theHistoManager->fillHisto("METpy",  "WZ", "", "", METy, weight);
            theHistoManager->fillHisto("MET"  ,  "WZ", "", "", MET, weight);
            theHistoManager->fillHisto("METphi", "WZ", "", "", METphi, weight);
            theHistoManager->fillHisto("METsum", "WZ", "", "", METsum, weight);
            theHistoManager->fillHisto("MHT",    "WZ", "", "", MHT, weight);
            if ( vSelectedJets.size() <= 3 ) theHistoManager->fillHisto("MET LD", "WZ", "", "", met_ld, weight);
            theHistoManager->fillHisto("lep1 Mt",     "WZ", "", "", lep1_mtw, weight);
            theHistoManager->fillHisto("lep1 dRmin",  "WZ", "", "", lep1_dr_min, weight);
            theHistoManager->fillHisto("lep2 dRmin",  "WZ", "", "", lep2_dr_min, weight);
            theHistoManager->fillHisto("lep Eta max", "WZ", "", "", lep_eta_max, weight);
            theHistoManager->fillHisto("jet dR av",   "WZ", "", "", jet_dr_av, weight);
            if ( lepw >= 0 ) {
                float mtw = sqrt( 2 * vSelectedLeptons.at(lepw).p4().Pt() * MET * (1 - cos( vSelectedLeptons.at(lepw).phi() - METphi )));
                theHistoManager->fillHisto("W Mt", "WZ", "", "", mtw, weight);
                theHistoManager->fillHisto("Z Pt", "WZ", "", "", zpt, weight);
            }
            TLorentzVector all_lep_invmass_p4;
            for (int i=0; i<vSelectedLeptons.size(); i++) {
                all_lep_invmass_p4 += vSelectedLeptons.at(i).p4();
            }
            float all_lep_invmass = all_lep_invmass_p4.M();
            float all_lep_sumofpt = sqrt( (lep_px*lep_px) + (lep_py*lep_py) );
            theHistoManager->fillHisto("inv.mass(l)", "WZ", "", "", all_lep_invmass, weight);
            theHistoManager->fillHisto("Pt Sum(l)",   "WZ", "", "", all_lep_sumofpt, weight);
            theHistoManager->fillHisto("LepId", "WZ", "", "", lepid3, weight);
            theHistoManager->fillHisto("HT",     "WZ", "", "", jetht, weight);
            theHistoManager->fillHisto("Mll min","WZ", "", "", Mllmin, weight);
            theHistoManager->fillHisto("best MZ","WZ", "", "", Mllbest, weight);
        }

        // WZrelaxed
        if ( is_3l_WZrel_CR ) {
            theHistoManager->fillHisto("nLep",   "WZrelaxed", "", "", vLeptons.size(), weight);
            theHistoManager->fillHisto("nFake",  "WZrelaxed", "", "", vFakeLeptons.size()+vLeptons.size(), weight);
            theHistoManager->fillHisto("lep1Pt", "WZrelaxed", "", "", vLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep2Pt", "WZrelaxed", "", "", vLeptons.at(1).pt(), weight);
            theHistoManager->fillHisto("lep3Pt", "WZrelaxed", "", "", vLeptons.at(2).pt(), weight);
            theHistoManager->fillHisto("lep1Eta", "WZrelaxed", "", "", vLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lep2Eta", "WZrelaxed", "", "", vLeptons.at(1).eta(), weight);
            theHistoManager->fillHisto("lep3Eta", "WZrelaxed", "", "", vLeptons.at(2).eta(), weight);
            theHistoManager->fillHisto("lepQ", "WZrelaxed", "", "", vLeptons.at(0).charge()+vLeptons.at(1).charge()+vLeptons.at(2).charge(), weight);
            int lepw = -1;
            float zpt = -1.;
            for (int i=0; i<vLeptons.size()-1; i++) {
                for (int j=i+1; j<vLeptons.size(); j++) {
                    if ( vLeptons.at(i).id() == -vLeptons.at(j).id()) {
                        float mz = ( vLeptons.at(i).p4() + vLeptons.at(j).p4() ).M();
                        if ( fabs(mz - 91.188) < 10. ) {
                            if ( i == 0 && j == 1 ) lepw = 2;
                            if ( i == 0 && j == 2 ) lepw = 1;    
                            if ( i == 1 && j == 2 ) lepw = 0;
                            zpt = ( vLeptons.at(i).p4() + vLeptons.at(j).p4() ).Pt();
                        }
                    }
                }
            }
            theHistoManager->fillHisto("nJets",     "WZrelaxed", "", "", vSelectedJets.size(), weight);
            theHistoManager->fillHisto("nLooseB",   "WZrelaxed", "", "", nLooseBJets, weight);
            theHistoManager->fillHisto("nMediumB",  "WZrelaxed", "", "", nMediumBJets, weight);
            for (unsigned int ijet=0; ijet < vJet->size() ; ijet++) {
                theHistoManager->fillHisto("JetPt",  "WZrelaxed", "", "", vJet->at(ijet).pt(), weight);
                theHistoManager->fillHisto("JetEta", "WZrelaxed", "", "", vJet->at(ijet).eta(), weight);
            }
            theHistoManager->fillHisto("METpx",  "WZrelaxed", "", "", METx, weight);
            theHistoManager->fillHisto("METpy",  "WZrelaxed", "", "", METy, weight);
            theHistoManager->fillHisto("MET"  ,  "WZrelaxed", "", "", MET, weight);
            theHistoManager->fillHisto("METphi", "WZrelaxed", "", "", METphi, weight);
            theHistoManager->fillHisto("METsum", "WZrelaxed", "", "", METsum, weight);
            theHistoManager->fillHisto("MHT",    "WZrelaxed", "", "", MHT, weight);
            if ( vSelectedJets.size() <= 3 ) theHistoManager->fillHisto("MET LD", "WZrelaxed", "", "", met_ld, weight);
            theHistoManager->fillHisto("lep1 Mt",     "WZrelaxed", "", "", lep1_mtw, weight);
            theHistoManager->fillHisto("lep1 dRmin",  "WZrelaxed", "", "", lep1_dr_min, weight);
            theHistoManager->fillHisto("lep2 dRmin",  "WZrelaxed", "", "", lep2_dr_min, weight);
            theHistoManager->fillHisto("lep Eta max", "WZrelaxed", "", "", lep_eta_max, weight);
            theHistoManager->fillHisto("jet dR av",   "WZrelaxed", "", "", jet_dr_av, weight);
            if ( lepw >= 0 ) {
                float mtw = sqrt( 2 * vLeptons.at(lepw).p4().Pt() * MET * (1 - cos( vLeptons.at(lepw).phi() - METphi )));
                theHistoManager->fillHisto("W Mt", "WZrelaxed", "", "", mtw, weight);
                theHistoManager->fillHisto("Z Pt", "WZrelaxed", "", "", zpt, weight);
            }
            TLorentzVector all_lep_invmass_p4;
            for (int i=0; i<vLeptons.size(); i++) {
                all_lep_invmass_p4 += vLeptons.at(i).p4();
            }
            float all_lep_invmass = all_lep_invmass_p4.M();
            float all_lep_sumofpt = sqrt( (lep_px*lep_px) + (lep_py*lep_py) );
            theHistoManager->fillHisto("inv.mass(l)", "WZrelaxed", "", "", all_lep_invmass, weight);
            theHistoManager->fillHisto("Pt Sum(l)",   "WZrelaxed", "", "", all_lep_sumofpt, weight);
            theHistoManager->fillHisto("LepId", "WZrelaxed", "", "", lepid3, weight);
            theHistoManager->fillHisto("HT",     "WZrelaxed", "", "", jetht, weight);
            theHistoManager->fillHisto("Mll min","WZrelaxed", "", "", Mllmin, weight);
            theHistoManager->fillHisto("best MZ","WZrelaxed", "", "", Mllbest, weight);
        }

        // TTZ
        if ( is_3l_TTZ_CR ) {
            theHistoManager->fillHisto("nLep",    "TTZ", "", "", vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("nFake",   "TTZ", "", "", vFakeLeptons.size()+vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("lep1Pt",  "TTZ", "", "", vSelectedLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep2Pt",  "TTZ", "", "", vSelectedLeptons.at(1).pt(), weight);
            theHistoManager->fillHisto("lep3Pt",  "TTZ", "", "", vSelectedLeptons.at(2).pt(), weight);
            theHistoManager->fillHisto("lep1Eta", "TTZ", "", "", vSelectedLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lep2Eta", "TTZ", "", "", vSelectedLeptons.at(1).eta(), weight);
            theHistoManager->fillHisto("lep3Eta", "TTZ", "", "", vSelectedLeptons.at(2).eta(), weight);
            if ( vSelectedLeptons.size() >= 4 ) {
                theHistoManager->fillHisto("lep4Pt",  "TTZ", "", "", vSelectedLeptons.at(3).pt(), weight);
                theHistoManager->fillHisto("lep4Eta", "TTZ", "", "", vSelectedLeptons.at(3).eta(), weight);
                theHistoManager->fillHisto("lepQ",    "TTZ", "", "", vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge()+vSelectedLeptons.at(2).charge()+vSelectedLeptons.at(3).charge(), weight);
            }
            else theHistoManager->fillHisto("lepQ", "TTZ", "", "", vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge()+vSelectedLeptons.at(2).charge(), weight);
            theHistoManager->fillHisto("nJets",     "TTZ", "", "", vSelectedJets.size(), weight);
            theHistoManager->fillHisto("nLooseB",   "TTZ", "", "", nLooseBJets, weight);
            theHistoManager->fillHisto("nMediumB",  "TTZ", "", "", nMediumBJets, weight);
            for (unsigned int ijet=0; ijet < vJet->size() ; ijet++) {
                theHistoManager->fillHisto("JetPt",  "TTZ", "", "", vJet->at(ijet).pt(), weight);
                theHistoManager->fillHisto("JetEta", "TTZ", "", "", vJet->at(ijet).eta(), weight);
            }
            theHistoManager->fillHisto("METpx",  "TTZ", "", "", METx, weight);
            theHistoManager->fillHisto("METpy",  "TTZ", "", "", METy, weight);
            theHistoManager->fillHisto("MET"  ,  "TTZ", "", "", MET, weight);
            theHistoManager->fillHisto("METphi", "TTZ", "", "", METphi, weight);
            theHistoManager->fillHisto("METsum", "TTZ", "", "", METsum, weight);
            theHistoManager->fillHisto("MHT",    "TTZ", "", "", MHT, weight);
            if ( vSelectedJets.size() <= 3 ) theHistoManager->fillHisto("MET LD", "TTZ", "", "", met_ld, weight);
            theHistoManager->fillHisto("lep1 Mt",     "TTZ", "", "", lep1_mtw, weight);
            theHistoManager->fillHisto("lep1 dRmin",  "TTZ", "", "", lep1_dr_min, weight);
            theHistoManager->fillHisto("lep2 dRmin",  "TTZ", "", "", lep2_dr_min, weight);
            theHistoManager->fillHisto("lep Eta max", "TTZ", "", "", lep_eta_max, weight);
            theHistoManager->fillHisto("jet dR av",   "TTZ", "", "", jet_dr_av, weight);
            theHistoManager->fillHisto("LepId", "TTZ", "", "", lepid3, weight);
            theHistoManager->fillHisto("HT",     "TTZ", "", "", jetht, weight);
            theHistoManager->fillHisto("Mll min","TTZ", "", "", Mllmin, weight);
            theHistoManager->fillHisto("best MZ","TTZ", "", "", Mllbest, weight);
        }

        // Zl
        if ( is_Zl_CR ) {
            theHistoManager->fillHisto("nLep",   "Zl", "", "", vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("nFake",  "Zl","", "",  vFakeLeptons.size()+vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("lep1Pt", "Zl", "", "", vSelectedLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep2Pt", "Zl", "", "", vSelectedLeptons.at(1).pt(), weight);
            theHistoManager->fillHisto("lep3Pt", "Zl", "", "", vFakeLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep1Eta", "Zl", "", "", vSelectedLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lep2Eta", "Zl", "", "", vSelectedLeptons.at(1).eta(), weight);
            theHistoManager->fillHisto("lep3Eta", "Zl", "", "", vFakeLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lepQ",    "Zl", "", "", vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge(), weight);
            theHistoManager->fillHisto("nJets",     "Zl", "", "", vSelectedJets.size(), weight);
            theHistoManager->fillHisto("nLooseB",   "Zl", "", "", nLooseBJets, weight);
            theHistoManager->fillHisto("nMediumB",  "Zl", "", "", nMediumBJets, weight);
            for (unsigned int ijet=0; ijet < vJet->size() ; ijet++) {
                theHistoManager->fillHisto("JetPt",  "Zl", "", "", vJet->at(ijet).pt(), weight);
                theHistoManager->fillHisto("JetEta", "Zl", "", "", vJet->at(ijet).eta(), weight);
            }
            theHistoManager->fillHisto("METpx",  "Zl", "", "", METx, weight);
            theHistoManager->fillHisto("METpy",  "Zl", "", "", METy, weight);
            theHistoManager->fillHisto("MET"  ,  "Zl", "", "", MET, weight);
            theHistoManager->fillHisto("METphi", "Zl", "", "", METphi, weight);
            theHistoManager->fillHisto("METsum", "Zl", "", "", METsum, weight);
            theHistoManager->fillHisto("MHT",    "Zl", "", "", MHT, weight);
            if ( vSelectedJets.size() <= 3 ) theHistoManager->fillHisto("MET LD", "Zl", "", "", met_ld, weight);
            theHistoManager->fillHisto("lep1 Mt",     "Zl", "", "", lep1_mtw, weight);
            theHistoManager->fillHisto("lep1 dRmin",  "Zl", "", "", lep1_dr_min, weight);
            theHistoManager->fillHisto("lep2 dRmin",  "Zl", "", "", lep2_dr_min, weight);
            theHistoManager->fillHisto("lep Eta max", "Zl", "", "", lep_eta_max, weight);
            theHistoManager->fillHisto("jet dR av",   "Zl", "", "", jet_dr_av, weight);
            theHistoManager->fillHisto("LepId", "Zl", "", "", lepid2, weight);
            theHistoManager->fillHisto("HT",     "Zl", "", "", jetht, weight);
            theHistoManager->fillHisto("Mll min","Zl", "", "", Mllmin, weight);
            theHistoManager->fillHisto("best MZ","Zl", "", "", Mllbest, weight);
        }

        // TTH2l
        if ( is_2lss_TTH_SR ) {
            theHistoManager->fillHisto("nLep",   "TTH2l", "", "", vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("nFake",  "TTH2l","", "",  vFakeLeptons.size()+vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("lep1Pt", "TTH2l", "", "", vSelectedLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep2Pt", "TTH2l", "", "", vSelectedLeptons.at(1).pt(), weight);
            theHistoManager->fillHisto("lep1Eta", "TTH2l", "", "", vSelectedLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lep2Eta", "TTH2l", "", "", vSelectedLeptons.at(1).eta(), weight);
            theHistoManager->fillHisto("lepQ", "TTH2l", "", "", vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge(), weight);
            theHistoManager->fillHisto("nJets",     "TTH2l", "", "", vSelectedJets.size(), weight);
            theHistoManager->fillHisto("nLooseB",   "TTH2l", "", "", nLooseBJets, weight);
            theHistoManager->fillHisto("nMediumB",  "TTH2l", "", "", nMediumBJets, weight);
            for (unsigned int ijet=0; ijet < vJet->size() ; ijet++) {
                theHistoManager->fillHisto("JetPt",  "TTH2l", "", "", vJet->at(ijet).pt(), weight);
                theHistoManager->fillHisto("JetEta", "TTH2l", "", "", vJet->at(ijet).eta(), weight);
            }
            theHistoManager->fillHisto("METpx",  "TTH2l", "", "", METx, weight);
            theHistoManager->fillHisto("METpy",  "TTH2l", "", "", METy, weight);
            theHistoManager->fillHisto("MET"  ,  "TTH2l", "", "", MET, weight);
            theHistoManager->fillHisto("METphi", "TTH2l", "", "", METphi, weight);
            theHistoManager->fillHisto("METsum", "TTH2l", "", "", METsum, weight);
            theHistoManager->fillHisto("MHT",    "TTH2l", "", "", MHT, weight);
            if ( vSelectedJets.size() <= 3 ) theHistoManager->fillHisto("MET LD", "TTH2l", "", "", met_ld, weight);
            theHistoManager->fillHisto("lep1 Mt",     "TTH2l", "", "", lep1_mtw, weight);
            theHistoManager->fillHisto("lep1 dRmin",  "TTH2l", "", "", lep1_dr_min, weight);
            theHistoManager->fillHisto("lep2 dRmin",  "TTH2l", "", "", lep2_dr_min, weight);
            theHistoManager->fillHisto("lep Eta max", "TTH2l", "", "", lep_eta_max, weight);
            theHistoManager->fillHisto("jet dR av",   "TTH2l", "", "", jet_dr_av, weight);
            theHistoManager->fillHisto("LepId", "TTH2l", "", "", lepid2, weight);
            theHistoManager->fillHisto("HT",     "TTH2l", "", "", jetht, weight);
            theHistoManager->fillHisto("Mll min","TTH2l", "", "", Mllmin, weight);
            theHistoManager->fillHisto("best MZ","TTH2l", "", "", Mllbest, weight);
        }

        // TTdilep
        if ( is_emu_TT_CR ) {
            theHistoManager->fillHisto("nLep",   "TTemu", "", "", vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("nFake",  "TTemu","", "",  vFakeLeptons.size()+vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("lep1Pt", "TTemu", "", "", vSelectedLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep2Pt", "TTemu", "", "", vSelectedLeptons.at(1).pt(), weight);
            theHistoManager->fillHisto("lep1Eta", "TTemu", "", "", vSelectedLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lep2Eta", "TTemu", "", "", vSelectedLeptons.at(1).eta(), weight);
            theHistoManager->fillHisto("lepQ", "TTemu", "", "", vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge(), weight);
            theHistoManager->fillHisto("nJets",     "TTemu", "", "", vSelectedJets.size(), weight);
            theHistoManager->fillHisto("nLooseB",   "TTemu", "", "", nLooseBJets, weight);
            theHistoManager->fillHisto("nMediumB",  "TTemu", "", "", nMediumBJets, weight);
            for (unsigned int ijet=0; ijet < vJet->size() ; ijet++) {
                theHistoManager->fillHisto("JetPt",  "TTemu", "", "", vJet->at(ijet).pt(), weight);
                theHistoManager->fillHisto("JetEta", "TTemu", "", "", vJet->at(ijet).eta(), weight);
            }
            theHistoManager->fillHisto("METpx",  "TTemu", "", "", METx, weight);
            theHistoManager->fillHisto("METpy",  "TTemu", "", "", METy, weight);
            theHistoManager->fillHisto("MET"  ,  "TTemu", "", "", MET, weight);
            theHistoManager->fillHisto("METphi", "TTemu", "", "", METphi, weight);
            theHistoManager->fillHisto("METsum", "TTemu", "", "", METsum, weight);
            theHistoManager->fillHisto("MHT",    "TTemu", "", "", MHT, weight);
            if ( vSelectedJets.size() <= 3 ) theHistoManager->fillHisto("MET LD", "TTemu", "", "", met_ld, weight);
            theHistoManager->fillHisto("lep1 Mt",     "TTemu", "", "", lep1_mtw, weight);
            theHistoManager->fillHisto("lep1 dRmin",  "TTemu", "", "", lep1_dr_min, weight);
            theHistoManager->fillHisto("lep2 dRmin",  "TTemu", "", "", lep2_dr_min, weight);
            theHistoManager->fillHisto("lep Eta max", "TTemu", "", "", lep_eta_max, weight);
            theHistoManager->fillHisto("jet dR av",   "TTemu", "", "", jet_dr_av, weight);
            theHistoManager->fillHisto("LepId", "TTemu", "", "", lepid2, weight);
            theHistoManager->fillHisto("HT",     "TTemu", "", "", jetht, weight);
            theHistoManager->fillHisto("Mll min","TTemu", "", "", Mllmin, weight);
            theHistoManager->fillHisto("best MZ","TTemu", "", "", Mllbest, weight);
        }

    }

}

void TTbarHiggsMultileptonAnalysis::TwoLeptonsSameSignSelection_TTH2l(int evt)
{

    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel",   "ttH2lss",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel",   "ttH2lss",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection #
    // ####################

    bool nLep               = ( vSelectedLeptons.size()     == 2 );
    if(!nLep)               return;

    bool leading_lep_pt     = ( vSelectedLeptons.at(0).pt() > 20 );
    if(!leading_lep_pt)     return;

    //bool following_lep_pt   = ( vSelectedLeptons.at(1).pt()               > 10 );
    bool following_lep_pt   = (  ( (abs(vSelectedLeptons.at(1).id()) == 11) && (vSelectedLeptons.at(1).pt() > 15) )
            || ( (abs(vSelectedLeptons.at(1).id()) == 13) && (vSelectedLeptons.at(1).pt() > 10) ) );
    if(!following_lep_pt)   return;

    bool passMll12Gt12      = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M()  > 12);
    if(!passMll12Gt12)      return;

    bool nJets              = ( vSelectedJets.size()        >= 4 );
    if(!nJets)              return;

    bool nLooseBtag         = ( nLooseBJets                 >= 2 );
    bool nMediumBtag        = ( nMediumBJets                >= 1 );
    if(!nLooseBtag && !nMediumBtag)      return;


    // ###############################
    // # Two leptons event selection #
    // ###############################

    bool same_sign      = ( vSelectedLeptons.at(0).charge() == vSelectedLeptons.at(1).charge() );
    if(!same_sign)      return;

    // ##########
    // # Z veto # here for leptons of same charge only !
    // ##########

    bool pass_Zveto = true;
    for(int i=0; i<vSelectedLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vSelectedLeptons.size(); j++)
        {
            if (  ( vSelectedLeptons.at(i).id() == vSelectedLeptons.at(j).id()                            )
                    && ( fabs( ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M() - 91.188) < 10 ) )
            { pass_Zveto = false ;}
        }
    }

    theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "ttH2lss",   "",   2, weight);

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "ttH2lss", _sampleName.Data(),   3, weight);

    // #################################
    // # b-tagging nominal reweighting #
    // #################################

    std::vector<double> jetPts;
    std::vector<double> jetEtas;
    std::vector<double> jetCSVs;
    std::vector<int>    jetFlavors;
    int iSys = 0;
    double wgt_csv, wgt_csv_def, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf, new_weight;

    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetPts.push_back(     vSelectedJets.at(i).pt()                );
        jetEtas.push_back(    vSelectedJets.at(i).eta()               );
        jetCSVs.push_back(    vSelectedJets.at(i).CSVv2()             );
        jetFlavors.push_back( vSelectedJets.at(i).jet_hadronFlavour() );
    }

    wgt_csv_def = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, iSys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);
    new_weight = weight * wgt_csv_def; // weight = weight * wgt_csv_def;

    // ##################################################################################################################################

    // ##################################
    // # b-tagging deriving systematics #
    // ##################################

    std::vector<double> weights_csv;
    double wgt_csv_def_sys = 0;

    for(int i=7; i<25; i++)
    {
        wgt_csv_def_sys = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, i, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf)/wgt_csv_def;
        weights_csv.push_back(wgt_csv_def_sys);
    }

    double min_weight_csv = *min_element(weights_csv.begin(),weights_csv.end());
    double max_weight_csv = *max_element(weights_csv.begin(),weights_csv.end());
    theHistoManager->fillHisto("WeightCSV_min",  "",   "ttH2l",   "", min_weight_csv, 1);
    theHistoManager->fillHisto("WeightCSV_max",  "",   "ttH2l",   "", max_weight_csv, 1);

    weight_csv_down = min_weight_csv;
    weight_csv_up   = max_weight_csv;    

    // ##################################################################################################################################

    if(  (abs(vSelectedLeptons.at(0).id()) == 11)
            && (abs(vSelectedLeptons.at(1).id()) == 11)
            && (pass_Zveto                            )
            && (met_ld                           > 0.2) // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "FullThreeLeptons",   "ttH2lee",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "ttH2lee",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "ttH2lee",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "FullThreeLeptons",   "ttH2lee",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "FullThreeLeptons",   "ttH2lee",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "FullThreeLeptons",   "ttH2lee",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "FullThreeLeptons",   "ttH2lee",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "FullThreeLeptons",   "ttH2lee",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "ttH2lee",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "ttH2lee",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "ttH2lee",   "", vSelectedLeptons.size()        , weight);
    }

    if(  ( (abs(vSelectedLeptons.at(0).id()) == 11) && (abs(vSelectedLeptons.at(1).id()) == 13) )
            || ( (abs(vSelectedLeptons.at(0).id()) == 13) && (abs(vSelectedLeptons.at(1).id()) == 11) ) // smarter way to do is probably exists
            && true                                      // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "FullThreeLeptons",   "ttH2lem",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "ttH2lem",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "ttH2lem",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "FullThreeLeptons",   "ttH2lem",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "FullThreeLeptons",   "ttH2lem",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "FullThreeLeptons",   "ttH2lem",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "FullThreeLeptons",   "ttH2lem",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "FullThreeLeptons",   "ttH2lem",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "ttH2lem",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "ttH2lem",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "ttH2lem",   "", vSelectedLeptons.size()        , weight);
    }

    if(  (abs(vSelectedLeptons.at(0).id()) == 13)
            && (abs(vSelectedLeptons.at(1).id()) == 13)
            && true                                      // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "FullThreeLeptons",   "ttH2lmm",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "ttH2lmm",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "ttH2lmm",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "FullThreeLeptons",   "ttH2lmm",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "FullThreeLeptons",   "ttH2lmm",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "FullThreeLeptons",   "ttH2lmm",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "FullThreeLeptons",   "ttH2lmm",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "FullThreeLeptons",   "ttH2lmm",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "ttH2lmm",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "ttH2lmm",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "ttH2lmm",   "", vSelectedLeptons.size()        , weight);
    }

    if (   (abs(vSelectedLeptons.at(0).id()) == 11)
            && (abs(vSelectedLeptons.at(1).id()) == 11)
            && ( !pass_Zveto || met_ld < 0.2)         ) return;

    is_2lss_TTH_SR = true;   

    // ####################################
    // #  ____  ____    ____  ____ _____  #
    // # |___ \|  _ \  | __ )|  _ \_   _| #
    // #   __) | | | | |  _ \| | | || |   #
    // #  / __/| |_| | | |_) | |_| || |   #
    // # |_____|____/  |____/|____/ |_|   #
    // #                                  #
    // ####################################

    // based on https://github.com/CERN-PH-CMG/cmgtools-lite/blob/0b47d4d1c50ea0e24ef0d9cf1c24c763e78c1bf0/TTHAnalysis/python/tools/kinMVA_2D_2lss_3l.py
    // and https://github.com/CERN-PH-CMG/cmgtools-lite/blob/0b47d4d1c50ea0e24ef0d9cf1c24c763e78c1bf0/TTHAnalysis/python/tools/eventVars_2lss.py

    // ======================================================================================================
    // variables against ttbar
    max_Lep_eta     = std::max( fabs(vSelectedLeptons.at(0).eta()), fabs(vSelectedLeptons.at(1).eta()) ) ; // ok

    //numJets_float   = vSelectedJets.size() ;                                                             // ok
    nJet25_Recl = vSelectedJets.size() ;     

    mindr_lep1_jet = 1000.;

    for (int i=0; i<vSelectedJets.size(); i++)
    {
        if( DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ) < mindr_lep1_jet )
        { mindr_lep1_jet = DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ); }            // ok
    }

    mindr_lep2_jet  = 1000. ;

    for (int i=0; i<vSelectedJets.size(); i++)
    {
        if( DeltaRLeptonJet( vSelectedLeptons.at(1), vSelectedJets.at(i) ) < mindr_lep2_jet )
        { mindr_lep2_jet = DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ); }            // ok
    }

    float met_max = 400;
    met             = std::min( vEvent->at(0).metpt(), met_max ) ;                                      // ok


    int njj = 0;
    float avg_dr_jet = 0.;
    for (int ijet=0; ijet < vSelectedJets.size()-1 ; ijet++) 
    {
        for (int kjet=ijet+1; kjet < vSelectedJets.size() ; kjet++) 
        {
            avg_dr_jet += GetDeltaR( vSelectedJets.at(ijet).eta(), vSelectedJets.at(ijet).phi(), vSelectedJets.at(kjet).eta(), vSelectedJets.at(kjet).phi() );
            njj++;
        }
    }
    if ( njj > 0 ) avg_dr_jet = avg_dr_jet / njj;                                                       // ok

    MT_met_lep1     = sqrt( 2 * vSelectedLeptons.at(0).p4().Pt() * vEvent->at(0).metpt() 
            * (1 - cos( vSelectedLeptons.at(0).phi() - vEvent->at(0).metphi() )));                      // ok

    signal_2lss_TT_MVA  = mva_2lss_tt->EvaluateMVA("BDTG method");

    theHistoManager->fillHisto("Signal_2lss_TT_MVA",                       "FinalCut", "ttH2lss",   "",  signal_2lss_TT_MVA,   weight);

    // ======================================================================================================
    // variables against ttV

    LepGood_conePt0 = vSelectedLeptons.at(0).pt() ;                                                    // not clear

    LepGood_conePt1 = vSelectedLeptons.at(1).pt() ;                                                    // not clear

    //mhtJet25_Recl   = sqrt( (jet_px*jet_px) + (jet_py*jet_py) );

    signal_2lss_TTV_MVA = mva_2lss_ttV->EvaluateMVA("BDTG method");

    theHistoManager->fillHisto("Signal_2lss_TTV_MVA",                      "FinalCut", "ttH2lss",   "",  signal_2lss_TTV_MVA,  weight);

    //std::cout << " signal 2lss TT MVA: "  << signal_2lss_TT_MVA
    //    << " signal 2lss TTV MVA: " << signal_2lss_TTV_MVA << std::endl;

    // ======================================================================================================

    //fillOutputTree();

    if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);
}

void TTbarHiggsMultileptonAnalysis::TwoLeptonsSameSignSelection_LepMVA_sideband(int evt)
{

    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel",   "LepMVA_2l_SB",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel",   "LepMVA_2l_SB",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection #
    // ####################

    bool nLep               = ( vSelectedLeptons.size()                   == 2 );
    if(!nLep)               return;

    bool leading_lep_pt     = ( vSelectedLeptons.at(0).pt()               > 20 );
    if(!leading_lep_pt)     return;

    //bool following_lep_pt   = ( vSelectedLeptons.at(1).pt()               > 10 );
    bool following_lep_pt   = (  ( (abs(vSelectedLeptons.at(1).id()) == 11) && (vSelectedLeptons.at(1).pt() > 15) )
            || ( (abs(vSelectedLeptons.at(1).id()) == 13) && (vSelectedLeptons.at(1).pt() > 10) ) );
    if(!following_lep_pt)   return;

    bool passMll12Gt12      = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M()  > 12);
    if(!passMll12Gt12)      return;

    bool nJets              = ( vSelectedJets.size()                      >= 4 );
    if(!nJets)              return;

    bool nLooseBtag         = ( nLooseBJets                               >= 2 );
    bool nMediumBtag        = ( nMediumBJets                              >= 1 );
    if(!nLooseBtag && !nMediumBtag)      return;


    // ###############################
    // # Two leptons event selection #
    // ###############################

    // ##########
    // # Z veto #
    // ##########

    bool pass_Zveto = true;
    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        for(int j=0; j<vSelectedLeptons.size(); j++)
        {
            if (  ( i                           != j                                                       )
                    && ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id()                            )
                    && ( fabs( ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M() - 91.188) < 10 ) )
            { pass_Zveto = false ;}
        }
    }

    theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "LepMVA_2l_SB",   "",   2, weight);

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "LepMVA_2l_SB", _sampleName.Data(),   3, weight);

    if(  (abs(vSelectedLeptons.at(0).id()) == 11)
            && (abs(vSelectedLeptons.at(1).id()) == 11)
            && (pass_Zveto                            )
            && (met_ld                           > 0.2) // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "FullThreeLeptons",   "LepMVA_2l_SB",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "FullThreeLeptons",   "LepMVA_2l_SB",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "FullThreeLeptons",   "LepMVA_2l_SB",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "FullThreeLeptons",   "LepMVA_2l_SB",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedLeptons.size()        , weight);
    }

    if(  ( (abs(vSelectedLeptons.at(0).id()) == 11) && (abs(vSelectedLeptons.at(1).id()) == 13) )
            || ( (abs(vSelectedLeptons.at(0).id()) == 13) && (abs(vSelectedLeptons.at(1).id()) == 11) ) // smarter way to do is probably exists
            && true                                      // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "FullThreeLeptons",   "LepMVA_2l_SB",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "FullThreeLeptons",   "LepMVA_2l_SB",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "FullThreeLeptons",   "LepMVA_2l_SB",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "FullThreeLeptons",   "LepMVA_2l_SB",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedLeptons.size()        , weight);
    }

    if(  (abs(vSelectedLeptons.at(0).id()) == 13)
            && (abs(vSelectedLeptons.at(1).id()) == 13)
            && true                                      // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "FullThreeLeptons",   "LepMVA_2l_SB",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "FullThreeLeptons",   "LepMVA_2l_SB",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "FullThreeLeptons",   "LepMVA_2l_SB",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "FullThreeLeptons",   "LepMVA_2l_SB",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedLeptons.size()        , weight);
    }

    is_2lss_LepMVA_SB = true;   

    //fillOutputTree();

    if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);
}

void TTbarHiggsMultileptonAnalysis::TwoLeptonsSameSignSelection_JetMultiplicity_sideband(int evt)
{

    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel",   "JM_2l_SB",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel",   "JM_2l_SB",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection #
    // ####################

    bool nLep               = ( vSelectedLeptons.size()                   == 2 );
    if(!nLep)               return;

    bool leading_lep_pt     = ( vSelectedLeptons.at(0).pt()               > 20 );
    if(!leading_lep_pt)     return;

    bool following_lep_pt   = (  ( (abs(vSelectedLeptons.at(1).id()) == 11) && (vSelectedLeptons.at(1).pt() > 15) )
            || ( (abs(vSelectedLeptons.at(1).id()) == 13) && (vSelectedLeptons.at(1).pt() > 10) ) );
    if(!following_lep_pt)   return;

    bool passMll12Gt12      = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M()  > 12);
    if(!passMll12Gt12)      return;

    bool nJets              = ( vSelectedJets.size()                      == 3 );
    if(!nJets)              return;

    bool nLooseBtag         = ( nLooseBJets                               >= 2 );
    bool nMediumBtag        = ( nMediumBJets                              >= 1 );
    if(!nLooseBtag && !nMediumBtag)      return;


    // ###############################
    // # Two leptons event selection #
    // ###############################

    // ##########
    // # Z veto #
    // ##########

    bool pass_Zveto = true;
    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        for(int j=0; j<vSelectedLeptons.size(); j++)
        {
            if (  ( i                           != j                                                       )
                    && ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id()                            )
                    && ( fabs( ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M() - 91.188) < 10 ) )
            { pass_Zveto = false ;}
        }
    }

    theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "JM_2l_SB",   "",   2, weight);

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "JM_2l_SB", _sampleName.Data(),   3, weight);

    if(  (abs(vSelectedLeptons.at(0).id()) == 11)
            && (abs(vSelectedLeptons.at(1).id()) == 11)
            && (pass_Zveto                            )
            && (met_ld                           > 0.2) // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "FullThreeLeptons",   "JM_2l_SB",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "FullThreeLeptons",   "JM_2l_SB",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "FullThreeLeptons",   "JM_2l_SB",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "FullThreeLeptons",   "JM_2l_SB",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedLeptons.size()        , weight);
    }

    if(  ( (abs(vSelectedLeptons.at(0).id()) == 11) && (abs(vSelectedLeptons.at(1).id()) == 13) )
            || ( (abs(vSelectedLeptons.at(0).id()) == 13) && (abs(vSelectedLeptons.at(1).id()) == 11) ) // smarter way to do is probably exists
            && true                                      // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "FullThreeLeptons",   "JM_2l_SB",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "FullThreeLeptons",   "JM_2l_SB",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "FullThreeLeptons",   "JM_2l_SB",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "FullThreeLeptons",   "JM_2l_SB",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedLeptons.size()        , weight);
    }

    if(  (abs(vSelectedLeptons.at(0).id()) == 13)
            && (abs(vSelectedLeptons.at(1).id()) == 13)
            && true                                      // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "FullThreeLeptons",   "JM_2l_SB",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "FullThreeLeptons",   "JM_2l_SB",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "FullThreeLeptons",   "JM_2l_SB",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "FullThreeLeptons",   "JM_2l_SB",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedLeptons.size()        , weight);
    }

    is_2lss_JM_SB = true;   

    //fillOutputTree();

    if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);
}

void TTbarHiggsMultileptonAnalysis::DiLeptonSelection_TT_CR(int evt)
{

    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel",   "TT_2l_CR",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel",   "TT_2l_CR",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection #
    // ####################

    bool nLep               = ( vSelectedLeptons.size()                   == 2 );
    if(!nLep)               return;

    bool leading_lep_pt     = ( vSelectedLeptons.at(0).pt()               > 20 );
    if(!leading_lep_pt)     return;

    bool following_lep_pt   = (  ( (abs(vSelectedLeptons.at(1).id()) == 11) && (vSelectedLeptons.at(1).pt() > 15) )
            || ( (abs(vSelectedLeptons.at(1).id()) == 13) && (vSelectedLeptons.at(1).pt() > 10) ) );
    if(!following_lep_pt)   return;

    bool passMll12Gt12      = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M()  > 12);
    if(!passMll12Gt12)      return;

    bool nJets              = ( vSelectedJets.size()                      >= 2 );
    if(!nJets)              return;

    bool nMediumBtag        = ( nMediumBJets                              >= 1 );
    if(!nMediumBtag)      return;

    // ###############################
    // #        e+mu- selection      #
    // ###############################

    if ( vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge() != 0 ) return;
    if ( vSelectedLeptons.at(0).id() == -vSelectedLeptons.at(1).id() ) return;

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if(  ( (abs(vSelectedLeptons.at(0).id()) == 11) && (abs(vSelectedLeptons.at(1).id()) == 13) )
            || ( (abs(vSelectedLeptons.at(0).id()) == 13) && (abs(vSelectedLeptons.at(1).id()) == 11) ) // smarter way to do this probably exists
            && (  vSelectedLeptons.at(0).charge()         == -vSelectedLeptons.at(1).charge()         )
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "FullThreeLeptons",   "TT_2l_CR",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "TT_2l_CR",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "TT_2l_CR",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "FullThreeLeptons",   "TT_2l_CR",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MetLD",                            "FullThreeLeptons",   "TT_2l_SB",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "FullThreeLeptons",   "TT_2l_CR",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "FullThreeLeptons",   "TT_2l_CR",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "TT_2l_CR",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "TT_2l_CR",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "TT_2l_CR",   "", vSelectedLeptons.size()        , weight);
    }

    is_emu_TT_CR = true;   

    //fillOutputTree();

    if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);
}



void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_TTH3l(int evt)
{

    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel",   "ttH3l",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel",   "ttH3l",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection #
    // ####################

    bool nLep             = ( vSelectedLeptons.size()                   >= 2 );
    if(!nLep)             return;

    bool leading_lep_pt   = ( vSelectedLeptons.at(0).pt()               > 20 );
    if(!leading_lep_pt)   return;

    bool following_lep_pt = ( vSelectedLeptons.at(1).pt()               > 10 );
    if(!following_lep_pt) return;

    bool passMll12Gt12    = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M()  > 12);
    if(!passMll12Gt12)    return;

    bool nJets            = ( vSelectedJets.size()                      >= 2 );
    if(!nJets)            return;

    bool nLooseBtag       = ( nLooseBJets                               >= 2 );
    bool nMediumBtag      = ( nMediumBJets                              >= 1 );
    if(!nLooseBtag && !nMediumBtag)      return;


    // #################################
    // # Three leptons event selection #
    // #################################

    nLep                  = ( vSelectedLeptons.size()				    >= 3 );
    //nLep                  = ( vSelectedLeptons.size()                 == 3 );
    if(!nLep)             return;

    bool third_lep_pt     = ( vSelectedLeptons.at(2).pt()               > 10 );
    if(!third_lep_pt)     return;

    theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "ttH3l",   "",    1, weight);

    passMll12Gt12  = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(2).p4() ).M()  > 12
            && ( vSelectedLeptons.at(1).p4() + vSelectedLeptons.at(2).p4() ).M()  > 12 );
    if(!passMll12Gt12)    return;

    // ##########
    // # Z veto # with tight or loose leptons ???
    // ##########

    bool pass_Zveto = true;
    for(int i=0; i<vSelectedLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vSelectedLeptons.size(); j++)
        {
            if (  ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id()                            )
                    && ( fabs( ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M() - 91.188) < 10 ) )
            { pass_Zveto = false ;}
        }
    }
    if(!pass_Zveto)       return;

    theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "ttH3l",   "",   2, weight);

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "ttH3l", _sampleName.Data(),   3, weight);
    if( vSelectedJets.size() >= 4  )  theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "ttH3l", _sampleName.Data(),   4, weight);

    // ########
    // # SFOS #
    // ########

    bool isSFOS = false;
    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        for(int j=0; j<vSelectedLeptons.size(); j++)
        {
            if (  ( i                           != j                            )
                    && ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id() ) )
            { isSFOS = true ;}
        }
    }

    if( met_ld               > 0.3 ) theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "ttH3l",   "",   6, weight);
    if( isSFOS                     ) theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "ttH3l",   "",   7, weight);

    if(vSelectedJets.size() < 4 && (met_ld < (0.2 + 0.1 * isSFOS)) ) return;

    int sum_charges = 0;
    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        sum_charges = sum_charges + vSelectedLeptons.at(i).charge();
    }

    int sum_charges_3l = 0;
    for(int i=0; i<3; i++)
    {
        sum_charges_3l = sum_charges_3l + vSelectedLeptons.at(i).charge();
    }
    if( fabs(sum_charges_3l) != 1 ) return;

    // #################################
    // # b-tagging nominal reweighting #
    // #################################

    std::vector<double> jetPts;
    std::vector<double> jetEtas;
    std::vector<double> jetCSVs;
    std::vector<int>    jetFlavors;
    int iSys = 0;
    double wgt_csv, wgt_csv_def, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf, new_weight;

    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetPts.push_back(     vSelectedJets.at(i).pt()                );
        jetEtas.push_back(    vSelectedJets.at(i).eta()               );
        jetCSVs.push_back(    vSelectedJets.at(i).CSVv2()             );
        jetFlavors.push_back( vSelectedJets.at(i).jet_hadronFlavour() );
    }

    wgt_csv_def = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, iSys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);
    new_weight = weight * wgt_csv_def; // weight = weight * wgt_csv_def;

    // ##################################################################################################################################

    // ##################################
    // # b-tagging deriving systematics #
    // ##################################

    std::vector<double> weights_csv;
    double wgt_csv_def_sys = 0;

    for(int i=7; i<25; i++)
    {
        wgt_csv_def_sys = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, i, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf)/wgt_csv_def;
        weights_csv.push_back(wgt_csv_def_sys);
    }

    double min_weight_csv = *min_element(weights_csv.begin(),weights_csv.end()); 
    double max_weight_csv = *max_element(weights_csv.begin(),weights_csv.end()); 
    theHistoManager->fillHisto("WeightCSV_min",  "",   "ttH2l",   "", min_weight_csv, 1);
    theHistoManager->fillHisto("WeightCSV_max",  "",   "ttH2l",   "", max_weight_csv, 1);

    weight_csv_down = min_weight_csv;
    weight_csv_up   = max_weight_csv;

    // ##################################################################################################################################

    theHistoManager->fillHisto("CutFlow",                          "FullThreeLeptons",   "ttH3l",   "", 8                              , weight);

    theHistoManager->fillHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "ttH3l",   "", vSelectedLeptons.at(0).pt()    , weight);
    theHistoManager->fillHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "ttH3l",   "", vSelectedLeptons.at(1).pt()    , weight);
    theHistoManager->fillHisto("ThirdLeptonPt",                    "FullThreeLeptons",   "ttH3l",   "", vSelectedLeptons.at(2).pt()    , weight);
    theHistoManager->fillHisto("MET",                              "FullThreeLeptons",   "ttH3l",   "", vEvent->at(0).metpt()          , weight);
    theHistoManager->fillHisto("MHT",                              "FullThreeLeptons",   "ttH3l",   "", MHT                            , weight);
    theHistoManager->fillHisto("MetLD",                            "FullThreeLeptons",   "ttH3l",   "", met_ld                         , weight);
    theHistoManager->fillHisto("TauMultiplicity",                  "FullThreeLeptons",   "ttH3l",   "", vSelectedTaus.size()           , weight);
    theHistoManager->fillHisto("JetMultiplicity",                  "FullThreeLeptons",   "ttH3l",   "", vSelectedJets.size()           , weight);
    theHistoManager->fillHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "ttH3l",   "", vSelectedBTagJets.size()       , weight);
    theHistoManager->fillHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "ttH3l",   "", vSelectedMediumBTagJets.size() , weight);
    theHistoManager->fillHisto("SumOfLeptonsCharges",              "FullThreeLeptons",   "ttH3l",   "", sum_charges                    , weight);
    theHistoManager->fillHisto("SumOfThreeLeptonsCharges",         "FullThreeLeptons",   "ttH3l",   "", sum_charges_3l                 , weight);
    theHistoManager->fillHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "ttH3l",   "", vSelectedLeptons.size()        , weight);

    is_3l_TTH_SR = true;   


    // ####################################
    // #  ____  ____    ____  ____ _____  #
    // # |___ \|  _ \  | __ )|  _ \_   _| #
    // #   __) | | | | |  _ \| | | || |   #
    // #  / __/| |_| | | |_) | |_| || |   #
    // # |_____|____/  |____/|____/ |_|   #
    // #                                  #
    // ####################################

    // based on https://github.com/CERN-PH-CMG/cmgtools-lite/blob/0b47d4d1c50ea0e24ef0d9cf1c24c763e78c1bf0/TTHAnalysis/python/tools/kinMVA_2D_2lss_3l.py

    // ======================================================================================================
    // variables against ttbar

    max_Lep_eta     = std::max( fabs(vSelectedLeptons.at(0).eta()), fabs(vSelectedLeptons.at(1).eta()) ) ;     // ok

    MT_met_lep1     = sqrt( 2 * vSelectedLeptons.at(0).p4().Pt() * vEvent->at(0).metpt() 
                      * (1 - cos( vSelectedLeptons.at(0).phi() - vEvent->at(0).metphi() )));                 // ok

    //numJets_float   = vSelectedJets.size() ;                                                                 // ok
    nJet25_Recl = vSelectedJets.size() ;

    mhtJet25_Recl   = sqrt( (jet_px*jet_px) + (jet_py*jet_py) );                                             // ok

    int njj = 0;
    float avg_dr_jet = 0.;
    for (int ijet=0; ijet < vSelectedJets.size()-1 ; ijet++)
    {
        for (int kjet=ijet+1; kjet < vSelectedJets.size() ; kjet++)
        {
            avg_dr_jet += GetDeltaR( vSelectedJets.at(ijet).eta(), vSelectedJets.at(ijet).phi(), vSelectedJets.at(kjet).eta(), vSelectedJets.at(kjet).phi() );
            njj++;
        }
    }
    if ( njj > 0 ) avg_dr_jet = avg_dr_jet / njj;                                                           // ok

    mindr_lep1_jet = 1000.;

    for (int i=0; i<vSelectedJets.size(); i++)
    {
        if( DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ) < mindr_lep1_jet )
        { mindr_lep1_jet = DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ); }                // ok
    }

    mindr_lep2_jet  = 1000. ;

    for (int i=0; i<vSelectedJets.size(); i++)
    {
        if( DeltaRLeptonJet( vSelectedLeptons.at(2), vSelectedJets.at(i) ) < mindr_lep2_jet )
        { mindr_lep2_jet = DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ); }               // ok
    }

    signal_3l_TT_MVA    = mva_3l_tt->EvaluateMVA("BDTG method");

    theHistoManager->fillHisto("Signal_3l_TT_MVA",                         "FinalCut",   "ttH3l",   "",  signal_3l_TT_MVA,   weight);

    // ======================================================================================================
    // variables against ttV

    met             = vEvent->at(0).metpt() ;

    LepGood_conePt0 = vSelectedLeptons.at(0).pt() ;

    LepGood_conePt1 = vSelectedLeptons.at(2).pt() ;

    signal_3l_TTV_MVA   = mva_3l_ttV->EvaluateMVA("BDTG method");

    theHistoManager->fillHisto("Signal_3l_TTV_MVA",                        "FinalCut",   "ttH3l",   "",  signal_3l_TTV_MVA,  weight);

    //std::cout << " signal 3l   TT MVA: "  << signal_3l_TT_MVA
    //          << " signal 3l   TTV MVA: " << signal_3l_TTV_MVA << std::endl;

    // ======================================================================================================

    //fillOutputTree();

    if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);
}

void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_CR_WZ(int evt)
{

    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel", "WZZZ_CR",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel", "WZZZ_CR",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection #
    // ####################

    bool nLep             = ( vSelectedLeptons.size()                   >= 2 );
    if(!nLep)             return;

    bool leading_lep_pt   = ( vSelectedLeptons.at(0).pt()               > 20 );
    if(!leading_lep_pt)   return;

    bool following_lep_pt = ( vSelectedLeptons.at(1).pt()               > 10 );
    if(!following_lep_pt) return;

    bool passMll12Gt12    = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M()  > 12);
    if(!passMll12Gt12)    return;

    bool nJets            = ( vSelectedJets.size()                      >= 2 );
    if(!nJets)            return;

    // #################################
    // # Three leptons event selection #
    // #################################

    nLep                  = ( vSelectedLeptons.size()		        >= 3 );
    if(!nLep)             return;

    bool third_lep_pt     = ( vSelectedLeptons.at(2).pt()               > 10 );
    if(!third_lep_pt)     return;

    passMll12Gt12  = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(2).p4() ).M()  > 12
            && ( vSelectedLeptons.at(1).p4() + vSelectedLeptons.at(2).p4() ).M()  > 12 );
    if(!passMll12Gt12)    return;

    theHistoManager->fillHisto("CutFlow",                           "noSel", "WZZZ_CR",   "",   2, weight);
    theHistoManager->fillHisto("CutFlow",                 "PassingTightMVA", "WZZZ_CR",   "",   1, weight);

    // ##########
    // # Z peak #
    // ##########

    bool pass_Zpeak = false;
    float Delta = 10, ZM = -1., Zpt = -1., Lep1Z = -1, Lep2Z = -1, LepW = -1;

    for(int i=0; i<vSelectedLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vSelectedLeptons.size(); j++)
        {
            if (  ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id()                               )
                    && ( fabs( ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M() - 91.188) < Delta ) )
            {
                Lep1Z      = i;
                Lep2Z      = j;
                ZM         = ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M();
                Zpt        = ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).Pt();
                Delta      = fabs(ZM - 91.188);
                pass_Zpeak = true;
            }
        }
    }
    if(!pass_Zpeak)       return;

    if( ( (Lep1Z == 0) && (Lep2Z == 1) ) || ( (Lep1Z == 1) && (Lep2Z == 0) ) ) LepW = 2;
    if( ( (Lep1Z == 0) && (Lep2Z == 2) ) || ( (Lep1Z == 2) && (Lep2Z == 0) ) ) LepW = 1;    
    if( ( (Lep1Z == 1) && (Lep2Z == 2) ) || ( (Lep1Z == 2) && (Lep2Z == 1) ) ) LepW = 0;

    float MTW = 0. ;
    if (LepW >=0) MTW = sqrt( 2 * vSelectedLeptons.at(LepW).p4().Pt() * vEvent->at(0).metpt() * (1 - cos( vSelectedLeptons.at(LepW).phi() - vEvent->at(0).metphi() )));

    theHistoManager->fillHisto("CutFlow",                           "noSel", "WZZZ_CR",   "",   3, weight);
    theHistoManager->fillHisto("CutFlow",                        "PassingZ", "WZZZ_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass",        "PassingZ", "WZZZ_CR",   "",  ZM, weight);

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                    "PassingMETLD", "WZZZ_CR",   "",   1, weight);
    if( vSelectedJets.size() >= 4  )  theHistoManager->fillHisto("CutFlow",                     "PassingJets", "WZZZ_CR",   "",   1, weight); 

    if( met_ld < 0.2 && vSelectedJets.size() < 4 ) return;
    theHistoManager->fillHisto("CutFlow",                 "PassingMETLDorJets", "WZZZ_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass", "PassingMETLDorJets", "WZZZ_CR",   "",  ZM, weight);

    // ########
    // # SFOS #
    // ########

    bool isSFOS = false;
    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        for(int j=0; j<vSelectedLeptons.size(); j++)
        {
            if (  ( i                           != j                            )
                    && ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id() ) )
            { isSFOS = true ;}
        }
    }

    theHistoManager->fillHisto("CutFlow",                 "PassingMETLDorSFOS", "WZZZ_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass", "PassingMETLDorSFOS", "WZZZ_CR",   "",  ZM, weight);

    if(vSelectedJets.size() < 4 && (met_ld < (0.2 + 0.1 * isSFOS)) ) return;

    // ##############
    // # b-jet veto #
    // ##############

    bool nLooseBtag       = ( nLooseBJets                               == 0 );
    bool nMediumBtag      = ( nMediumBJets                              == 0 );
    if(!nLooseBtag || !nMediumBtag)      return;

    theHistoManager->fillHisto("CutFlow",                 "PassingbJetsVeto", "WZZZ_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass", "PassingbJetsVeto", "WZZZ_CR",   "",  ZM, weight);

    int sum_charges_3l = 0;
    for(int i=0; i<3; i++)
    {
        sum_charges_3l = sum_charges_3l + vSelectedLeptons.at(i).charge();
    }
    if( fabs(sum_charges_3l) != 1 ) return;

    // Building arbitrary variables with all leptons...

    TLorentzVector all_lep_invmass_p4;

    int   all_lep_sumofcharges = 0;
    float all_lep_invmass      = 0;
    float all_lep_sumofpt      = 0;

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        all_lep_sumofcharges = all_lep_sumofcharges + vSelectedLeptons.at(i).charge();
        all_lep_invmass_p4   = all_lep_invmass_p4   + vSelectedLeptons.at(i).p4();
    }

    all_lep_invmass = all_lep_invmass_p4.M();
    all_lep_sumofpt = sqrt( (lepton_px*lepton_px) + (lepton_py*lepton_py) );

    // Building arbitrary variables with remaining (after Z peak determination) leptons...

    TLorentzVector rem_lep_invmass_p4;

    int   rem_lep_sumofcharges = 0;
    float rem_lep_invmass      = 0;
    float rem_lep_sumofpt      = 0;
    float rem_lep_px           = 0;
    float rem_lep_py           = 0;


    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        if( i!=Lep1Z && i!=Lep2Z)
        {
            rem_lep_sumofcharges = rem_lep_sumofcharges + vSelectedLeptons.at(i).charge();
            rem_lep_px           = rem_lep_px           + vSelectedLeptons.at(i).p4().Px();
            rem_lep_py           = rem_lep_py           + vSelectedLeptons.at(i).p4().Py();
            rem_lep_invmass_p4   = rem_lep_invmass_p4   + vSelectedLeptons.at(i).p4();
        }
    }

    rem_lep_invmass = rem_lep_invmass_p4.M();
    rem_lep_sumofpt = sqrt( (rem_lep_px*rem_lep_px) + (rem_lep_py*rem_lep_py) );

    theHistoManager->fillHisto("ZCandidateInvariantMass",        "FinalCut", "WZZZ_CR",   "",   ZM                    , weight);
    theHistoManager->fillHisto("ZCandidateTransverseMomentum",   "FinalCut", "WZZZ_CR",   "",   Zpt                   , weight);
    theHistoManager->fillHisto("MET",                            "FinalCut", "WZZZ_CR",   "",   vEvent->at(0).metpt() , weight);
    theHistoManager->fillHisto("MHT",                            "FinalCut", "WZZZ_CR",   "",   MHT                   , weight);
    theHistoManager->fillHisto("MetLD",                          "FinalCut", "WZZZ_CR",   "",   met_ld                , weight);
    theHistoManager->fillHisto("TauMultiplicity",                "FinalCut", "WZZZ_CR",   "",   vSelectedTaus.size()  , weight);
    theHistoManager->fillHisto("JetMultiplicity",                "FinalCut", "WZZZ_CR",   "",   vSelectedJets.size()  , weight);

    theHistoManager->fillHisto("MTW",                            "FinalCut", "WZZZ_CR",   "",   MTW                   , weight);
    theHistoManager->fillHisto("InvariantMassOfSelectedLeptons", "FinalCut", "WZZZ_CR",   "",   all_lep_invmass       , weight);
    theHistoManager->fillHisto("SumOfSelectedLeptonsCharges",    "FinalCut", "WZZZ_CR",   "",   all_lep_sumofcharges  , weight);
    theHistoManager->fillHisto("SumVecPtSelectedLeptons",        "FinalCut", "WZZZ_CR",   "",   all_lep_sumofpt       , weight);

    theHistoManager->fillHisto("InvMassRemainingLepton",         "FinalCut", "WZZZ_CR",   "",   rem_lep_invmass       , weight);
    theHistoManager->fillHisto("SumOfRemainingLeptonsCharges",   "FinalCut", "WZZZ_CR",   "",   rem_lep_sumofcharges  , weight);
    theHistoManager->fillHisto("SumVecPtRemainingLeptons",       "FinalCut", "WZZZ_CR",   "",   rem_lep_sumofpt       , weight);

    theHistoManager->fillHisto2D("InvMassLastLeptonVSZMass",     "FinalCut", "WZZZ_CR",   "",   rem_lep_invmass,    ZM, weight);
    theHistoManager->fillHisto2D("SumPtLepVSZMass",              "FinalCut", "WZZZ_CR",   "",   rem_lep_sumofpt,    ZM, weight);
    theHistoManager->fillHisto2D("METLDVSZMass",                 "FinalCut", "WZZZ_CR",   "",   met_ld,             ZM, weight);

    is_3l_WZ_CR = true;
}

void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_CR_WZrelaxed(int evt)
{

    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel", "WZrel_CR",   "",    vLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel", "WZrel_CR",   "",    vLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection #
    // ####################

    bool nLep             = ( vLeptons.size()                   >= 2 );
    if(!nLep)             return;

    bool leading_lep_pt   = ( vLeptons.at(0).pt()               > 20 );
    if(!leading_lep_pt)   return;

    bool following_lep_pt = ( vLeptons.at(1).pt()               > 10 );
    if(!following_lep_pt) return;

    bool passMll12Gt12    = ( ( vLeptons.at(0).p4() + vLeptons.at(1).p4() ).M()  > 12);
    if(!passMll12Gt12)    return;

    bool nJets            = ( vSelectedJets.size()                      >= 2 );
    if(!nJets)            return;

    // #################################
    // # Three leptons event selection #
    // #################################

    nLep                  = ( vLeptons.size()				    >= 3 );
    if(!nLep)             return;

    bool third_lep_pt     = ( vLeptons.at(2).pt()               > 10 );
    if(!third_lep_pt)     return;

    passMll12Gt12  = ( ( vLeptons.at(0).p4() + vLeptons.at(2).p4() ).M()  > 12
            && ( vLeptons.at(1).p4() + vLeptons.at(2).p4() ).M()  > 12 );
    if(!passMll12Gt12)    return;

    theHistoManager->fillHisto("CutFlow",                           "noSel", "WZrel_CR",   "",   2, weight);
    theHistoManager->fillHisto("CutFlow",                 "PassingTightMVA", "WZrel_CR",   "",   1, weight);

    // ##########
    // # Z peak #
    // ##########

    bool pass_Zpeak = false;
    float Delta = 10, ZM = -1., Zpt = -1., Lep1Z = -1, Lep2Z = -1;

    for(int i=0; i<vLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vLeptons.size(); j++)
        {
            if (  ( vLeptons.at(i).id() == -vLeptons.at(j).id()                               )
                    && ( fabs( ( vLeptons.at(i).p4() + vLeptons.at(j).p4() ).M() - 91.188) < Delta ) )
            {
                Lep1Z      = i;
                Lep2Z      = j;
                ZM         = ( vLeptons.at(i).p4() + vLeptons.at(j).p4() ).M();
                Zpt        = ( vLeptons.at(i).p4() + vLeptons.at(j).p4() ).Pt();
                Delta      = fabs(ZM - 91.188);
                pass_Zpeak = true;
            }
        }
    }
    if(!pass_Zpeak)       return;

    theHistoManager->fillHisto("CutFlow",                           "noSel", "WZrel_CR",   "",   3, weight);
    theHistoManager->fillHisto("CutFlow",                        "PassingZ", "WZrel_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass",        "PassingZ", "WZrel_CR",   "",  ZM, weight);

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vLeptons.size(); i++)
    {
        lepton_px = lepton_px + vLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                    "PassingMETLD", "WZrel_CR",   "",   1, weight);
    if( vSelectedJets.size() >= 4  )  theHistoManager->fillHisto("CutFlow",                     "PassingJets", "WZrel_CR",   "",   1, weight); 

    if( met_ld < 0.2 && vSelectedJets.size() < 4 ) return;
    theHistoManager->fillHisto("CutFlow",                 "PassingMETLDorJets", "WZrel_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass", "PassingMETLDorJets", "WZrel_CR",   "",  ZM, weight);

    // ########
    // # SFOS #
    // ########

    bool isSFOS = false;
    for(int i=0; i<vLeptons.size(); i++)
    {
        for(int j=0; j<vLeptons.size(); j++)
        {
            if (  ( i                   != j                            )
                    && ( vLeptons.at(i).id() == -vLeptons.at(j).id() ) )
            { isSFOS = true ;}
        }
    }

    theHistoManager->fillHisto("CutFlow",                 "PassingMETLDorSFOS", "WZrel_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass", "PassingMETLDorSFOS", "WZrel_CR",   "",  ZM, weight);

    if(vSelectedJets.size() < 4 && (met_ld < (0.2 + 0.1 * isSFOS)) ) return;

    // ##############
    // # b-jet veto #
    // ##############

    bool nMediumBtag      = ( nMediumBJets                              == 0 );
    if(!nMediumBtag)      return;

    theHistoManager->fillHisto("CutFlow",                 "PassingbJetsVeto", "WZrel_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass", "PassingbJetsVeto", "WZrel_CR",   "",  ZM, weight);

    int sum_charges_3l = 0;
    for(int i=0; i<3; i++)
    {
        sum_charges_3l = sum_charges_3l + vLeptons.at(i).charge();
    }
    if( fabs(sum_charges_3l) != 1 ) return;

    // Building arbitrary variables with all leptons...

    TLorentzVector all_lep_invmass_p4;

    int   all_lep_sumofcharges = 0;
    float all_lep_invmass      = 0;
    float all_lep_sumofpt      = 0;

    for(int i=0; i<vLeptons.size(); i++)
    {
        all_lep_sumofcharges = all_lep_sumofcharges + vLeptons.at(i).charge();
        all_lep_invmass_p4   = all_lep_invmass_p4   + vLeptons.at(i).p4();
    }

    all_lep_invmass = all_lep_invmass_p4.M();
    all_lep_sumofpt = sqrt( (lepton_px*lepton_px) + (lepton_py*lepton_py) );

    // Building arbitrary variables with remaining (after Z peak determination) leptons...

    TLorentzVector rem_lep_invmass_p4;

    int   rem_lep_sumofcharges = 0;
    float rem_lep_invmass      = 0;
    float rem_lep_sumofpt      = 0;
    float rem_lep_px           = 0;
    float rem_lep_py           = 0;


    for(int i=0; i<vLeptons.size(); i++)
    {
        if( i!=Lep1Z && i!=Lep2Z)
        {
            rem_lep_sumofcharges = rem_lep_sumofcharges + vLeptons.at(i).charge();
            rem_lep_px           = rem_lep_px           + vLeptons.at(i).p4().Px();
            rem_lep_py           = rem_lep_py           + vLeptons.at(i).p4().Py();
            rem_lep_invmass_p4   = rem_lep_invmass_p4   + vLeptons.at(i).p4();
        }
    }

    rem_lep_invmass = rem_lep_invmass_p4.M();
    rem_lep_sumofpt = sqrt( (rem_lep_px*rem_lep_px) + (rem_lep_py*rem_lep_py) );

    theHistoManager->fillHisto("ZCandidateInvariantMass",        "FinalCut", "WZrel_CR",   "",   ZM                    , weight);
    theHistoManager->fillHisto("ZCandidateTransverseMomentum",   "FinalCut", "WZrel_CR",   "",   Zpt                   , weight);
    theHistoManager->fillHisto("MET",                            "FinalCut", "WZrel_CR",   "",   vEvent->at(0).metpt() , weight);
    theHistoManager->fillHisto("MHT",                            "FinalCut", "WZrel_CR",   "",   MHT                   , weight);
    theHistoManager->fillHisto("MetLD",                          "FinalCut", "WZrel_CR",   "",   met_ld                , weight);
    theHistoManager->fillHisto("TauMultiplicity",                "FinalCut", "WZrel_CR",   "",   vSelectedTaus.size()  , weight);
    theHistoManager->fillHisto("JetMultiplicity",                "FinalCut", "WZrel_CR",   "",   vSelectedJets.size()  , weight);

    theHistoManager->fillHisto("InvariantMassOfSelectedLeptons", "FinalCut", "WZrel_CR",   "",   all_lep_invmass       , weight);
    theHistoManager->fillHisto("SumOfSelectedLeptonsCharges",    "FinalCut", "WZrel_CR",   "",   all_lep_sumofcharges  , weight);
    theHistoManager->fillHisto("SumVecPtSelectedLeptons",        "FinalCut", "WZrel_CR",   "",   all_lep_sumofpt       , weight);

    theHistoManager->fillHisto("InvMassRemainingLepton",         "FinalCut", "WZrel_CR",   "",   rem_lep_invmass       , weight);
    theHistoManager->fillHisto("SumOfRemainingLeptonsCharges",   "FinalCut", "WZrel_CR",   "",   rem_lep_sumofcharges  , weight);
    theHistoManager->fillHisto("SumVecPtRemainingLeptons",       "FinalCut", "WZrel_CR",   "",   rem_lep_sumofpt       , weight);

    theHistoManager->fillHisto2D("InvMassLastLeptonVSZMass",     "FinalCut", "WZrel_CR",   "",   rem_lep_invmass,    ZM, weight);
    theHistoManager->fillHisto2D("SumPtLepVSZMass",              "FinalCut", "WZrel_CR",   "",   rem_lep_sumofpt,    ZM, weight);
    theHistoManager->fillHisto2D("METLDVSZMass",                 "FinalCut", "WZrel_CR",   "",   met_ld,             ZM, weight);

    is_3l_WZrel_CR = true;
}

void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_TTZ(int evt)
{

    theHistoManager->fillHisto2D("LeptonsVsJets",           "noSel", "TTZ_CR",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("LeptonsVsBJets",          "noSel", "TTZ_CR",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection # with tight leptons !
    // ####################

    bool nLep             = ( vSelectedLeptons.size()                   >= 2 );
    if(!nLep)             return;

    bool leading_lep_pt   = ( vSelectedLeptons.at(0).pt()               > 20 );
    if(!leading_lep_pt)   return;

    bool following_lep_pt = ( vSelectedLeptons.at(1).pt()               > 10 );
    if(!following_lep_pt) return;

    bool passMll12Gt12    = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M()  > 12);
    if(!passMll12Gt12)    return;

    bool nJets            = ( vSelectedJets.size()                      >= 2 );
    if(!nJets)            return;

    // #################################
    // # Three leptons event selection #
    // #################################

    nLep                  = ( vSelectedLeptons.size()			>= 3 );
    if(!nLep)             return;

    bool third_lep_pt     = ( vSelectedLeptons.at(2).pt()               > 10 );
    if(!third_lep_pt)     return;

    passMll12Gt12  = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(2).p4() ).M()  > 12
            && ( vSelectedLeptons.at(1).p4() + vSelectedLeptons.at(2).p4() ).M()  > 12 );
    if(!passMll12Gt12)    return;

    theHistoManager->fillHisto("CutFlow",                           "noSel", "TTZ_CR",   "",   2, weight);
    theHistoManager->fillHisto("CutFlow",                 "PassingTightMVA", "TTZ_CR",   "",   1, weight);

    // ##########
    // # Z peak #
    // ##########

    bool pass_Zpeak = false;
    float Delta = 10, ZM = -1., Zpt = -1., Lep1Z = -1, Lep2Z = -1;

    for(int i=0; i<vSelectedLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vSelectedLeptons.size(); j++)
        {
            if (  ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id()                               )
                    && ( fabs( ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M() - 91.188) < Delta ) )
            {
                Lep1Z      = i;
                Lep2Z      = j;
                ZM         = ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M();
                Zpt        = ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).Pt();
                Delta      = fabs(ZM - 91.188);
                pass_Zpeak = true;
            }
        }
    }
    if(!pass_Zpeak)       return;

    theHistoManager->fillHisto("CutFlow",                           "noSel", "TTZ_CR",   "",   3, weight);
    theHistoManager->fillHisto("CutFlow",                        "PassingZ", "TTZ_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass",        "PassingZ", "TTZ_CR",   "",  ZM, weight);

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                    "PassingMETLD", "TTZ_CR",   "",   1, weight);
    if( vSelectedJets.size() >= 4  )  theHistoManager->fillHisto("CutFlow",                     "PassingJets", "TTZ_CR",   "",   1, weight); 

    if( met_ld < 0.2 && vSelectedJets.size() < 4 ) return;
    theHistoManager->fillHisto("CutFlow",                 "PassingMETLDorJets", "TTZ_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass", "PassingMETLDorJets", "TTZ_CR",   "",  ZM, weight);

    // ########
    // # SFOS #
    // ########

    bool isSFOS = false;
    for(int i=0; i<vSelectedLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vSelectedLeptons.size(); j++)
        {
            if (  ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id() ) )
            { isSFOS = true ;}
        }
    }

    theHistoManager->fillHisto("CutFlow",                 "PassingMETLDorSFOS", "TTZ_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass", "PassingMETLDorSFOS", "TTZ_CR",   "",  ZM, weight);

    if(vSelectedJets.size() < 4 && (met_ld < (0.2 + 0.1 * isSFOS)) ) return;

    // #####################
    // # b-jet requirement #
    // #####################

    bool nLooseBtag       = ( nLooseBJets                               >= 2 );
    bool nMediumBtag      = ( nMediumBJets                              >= 1 );
    if(!nLooseBtag || !nMediumBtag)      return;

    theHistoManager->fillHisto("CutFlow",                 "PassingbJetsVeto", "TTZ_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass", "PassingbJetsVeto", "TTZ_CR",   "",  ZM, weight);

    int sum_charges_3l = 0;
    for(int i=0; i<3; i++)
    {
        sum_charges_3l = sum_charges_3l + vSelectedLeptons.at(i).charge();
    }
    if( fabs(sum_charges_3l) != 1 ) return;

    // Building arbitrary variables with all leptons...

    TLorentzVector all_lep_invmass_p4;

    int   all_lep_sumofcharges = 0;
    float all_lep_invmass      = 0;
    float all_lep_sumofpt      = 0;

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        all_lep_sumofcharges = all_lep_sumofcharges + vSelectedLeptons.at(i).charge();
        all_lep_invmass_p4   = all_lep_invmass_p4   + vSelectedLeptons.at(i).p4();
    }

    all_lep_invmass = all_lep_invmass_p4.M();
    all_lep_sumofpt = sqrt( (lepton_px*lepton_px) + (lepton_py*lepton_py) );

    // Building arbitrary variables with remaining (after Z peak determination) leptons...

    TLorentzVector rem_lep_invmass_p4;

    int   rem_lep_sumofcharges = 0;
    float rem_lep_invmass      = 0;
    float rem_lep_sumofpt      = 0;
    float rem_lep_px           = 0;
    float rem_lep_py           = 0;


    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        if( i!=Lep1Z && i!=Lep2Z)
        {
            rem_lep_sumofcharges = rem_lep_sumofcharges + vSelectedLeptons.at(i).charge();
            rem_lep_px           = rem_lep_px           + vSelectedLeptons.at(i).p4().Px();
            rem_lep_py           = rem_lep_py           + vSelectedLeptons.at(i).p4().Py();
            rem_lep_invmass_p4   = rem_lep_invmass_p4   + vSelectedLeptons.at(i).p4();
        }
    }

    rem_lep_invmass = rem_lep_invmass_p4.M();
    rem_lep_sumofpt = sqrt( (rem_lep_px*rem_lep_px) + (rem_lep_py*rem_lep_py) );

    theHistoManager->fillHisto("LeadingLeptonPt",                "FinalCut", "TTZ_CR",   "",   vSelectedLeptons.at(0).pt()   , weight);
    theHistoManager->fillHisto("SubLeadingLeptonPt",             "FinalCut", "TTZ_CR",   "",   vSelectedLeptons.at(1).pt()   , weight);    

    theHistoManager->fillHisto("ZCandidateInvariantMass",        "FinalCut", "TTZ_CR",   "",   ZM                    , weight);
    theHistoManager->fillHisto("ZCandidateTransverseMomentum",   "FinalCut", "TTZ_CR",   "",   Zpt                   , weight);
    theHistoManager->fillHisto("MET",                            "FinalCut", "TTZ_CR",   "",   vEvent->at(0).metpt() , weight);
    theHistoManager->fillHisto("MHT",                            "FinalCut", "TTZ_CR",   "",   MHT                   , weight);
    theHistoManager->fillHisto("MetLD",                          "FinalCut", "TTZ_CR",   "",   met_ld                , weight);
    theHistoManager->fillHisto("TauMultiplicity",                "FinalCut", "TTZ_CR",   "",   vSelectedTaus.size()  , weight);
    theHistoManager->fillHisto("JetMultiplicity",                "FinalCut", "TTZ_CR",   "",   vSelectedJets.size()  , weight);

    theHistoManager->fillHisto("InvariantMassOfSelectedLeptons", "FinalCut", "TTZ_CR",   "",   all_lep_invmass       , weight);
    theHistoManager->fillHisto("SumOfSelectedLeptonsCharges",    "FinalCut", "TTZ_CR",   "",   all_lep_sumofcharges  , weight);
    theHistoManager->fillHisto("SumVecPtSelectedLeptons",        "FinalCut", "TTZ_CR",   "",   all_lep_sumofpt       , weight);

    theHistoManager->fillHisto("InvMassRemainingLepton",         "FinalCut", "TTZ_CR",   "",   rem_lep_invmass       , weight);
    theHistoManager->fillHisto("SumOfRemainingLeptonsCharges",   "FinalCut", "TTZ_CR",   "",   rem_lep_sumofcharges  , weight);
    theHistoManager->fillHisto("SumVecPtRemainingLeptons",       "FinalCut", "TTZ_CR",   "",   rem_lep_sumofpt       , weight);

    theHistoManager->fillHisto2D("InvMassLastLeptonVSZMass",     "FinalCut", "TTZ_CR",   "",   rem_lep_invmass,    ZM, weight);
    theHistoManager->fillHisto2D("SumPtLepVSZMass",              "FinalCut", "TTZ_CR",   "",   rem_lep_sumofpt,    ZM, weight);
    theHistoManager->fillHisto2D("METLDVSZMass",                 "FinalCut", "TTZ_CR",   "",   met_ld,             ZM, weight);

    if(vSelectedJets.size() >= 4)
    {
        theHistoManager->fillHisto("LeadingLeptonPt",                "FinalCut4j", "TTZ_CR",   "",   vSelectedLeptons.at(0).pt()  , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",             "FinalCut4j", "TTZ_CR",   "",   vSelectedLeptons.at(1).pt()  , weight);

        theHistoManager->fillHisto("MET",                            "FinalCut4j", "TTZ_CR",   "",   vEvent->at(0).metpt() , weight);
        theHistoManager->fillHisto("JetMultiplicity",                "FinalCut4j", "TTZ_CR",   "",   vSelectedJets.size()  , weight);
        theHistoManager->fillHisto("ZCandidateInvariantMass",        "FinalCut4j", "TTZ_CR",   "",   ZM                    , weight);
    }

    is_3l_TTZ_CR = true;
}

void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_CR_Zl(int evt)
{
    bool nLep     = ( vSelectedLeptons.size() == 2 ); 
    bool nLepFake = ( vFakeLeptons.size() == 1 );
    bool nJets    = ( vSelectedNonBTagJets.size() + vSelectedBTagJets.size() >= 2 );

    float MZ = -1.;

    bool pass_OSSF = false;
    if (nLep && vSelectedLeptons.at(0).id() == -vSelectedLeptons.at(1).id() &&  fabs( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M() - 91.188 ) < 10  )
    {
        MZ = ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M();
        pass_OSSF = true;
    }

    // common selection
    bool leading_lep_pt = 0;
    if (nLep) leading_lep_pt = ( vSelectedLeptons.at(0).pt() > 20 );
    bool following_lep_pt = 0;
    if (nLep) following_lep_pt = ( vSelectedLeptons.at(1).pt() > 10 );

    bool passMll12Gt12 = 0;
    if (nLep) passMll12Gt12  = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M()  > 12 );

    bool nLooseBtag  = ( nLooseBJets  >= 2 ); 

    /* std::cout << "nLep: "    << nLep
       << " nLepFake: "	     << nLepFake 
       << " nJets: "	     << nJets 
       << " pass_OSSF: "	     << pass_OSSF 
       << " leading_lep_pt: "   << leading_lep_pt 
       << " following_lep_pt: " << following_lep_pt
       << " passMll12Gt12: "    << passMll12Gt12 
       << " nLooseBtag: "       << nLooseBtag << std::endl; */


    if ( nLep && nLepFake && nJets && pass_OSSF && leading_lep_pt && following_lep_pt && passMll12Gt12 && nLooseBtag )
    {        
        is_Zl_CR = true;
        theHistoManager->fillHisto("ZCandidateInvariantMass", "CR_Zl", "", "", MZ, weight*weight_PV*mc_weight);
        //fillOutputTree();
    }

}

bool TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_TTH3l_MC() 
{ 
    bool sel_MC = true;

    //Check decays 
    if (!((vTruth->at(0).ttbar_decay()==1 && vTruth->at(0).boson_decay()==0) || //ttH
                (vTruth->at(0).ttbar_decay()==2 && vTruth->at(0).boson_decay()==1) || //ttH
                (vTruth->at(0).ttbar_decay()==1 && vTruth->at(0).boson_decay()==2) || //tt semi-lep, ttZ
                (vTruth->at(0).ttbar_decay()==2 && vTruth->at(0).boson_decay()==3)    //tt di-lep, ttW
         )) 
    { 
        sel_MC = false; 
        return sel_MC;}

        // 
        if (vTruth->at(0).JetsHighestPt_phi().size()<2 || vTruth->at(0).JetsClosestMw_phi().size()<2 || vTruth->at(0).JetsLowestMjj_phi().size()<2) 
        { 
            sel_MC = false; 
            return sel_MC;}

            //pt, eta of leptons	
            if (!(vTruth->at(0).Leptons_pt().at(0)> 10 && 
                        vTruth->at(0).Leptons_pt().at(1)> 10 &&
                        vTruth->at(0).Leptons_pt().at(2)> 10   )) sel_MC = false; 


            if (!(fabs(vTruth->at(0).Leptons_eta().at(0)) <2.5 && 
                        fabs(vTruth->at(0).Leptons_eta().at(1)) <2.5 &&
                        fabs(vTruth->at(0).Leptons_eta().at(2)) <2.5   )) sel_MC = false; 


            //lead. lepton
            if (!(vTruth->at(0).Leptons_pt().at(0) > 20 || 
                        vTruth->at(0).Leptons_pt().at(1) > 20 ||
                        vTruth->at(0).Leptons_pt().at(2) > 20  )) sel_MC = false; 


            //SFOS && M(ll) not in 81-101 ??? 
            int SFOSpair = -1;
            if ((vTruth->at(0).Leptons_id().at(0)==-vTruth->at(0).Leptons_id().at(1))) SFOSpair = 0;
            if ((vTruth->at(0).Leptons_id().at(1)==-vTruth->at(0).Leptons_id().at(2))) SFOSpair = 1;
            if ((vTruth->at(0).Leptons_id().at(0)==-vTruth->at(0).Leptons_id().at(2))) SFOSpair = 2;


            TLorentzVector Lep1;
            Lep1.SetPtEtaPhiE(vTruth->at(0).Leptons_pt().at(0),  vTruth->at(0).Leptons_eta().at(0), vTruth->at(0).Leptons_phi().at(0), vTruth->at(0).Leptons_E().at(0));
            TLorentzVector Lep2;
            Lep2.SetPtEtaPhiE(vTruth->at(0).Leptons_pt().at(1),  vTruth->at(0).Leptons_eta().at(1), vTruth->at(0).Leptons_phi().at(1), vTruth->at(0).Leptons_E().at(1));
            TLorentzVector Lep3;
            Lep3.SetPtEtaPhiE(vTruth->at(0).Leptons_pt().at(2),  vTruth->at(0).Leptons_eta().at(2), vTruth->at(0).Leptons_phi().at(2), vTruth->at(0).Leptons_E().at(2));


            if ( !(( Lep1+Lep2 ).M()  > 12 && ( Lep1+Lep3 ).M()  > 12 && ( Lep2+Lep3 ).M()  > 12 )) sel_MC = false;

            if ( (SFOSpair == 0 && fabs( (Lep1+Lep2 ).M()-91.188 ) < 10. ) || 
                    (SFOSpair == 1 && fabs( (Lep2+Lep3 ).M()-91.188 ) < 10. ) ||  
                    (SFOSpair == 2 && fabs( (Lep1+Lep3 ).M()-91.188 ) < 10. )    ) sel_MC = false;

            return sel_MC;

}

void TTbarHiggsMultileptonAnalysis::initializeOutputTree()
{

    outputfile->cd();
    tOutput = new TTree("Tree", "Tree");

    tOutput->Branch("mc_event",&mc_event,"mc_event/I");
    tOutput->Branch("weight",&weight,"weight/F");
    tOutput->Branch("mc_weight",&mc_weight,"mc_weight/F"); 
    tOutput->Branch("weight_scale_muF0p5",&weight_scale_muF0p5,"weight_scale_muF0p5/F");
    tOutput->Branch("weight_scale_muF2",&weight_scale_muF2,"weight_scale_muF2/F");
    tOutput->Branch("weight_scale_muR0p5",&weight_scale_muR0p5,"weight_scale_muR0p5/F");
    tOutput->Branch("weight_scale_muR2",&weight_scale_muR2,"weight_scale_muR2/F");
    tOutput->Branch("weight_csv_down",&weight_csv_down,"weight_csv_down/F");
    tOutput->Branch("weight_csv_up",&weight_csv_up,"weight_csv_up/F");

    tOutput->Branch("PV_weight",&weight_PV,"PV_weight/F");
    tOutput->Branch("mc_3l_category",&mc_3l_category,"mc_3l_category/I");
    tOutput->Branch("mc_ttbar_decay",&mc_ttbar_decay,"mc_ttbar_decay/I");
    tOutput->Branch("mc_boson_decay",&mc_boson_decay,"mc_boson_decay/I");
    tOutput->Branch("mc_ttZhypAllowed",&mc_ttZhypAllowed,"mc_ttZhypAllowed/I");
    tOutput->Branch("mc_nJets25",&mc_nJets25,"mc_nJets25/I");
    tOutput->Branch("mc_nBtagJets25",&mc_nBtagJets25,"mc_nBtagJets25/I");
    tOutput->Branch("mc_nMediumBtagJets25",&mc_nMediumBtagJets25,"mc_nMediumBtagJets25/I");
    tOutput->Branch("mc_nNonBtagJets25",&mc_nNonBtagJets25,"mc_nNonBtagJets25/I");

    tOutput->Branch("catJets",&catJets,"catJets/I");

    //tOutput->Branch("is_2lss_TTH_SR",&is_2lss_TTH_SR,"is_2lss_TTH_SR/B");
    //tOutput->Branch("is_2lss_JM_SB",&is_2lss_JM_SB,"is_2lss_JM_SB/B");
    //tOutput->Branch("is_2lss_LepMVA_SB",&is_2lss_LepMVA_SB,"is_2lss_LepMVA_SB/B");

    tOutput->Branch("is_2lss_TTH_SR",&is_2lss_TTH_SR,"is_2lss_TTH_SR/B");
    tOutput->Branch("is_3l_TTH_SR",&is_3l_TTH_SR,"is_3l_TTH_SR/B");

    tOutput->Branch("is_emu_TT_CR",&is_emu_TT_CR,"is_emu_TT_CR/B");
    //tOutput->Branch("is_3l_WZ_CR",&is_3l_WZ_CR,"is_3l_WZ_CR/B");
    tOutput->Branch("is_3l_TTZ_CR",&is_3l_TTZ_CR,"is_3l_TTZ_CR/B");
    tOutput->Branch("is_Zl_CR",&is_Zl_CR,"is_Zl_CR/B");

    tOutput->Branch("is_trigger",&is_trigger,"is_trigger/B");

    tOutput->Branch("max_Lep_eta", &max_Lep_eta, "max_Lep_eta/F");
    tOutput->Branch("MT_met_lep1",&MT_met_lep1,"MT_met_lep1/F");
    tOutput->Branch("nJet25_Recl",&nJet25_Recl,"nJet25_Recl/F");
    tOutput->Branch("mindr_lep1_jet",&mindr_lep1_jet,"mindr_lep1_jet/F");
    tOutput->Branch("mindr_lep2_jet",&mindr_lep2_jet,"mindr_lep2_jet/F");
    tOutput->Branch("LepGood_conePt0",&LepGood_conePt0,"LepGood_conePt0/F");
    tOutput->Branch("LepGood_conePt1",&LepGood_conePt1,"LepGood_conePt1/F");

    tOutput->Branch("signal_2lss_TT_MVA",&signal_2lss_TT_MVA,"signal_2lss_TT_MVA/F");
    tOutput->Branch("signal_2lss_TTV_MVA",&signal_2lss_TTV_MVA,"signal_2lss_TTV_MVA/F");

    tOutput->Branch("signal_3l_TT_MVA",&signal_3l_TT_MVA,"signal_3l_TT_MVA/F");
    tOutput->Branch("signal_3l_TTV_MVA",&signal_3l_TTV_MVA,"signal_3l_TTV_MVA/F");

    tOutput->Branch("multilepton_Lepton1_Id",&multilepton_Lepton1_Id,"multilepton_Lepton1_Id/I");
    tOutput->Branch("multilepton_Lepton1_P4","TLorentzVector",&multilepton_Lepton1_P4);
    tOutput->Branch("multilepton_Lepton2_Id",&multilepton_Lepton2_Id,"multilepton_Lepton2_Id/I");
    tOutput->Branch("multilepton_Lepton2_P4","TLorentzVector",&multilepton_Lepton2_P4);
    tOutput->Branch("multilepton_Lepton3_Id",&multilepton_Lepton3_Id,"multilepton_Lepton3_Id/I");
    tOutput->Branch("multilepton_Lepton3_P4","TLorentzVector",&multilepton_Lepton3_P4);
    tOutput->Branch("multilepton_Lepton4_Id",&multilepton_Lepton4_Id,"multilepton_Lepton4_Id/I");
    tOutput->Branch("multilepton_Lepton4_P4","TLorentzVector",&multilepton_Lepton4_P4);

    tOutput->Branch("multilepton_Bjet1_Id",&multilepton_Bjet1_Id,"multilepton_Bjet1_Id/I");
    tOutput->Branch("multilepton_Bjet1_P4","TLorentzVector",&multilepton_Bjet1_P4);
    tOutput->Branch("multilepton_Bjet1_CSV",&multilepton_Bjet1_CSV,"multilepton_Bjet1_CSV/F");
    tOutput->Branch("multilepton_Bjet1_JEC_Up",&multilepton_Bjet1_JEC_Up,"multilepton_Bjet1_JEC_Up/F");
    tOutput->Branch("multilepton_Bjet1_JEC_Down",&multilepton_Bjet1_JEC_Down,"multilepton_Bjet1_JEC_Down/F");
    tOutput->Branch("multilepton_Bjet1_JER_Up",&multilepton_Bjet1_JER_Up,"multilepton_Bjet1_JER_Up/F");
    tOutput->Branch("multilepton_Bjet1_JER_Down",&multilepton_Bjet1_JER_Down,"multilepton_Bjet1_JER_Down/F");

    tOutput->Branch("multilepton_Bjet2_Id",&multilepton_Bjet2_Id,"multilepton_Bjet2_Id/I");
    tOutput->Branch("multilepton_Bjet2_P4","TLorentzVector",&multilepton_Bjet2_P4);
    tOutput->Branch("multilepton_Bjet2_CSV",&multilepton_Bjet2_CSV,"multilepton_Bjet2_CSV/F");
    tOutput->Branch("multilepton_Bjet2_JEC_Up",&multilepton_Bjet2_JEC_Up,"multilepton_Bjet2_JEC_Up/F");
    tOutput->Branch("multilepton_Bjet2_JEC_Down",&multilepton_Bjet2_JEC_Down,"multilepton_Bjet2_JEC_Down/F");
    tOutput->Branch("multilepton_Bjet2_JER_Up",&multilepton_Bjet2_JER_Up,"multilepton_Bjet2_JER_Up/F");
    tOutput->Branch("multilepton_Bjet2_JER_Down",&multilepton_Bjet2_JER_Down,"multilepton_Bjet2_JER_Down/F");

    tOutput->Branch("multilepton_JetHighestPt1_Id",&multilepton_JetHighestPt1_Id,"multilepton_JetHighestPt1_Id/I");
    tOutput->Branch("multilepton_JetHighestPt1_P4","TLorentzVector",&multilepton_JetHighestPt1_P4);
    tOutput->Branch("multilepton_JetHighestPt1_CSV",&multilepton_JetHighestPt1_CSV,"multilepton_JetHighestPt1_CSV/F");
    tOutput->Branch("multilepton_JetHighestPt1_JEC_Up",&multilepton_JetHighestPt1_JEC_Up,"multilepton_JetHighestPt1_JEC_Up/F");
    tOutput->Branch("multilepton_JetHighestPt1_JEC_Down",&multilepton_JetHighestPt1_JEC_Down,"multilepton_JetHighestPt1_JEC_Down/F");
    tOutput->Branch("multilepton_JetHighestPt1_JER_Up",&multilepton_JetHighestPt1_JER_Up,"multilepton_JetHighestPt1_JER_Up/F");
    tOutput->Branch("multilepton_JetHighestPt1_JER_Down",&multilepton_JetHighestPt1_JER_Down,"multilepton_JetHighestPt1_JER_Down/F");

    tOutput->Branch("multilepton_JetHighestPt2_Id",&multilepton_JetHighestPt2_Id,"multilepton_JetHighestPt2_Id/I");
    tOutput->Branch("multilepton_JetHighestPt2_P4","TLorentzVector",&multilepton_JetHighestPt2_P4);
    tOutput->Branch("multilepton_JetHighestPt2_CSV",&multilepton_JetHighestPt2_CSV,"multilepton_JetHighestPt2_CSV/F");
    tOutput->Branch("multilepton_JetHighestPt2_JEC_Up",&multilepton_JetHighestPt2_JEC_Up,"multilepton_JetHighestPt2_JEC_Up/F");
    tOutput->Branch("multilepton_JetHighestPt2_JEC_Down",&multilepton_JetHighestPt2_JEC_Down,"multilepton_JetHighestPt2_JEC_Down/F");
    tOutput->Branch("multilepton_JetHighestPt2_JER_Up",&multilepton_JetHighestPt2_JER_Up,"multilepton_JetHighestPt2_JER_Up/F");
    tOutput->Branch("multilepton_JetHighestPt2_JER_Down",&multilepton_JetHighestPt2_JER_Down,"multilepton_JetHighestPt2_JER_Down/F");

    tOutput->Branch("multilepton_JetClosestMw1_Id",&multilepton_JetClosestMw1_Id,"multilepton_JetClosestMw1_Id/I");
    tOutput->Branch("multilepton_JetClosestMw1_P4","TLorentzVector",&multilepton_JetClosestMw1_P4);
    tOutput->Branch("multilepton_JetClosestMw1_CSV",&multilepton_JetClosestMw1_CSV,"multilepton_JetClosestMw1_CSV/F");
    tOutput->Branch("multilepton_JetClosestMw1_JEC_Up",&multilepton_JetClosestMw1_JEC_Up,"multilepton_JetClosestMw1_JEC_Up/F");
    tOutput->Branch("multilepton_JetClosestMw1_JEC_Down",&multilepton_JetClosestMw1_JEC_Down,"multilepton_JetClosestMw1_JEC_Down/F");
    tOutput->Branch("multilepton_JetClosestMw1_JER_Up",&multilepton_JetClosestMw1_JER_Up,"multilepton_JetClosestMw1_JER_Up/F");
    tOutput->Branch("multilepton_JetClosestMw1_JER_Down",&multilepton_JetClosestMw1_JER_Down,"multilepton_JetClosestMw1_JER_Down/F");

    tOutput->Branch("multilepton_JetClosestMw2_Id",&multilepton_JetClosestMw2_Id,"multilepton_JetClosestMw2_Id/I");
    tOutput->Branch("multilepton_JetClosestMw2_P4","TLorentzVector",&multilepton_JetClosestMw2_P4);
    tOutput->Branch("multilepton_JetClosestMw2_CSV",&multilepton_JetClosestMw2_CSV,"multilepton_JetClosestMw2_CSV/F");
    tOutput->Branch("multilepton_JetClosestMw2_JEC_Up",&multilepton_JetClosestMw2_JEC_Up,"multilepton_JetClosestMw2_JEC_Up/F");
    tOutput->Branch("multilepton_JetClosestMw2_JEC_Down",&multilepton_JetClosestMw2_JEC_Down,"multilepton_JetClosestMw2_JEC_Down/F");
    tOutput->Branch("multilepton_JetClosestMw2_JER_Up",&multilepton_JetClosestMw2_JER_Up,"multilepton_JetClosestMw2_JER_Up/F");
    tOutput->Branch("multilepton_JetClosestMw2_JER_Down",&multilepton_JetClosestMw2_JER_Down,"multilepton_JetClosestMw2_JER_Down/F");

    tOutput->Branch("multilepton_JetLowestMjj1_Id",&multilepton_JetLowestMjj1_Id,"multilepton_JetLowestMjj1_Id/I");
    tOutput->Branch("multilepton_JetLowestMjj1_P4","TLorentzVector",&multilepton_JetLowestMjj1_P4);
    tOutput->Branch("multilepton_JetLowestMjj1_CSV",&multilepton_JetLowestMjj1_CSV,"multilepton_JetLowestMjj1_CSV/F");
    tOutput->Branch("multilepton_JetLowestMjj1_JEC_Up",&multilepton_JetLowestMjj1_JEC_Up,"multilepton_JetLowestMjj1_JEC_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj1_JEC_Down",&multilepton_JetLowestMjj1_JEC_Down,"multilepton_JetLowestMjj1_JEC_Down/F");
    tOutput->Branch("multilepton_JetLowestMjj1_JER_Up",&multilepton_JetLowestMjj1_JER_Up,"multilepton_JetLowestMjj1_JER_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj1_JER_Down",&multilepton_JetLowestMjj1_JER_Down,"multilepton_JetLowestMjj1_JER_Down/F");

    tOutput->Branch("multilepton_JetLowestMjj2_Id",&multilepton_JetLowestMjj2_Id,"multilepton_JetLowestMjj2_Id/I");
    tOutput->Branch("multilepton_JetLowestMjj2_P4","TLorentzVector",&multilepton_JetLowestMjj2_P4);
    tOutput->Branch("multilepton_JetLowestMjj2_CSV",&multilepton_JetLowestMjj2_CSV,"multilepton_JetLowestMjj2_CSV/F");
    tOutput->Branch("multilepton_JetLowestMjj2_JEC_Up",&multilepton_JetLowestMjj2_JEC_Up,"multilepton_JetLowestMjj2_JEC_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj2_JEC_Down",&multilepton_JetLowestMjj2_JEC_Down,"multilepton_JetLowestMjj2_JEC_Down/F");
    tOutput->Branch("multilepton_JetLowestMjj2_JER_Up",&multilepton_JetLowestMjj2_JER_Up,"multilepton_JetLowestMjj2_JER_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj2_JER_Down",&multilepton_JetLowestMjj2_JER_Down,"multilepton_JetLowestMjj2_JER_Down/F");

    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_Id",&multilepton_JetHighestPt1_2ndPair_Id,"multilepton_JetHighestPt1_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_P4","TLorentzVector",&multilepton_JetHighestPt1_2ndPair_P4);
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_CSV",&multilepton_JetHighestPt1_2ndPair_CSV,"multilepton_JetHighestPt1_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_JEC_Up",&multilepton_JetHighestPt1_2ndPair_JEC_Up,"multilepton_JetHighestPt1_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_JEC_Down",&multilepton_JetHighestPt1_2ndPair_JEC_Down,"multilepton_JetHighestPt1_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_JER_Up",&multilepton_JetHighestPt1_2ndPair_JER_Up,"multilepton_JetHighestPt1_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_JER_Down",&multilepton_JetHighestPt1_2ndPair_JER_Down,"multilepton_JetHighestPt1_2ndPair_JER_Down/F");

    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_Id",&multilepton_JetHighestPt2_2ndPair_Id,"multilepton_JetHighestPt2_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_P4","TLorentzVector",&multilepton_JetHighestPt2_2ndPair_P4);
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_CSV",&multilepton_JetHighestPt2_2ndPair_CSV,"multilepton_JetHighestPt2_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_JEC_Up",&multilepton_JetHighestPt2_2ndPair_JEC_Up,"multilepton_JetHighestPt2_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_JEC_Down",&multilepton_JetHighestPt2_2ndPair_JEC_Down,"multilepton_JetHighestPt2_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_JER_Up",&multilepton_JetHighestPt2_2ndPair_JER_Up,"multilepton_JetHighestPt2_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_JER_Down",&multilepton_JetHighestPt2_2ndPair_JER_Down,"multilepton_JetHighestPt2_2ndPair_JER_Down/F");

    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_Id",&multilepton_JetClosestMw1_2ndPair_Id,"multilepton_JetClosestMw1_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_P4","TLorentzVector",&multilepton_JetClosestMw1_2ndPair_P4);
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_CSV",&multilepton_JetClosestMw1_2ndPair_CSV,"multilepton_JetClosestMw1_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_JEC_Up",&multilepton_JetClosestMw1_2ndPair_JEC_Up,"multilepton_JetClosestMw1_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_JEC_Down",&multilepton_JetClosestMw1_2ndPair_JEC_Down,"multilepton_JetClosestMw1_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_JER_Up",&multilepton_JetClosestMw1_2ndPair_JER_Up,"multilepton_JetClosestMw1_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_JER_Down",&multilepton_JetClosestMw1_2ndPair_JER_Down,"multilepton_JetClosestMw1_2ndPair_JER_Down/F");

    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_Id",&multilepton_JetClosestMw2_2ndPair_Id,"multilepton_JetClosestMw2_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_P4","TLorentzVector",&multilepton_JetClosestMw2_2ndPair_P4);
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_CSV",&multilepton_JetClosestMw2_2ndPair_CSV,"multilepton_JetClosestMw2_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_JEC_Up",&multilepton_JetClosestMw2_2ndPair_JEC_Up,"multilepton_JetClosestMw2_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_JEC_Down",&multilepton_JetClosestMw2_2ndPair_JEC_Down,"multilepton_JetClosestMw2_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_JER_Up",&multilepton_JetClosestMw2_2ndPair_JER_Up,"multilepton_JetClosestMw2_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_JER_Down",&multilepton_JetClosestMw2_2ndPair_JER_Down,"multilepton_JetClosestMw2_2ndPair_JER_Down/F");

    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_Id",&multilepton_JetLowestMjj1_2ndPair_Id,"multilepton_JetLowestMjj1_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_P4","TLorentzVector",&multilepton_JetLowestMjj1_2ndPair_P4);
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_CSV",&multilepton_JetLowestMjj1_2ndPair_CSV,"multilepton_JetLowestMjj1_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_JEC_Up",&multilepton_JetLowestMjj1_2ndPair_JEC_Up,"multilepton_JetLowestMjj1_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_JEC_Down",&multilepton_JetLowestMjj1_2ndPair_JEC_Down,"multilepton_JetLowestMjj1_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_JER_Up",&multilepton_JetLowestMjj1_2ndPair_JER_Up,"multilepton_JetLowestMjj1_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_JER_Down",&multilepton_JetLowestMjj1_2ndPair_JER_Down,"multilepton_JetLowestMjj1_2ndPair_JER_Down/F");

    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_Id",&multilepton_JetLowestMjj2_2ndPair_Id,"multilepton_JetLowestMjj2_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_P4","TLorentzVector",&multilepton_JetLowestMjj2_2ndPair_P4);
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_CSV",&multilepton_JetLowestMjj2_2ndPair_CSV,"multilepton_JetLowestMjj2_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_JEC_Up",&multilepton_JetLowestMjj2_2ndPair_JEC_Up,"multilepton_JetLowestMjj2_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_JEC_Down",&multilepton_JetLowestMjj2_2ndPair_JEC_Down,"multilepton_JetLowestMjj2_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_JER_Up",&multilepton_JetLowestMjj2_2ndPair_JER_Up,"multilepton_JetLowestMjj2_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_JER_Down",&multilepton_JetLowestMjj2_2ndPair_JER_Down,"multilepton_JetLowestMjj2_2ndPair_JER_Down/F");

    tOutput->Branch("multilepton_mET","TLorentzVector",&multilepton_mET);
    tOutput->Branch("multilepton_mETcov00",&multilepton_mETcov00,"multilepton_mETcov00/D");
    tOutput->Branch("multilepton_mETcov01",&multilepton_mETcov01,"multilepton_mETcov01/D");
    tOutput->Branch("multilepton_mETcov10",&multilepton_mETcov10,"multilepton_mETcov10/D");
    tOutput->Branch("multilepton_mETcov11",&multilepton_mETcov11,"multilepton_mETcov11/D");
    tOutput->Branch("multilepton_mHT",&multilepton_mHT,"multilepton_mHT/F");
    tOutput->Branch("multilepton_Ptot","TLorentzVector",&multilepton_Ptot);

    return;
}


void TTbarHiggsMultileptonAnalysis::fillOutputTree(){

    bool is2lss=false, is3l=false, is4l=false;

    int tot_charge = 0;
    int tot_id = 0;
    if (vSelectedLeptons.size()>=4) {
        for (unsigned int i=0; i<4; i++) {
            tot_charge += vSelectedLeptons.at(i).charge();
	    tot_id += vSelectedLeptons.at(i).id();
        }
    }
    if ((is_3l_TTH_SR || is_3l_TTZ_CR) && vSelectedLeptons.size()>=4 && tot_charge==0 && tot_id==0) is4l = true;
    else if ( ((is_3l_TTH_SR || is_3l_TTZ_CR) && vSelectedLeptons.size()==3)
	 || (is_Zl_CR && vSelectedLeptons.size() == 2 &&  vFakeLeptons.size() == 1) 
	 || ((is_3l_TTH_SR || is_3l_TTZ_CR) && vSelectedLeptons.size()>=4 && (tot_charge!=0 || tot_id!=0))) 
	is3l = true;
    else if ( is_2lss_TTH_SR && vSelectedLeptons.size()==2 && vSelectedLeptons.at(0).charge()==vSelectedLeptons.at(1).charge()) is2lss = true;
    if (!is2lss && !is3l && !is4l) return;

    if (vSelectedJets.size()<2) return;
    if (!(vSelectedBTagJets.size()>=2 || (vSelectedMediumBTagJets.size()==1))) return; 

    //if (vSelectedLeptons.size()<2) return; // 2lss only at the moment

    //std::cout << "lept="<<vSelectedLeptons.size()<<" fake="<<vFakeLeptons.size()<<std::endl;
    //std::cout << "btag="<<vSelectedBTagJets.size()<<" nonbtag="<<vSelectedNonBTagJets.size()<<std::endl;

    //Choosing 2 b-jets
    bool doSelectOnlyBjets = false;

    TLorentzVector Bjet1, Bjet2;
    int ib1=-1, ib2=-1;
    selectBjets("HighestBtagDiscrim", &ib1, &ib2, doSelectOnlyBjets);
    if (ib1!=-1) Bjet1.SetPtEtaPhiE(vSelectedJets.at(ib1).pt(), vSelectedJets.at(ib1).eta(), vSelectedJets.at(ib1).phi(), vSelectedJets.at(ib1).E());
    if (ib2!=-1) Bjet2.SetPtEtaPhiE(vSelectedJets.at(ib2).pt(), vSelectedJets.at(ib2).eta(), vSelectedJets.at(ib2).phi(), vSelectedJets.at(ib2).E());

    //2lss
    if (is2lss && ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2>=4) catJets = kCat_2lss_2b_4j;
    else if (is2lss && ib1!=-1 && ib2==-1 && vSelectedJets.size()-1>=4) catJets = kCat_2lss_1b_4j;
    else if (is2lss && ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2==3) catJets = kCat_2lss_2b_3j;
    else if (is2lss && ib1!=-1 && ib2==-1 && vSelectedJets.size()-1==3) catJets = kCat_2lss_1b_3j;
    else if (is2lss && ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2==2) catJets = kCat_2lss_2b_2j;
    //4l 
    else if (is4l && ib1!=-1 && ib2!=-1) catJets = kCat_4l_2b;
    else if (is4l && ib1!=-1 && ib2==-1) catJets = kCat_4l_1b;
    //3l
    else if (is3l && ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2>=2) catJets = kCat_3l_2b_2j;
    else if (is3l && ib1!=-1 && ib2==-1 && vSelectedJets.size()-1>=2) catJets = kCat_3l_1b_2j;
    else if (is3l && ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2==1) catJets = kCat_3l_2b_1j;
    else if (is3l && ib1!=-1 && ib2==-1 && vSelectedJets.size()-1==1) catJets = kCat_3l_1b_1j;
    else if (is3l && ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2==0) catJets = kCat_3l_2b_0j;
    else catJets = -1;

    //std::cout << "catJets="<<catJets<<std::endl;

    multilepton_Lepton1_Id = -999;
    multilepton_Lepton2_Id = -999;
    multilepton_Lepton3_Id = -999;
    multilepton_Lepton4_Id = -999;

    if (vSelectedLeptons.size()>=2){
        multilepton_Lepton1_P4 = vSelectedLeptons.at(0).p4();
        multilepton_Lepton1_Id = vSelectedLeptons.at(0).id();
        multilepton_Lepton2_P4 = vSelectedLeptons.at(1).p4();
        multilepton_Lepton2_Id = vSelectedLeptons.at(1).id();
    }

    if (vSelectedLeptons.size()>=3)
    {
        multilepton_Lepton3_P4 = vSelectedLeptons.at(2).p4();
        multilepton_Lepton3_Id = vSelectedLeptons.at(2).id();
    }
    else if (is_Zl_CR && vSelectedLeptons.size()==2 && vFakeLeptons.size()==1)
    {
        multilepton_Lepton3_P4 = vFakeLeptons.at(0).p4();
        multilepton_Lepton3_Id = vFakeLeptons.at(0).id();
    }

    if (vSelectedLeptons.size()>=4 && tot_charge==0)
    {
        multilepton_Lepton4_P4 = vSelectedLeptons.at(3).p4();
        multilepton_Lepton4_Id = vSelectedLeptons.at(3).id();
    }

    multilepton_Bjet1_Id = -999;
    if (ib1!=-1){
	FillJetInfoOutputTree(&multilepton_Bjet1_Id, 5, &multilepton_Bjet1_P4, Bjet1, &multilepton_Bjet1_CSV, vSelectedJets.at(ib1).CSVv2(), &multilepton_Bjet1_JEC_Up, &multilepton_Bjet1_JEC_Down, vSelectedJets.at(ib1).JES_uncert(), &multilepton_Bjet1_JER_Up, &multilepton_Bjet1_JER_Down, vSelectedJets.at(ib1).pt_JER(), vSelectedJets.at(ib1).pt_JER_up(), vSelectedJets.at(ib1).pt_JER_down());
        //multilepton_Bjet1_P4 = Bjet1;
        //multilepton_Bjet1_Id = 5;
    }
    multilepton_Bjet2_Id = -999;
    if (ib2!=-1){
	FillJetInfoOutputTree(&multilepton_Bjet2_Id, 5, &multilepton_Bjet2_P4, Bjet2, &multilepton_Bjet2_CSV, vSelectedJets.at(ib2).CSVv2(), &multilepton_Bjet2_JEC_Up, &multilepton_Bjet2_JEC_Down, vSelectedJets.at(ib2).JES_uncert(), &multilepton_Bjet2_JER_Up, &multilepton_Bjet2_JER_Down, vSelectedJets.at(ib2).pt_JER(), vSelectedJets.at(ib2).pt_JER_up(), vSelectedJets.at(ib2).pt_JER_down());
        //multilepton_Bjet2_P4 = Bjet2;
        //multilepton_Bjet2_Id = 5;
    }

    //Choose 2 jets
    TLorentzVector Pjet1, Pjet2;
    float pt_max=0, pt_max2=0; int ij1=-1, ij2=-1;
    float diffmass_min = 10000, mass_min = 10000; int ik1=-1, ik2=-1, il1=-1, il2=-1;
    for (unsigned int ij=0; ij<vSelectedJets.size(); ij++){
        if (ij==ib1 || ij==ib2) continue;
        if (vSelectedJets.at(ij).pt() > pt_max ) {
            pt_max2 = pt_max;
            ij2 = ij1;
            pt_max = vSelectedJets.at(ij).pt();
            ij1 = ij;
        } 
        if (vSelectedJets.at(ij).pt() < pt_max && vSelectedJets.at(ij).pt() > pt_max2){
            pt_max2 = vSelectedJets.at(ij).pt(); 
            ij2 = ij; 
        } 
        for (unsigned int ik=0; ik<vSelectedJets.size(); ik++){
            if (ik==ij) continue;
            if (ik==ib1 || ik==ib2) continue;
            Pjet1.SetPtEtaPhiE(vSelectedJets.at(ij).pt(), vSelectedJets.at(ij).eta(), vSelectedJets.at(ij).phi(), vSelectedJets.at(ij).E());
            Pjet2.SetPtEtaPhiE(vSelectedJets.at(ik).pt(), vSelectedJets.at(ik).eta(), vSelectedJets.at(ik).phi(), vSelectedJets.at(ik).E()); 
            if (TMath::Abs((Pjet1+Pjet2).M()-80.419)<diffmass_min){
                ik1=ij;
                ik2=ik;
                diffmass_min = TMath::Abs((Pjet1+Pjet2).M()-80.419);
            } 
            if ((Pjet1+Pjet2).M()<mass_min){
                il1=ij;
                il2=ik;
                mass_min = (Pjet1+Pjet2).M();
            } 
        } 
    }  

    //Choose 2 more jets
    int io1=-1, io2=-1, ip1=-1, ip2=-1, im1=-1, im2=-1;
    diffmass_min = 10000, mass_min = 10000, pt_max2 = 0, pt_max = 0;
    for (unsigned int im=0; im<vSelectedJets.size(); im++){
        if (im==ib1 || im==ib2 || im==ik1 || im==ik2) continue;
        if (vSelectedJets.at(im).pt() > pt_max ) {
            pt_max2 = pt_max;
            im2 = im1;
            pt_max = vSelectedJets.at(im).pt();
            im1 = im;
        }
        if (vSelectedJets.at(im).pt() < pt_max && vSelectedJets.at(im).pt() > pt_max2){
            pt_max2 = vSelectedJets.at(im).pt();
            im2 = im;
        }
        for (unsigned int in=0; in<vSelectedJets.size(); in++){
            if (in==ib1 || in==ib2 || in==ik1 || in==ik2 || in==im) continue;
            Pjet1.SetPtEtaPhiE(vSelectedJets.at(im).pt(), vSelectedJets.at(im).eta(), vSelectedJets.at(im).phi(), vSelectedJets.at(im).E());
            Pjet2.SetPtEtaPhiE(vSelectedJets.at(in).pt(), vSelectedJets.at(in).eta(), vSelectedJets.at(in).phi(), vSelectedJets.at(in).E());
            if (TMath::Abs((Pjet1+Pjet2).M()-80.419)<diffmass_min){
                io1=im;
                io2=in;
                diffmass_min = TMath::Abs((Pjet1+Pjet2).M()-80.419);
            }
            if ((Pjet1+Pjet2).M()<mass_min){
                ip1=im;
                ip2=in;
                mass_min = (Pjet1+Pjet2).M();
            }
        }
    }


    multilepton_JetHighestPt1_Id = -999;
    multilepton_JetHighestPt2_Id = -999;
    multilepton_JetClosestMw1_Id = -999;
    multilepton_JetClosestMw2_Id = -999;
    multilepton_JetLowestMjj1_Id = -999;
    multilepton_JetLowestMjj2_Id = -999;

    multilepton_JetHighestPt1_2ndPair_Id = -999;
    multilepton_JetHighestPt2_2ndPair_Id = -999;
    multilepton_JetClosestMw1_2ndPair_Id = -999;
    multilepton_JetClosestMw2_2ndPair_Id = -999;
    multilepton_JetLowestMjj1_2ndPair_Id = -999;
    multilepton_JetLowestMjj2_2ndPair_Id = -999;

    TLorentzVector Jet1, Jet2;

    if (ij1!=-1 && ij2==-1){
	Jet1.SetPtEtaPhiE(vSelectedJets.at(ij1).pt(), vSelectedJets.at(ij1).eta(), vSelectedJets.at(ij1).phi(), vSelectedJets.at(ij1).E());
	FillJetInfoOutputTree(&multilepton_JetHighestPt1_Id, 1, &multilepton_JetHighestPt1_P4, Jet1, &multilepton_JetHighestPt1_CSV, vSelectedJets.at(ij1).CSVv2(), &multilepton_JetHighestPt1_JEC_Up, &multilepton_JetHighestPt1_JEC_Down, vSelectedJets.at(ij1).JES_uncert(), &multilepton_JetHighestPt1_JER_Up, &multilepton_JetHighestPt1_JER_Down, vSelectedJets.at(ij1).pt_JER(), vSelectedJets.at(ij1).pt_JER_up(), vSelectedJets.at(ij1).pt_JER_down());
        //multilepton_JetHighestPt1_Id = 1;
        //multilepton_JetHighestPt1_P4.SetPtEtaPhiE(vSelectedJets.at(ij1).pt(), vSelectedJets.at(ij1).eta(), vSelectedJets.at(ij1).phi(), vSelectedJets.at(ij1).E());
    }   
    if (ij1!=-1 && ij2!=-1) {
	Jet1.SetPtEtaPhiE(vSelectedJets.at(ij1).pt(), vSelectedJets.at(ij1).eta(), vSelectedJets.at(ij1).phi(), vSelectedJets.at(ij1).E());
	Jet2.SetPtEtaPhiE(vSelectedJets.at(ij2).pt(), vSelectedJets.at(ij2).eta(), vSelectedJets.at(ij2).phi(), vSelectedJets.at(ij2).E());
	FillJetInfoOutputTree(&multilepton_JetHighestPt1_Id, 1, &multilepton_JetHighestPt1_P4, Jet1, &multilepton_JetHighestPt1_CSV, vSelectedJets.at(ij1).CSVv2(), &multilepton_JetHighestPt1_JEC_Up, &multilepton_JetHighestPt1_JEC_Down, vSelectedJets.at(ij1).JES_uncert(), &multilepton_JetHighestPt1_JER_Up, &multilepton_JetHighestPt1_JER_Down, vSelectedJets.at(ij1).pt_JER(), vSelectedJets.at(ij1).pt_JER_up(), vSelectedJets.at(ij1).pt_JER_down());
	FillJetInfoOutputTree(&multilepton_JetHighestPt2_Id, 1, &multilepton_JetHighestPt2_P4, Jet2, &multilepton_JetHighestPt2_CSV, vSelectedJets.at(ij2).CSVv2(), &multilepton_JetHighestPt2_JEC_Up, &multilepton_JetHighestPt2_JEC_Down, vSelectedJets.at(ij2).JES_uncert(), &multilepton_JetHighestPt2_JER_Up, &multilepton_JetHighestPt2_JER_Down, vSelectedJets.at(ij2).pt_JER(), vSelectedJets.at(ij2).pt_JER_up(), vSelectedJets.at(ij2).pt_JER_down());	
        //multilepton_JetHighestPt1_Id = 1;
        //multilepton_JetHighestPt2_Id = 1;
        //multilepton_JetHighestPt1_P4.SetPtEtaPhiE(vSelectedJets.at(ij1).pt(), vSelectedJets.at(ij1).eta(), vSelectedJets.at(ij1).phi(), vSelectedJets.at(ij1).E());
        //multilepton_JetHighestPt2_P4.SetPtEtaPhiE(vSelectedJets.at(ij2).pt(), vSelectedJets.at(ij2).eta(), vSelectedJets.at(ij2).phi(), vSelectedJets.at(ij2).E());
    }
    if (ik1!=-1 && ik2!=-1){
	Jet1.SetPtEtaPhiE(vSelectedJets.at(ik1).pt(), vSelectedJets.at(ik1).eta(), vSelectedJets.at(ik1).phi(), vSelectedJets.at(ik1).E());
	Jet2.SetPtEtaPhiE(vSelectedJets.at(ik2).pt(), vSelectedJets.at(ik2).eta(), vSelectedJets.at(ik2).phi(), vSelectedJets.at(ik2).E());
	FillJetInfoOutputTree(&multilepton_JetClosestMw1_Id, 2, &multilepton_JetClosestMw1_P4, Jet1, &multilepton_JetClosestMw1_CSV, vSelectedJets.at(ik1).CSVv2(), &multilepton_JetClosestMw1_JEC_Up, &multilepton_JetClosestMw1_JEC_Down, vSelectedJets.at(ik1).JES_uncert(), &multilepton_JetClosestMw1_JER_Up, &multilepton_JetClosestMw1_JER_Down, vSelectedJets.at(ik1).pt_JER(), vSelectedJets.at(ik1).pt_JER_up(), vSelectedJets.at(ik1).pt_JER_down());
        FillJetInfoOutputTree(&multilepton_JetClosestMw2_Id, 2, &multilepton_JetClosestMw2_P4, Jet2, &multilepton_JetClosestMw2_CSV, vSelectedJets.at(ik2).CSVv2(), &multilepton_JetClosestMw2_JEC_Up, &multilepton_JetClosestMw2_JEC_Down, vSelectedJets.at(ik2).JES_uncert(), &multilepton_JetClosestMw2_JER_Up, &multilepton_JetClosestMw2_JER_Down, vSelectedJets.at(ik2).pt_JER(), vSelectedJets.at(ik2).pt_JER_up(), vSelectedJets.at(ik2).pt_JER_down());
        //multilepton_JetClosestMw1_Id = 2;
        //multilepton_JetClosestMw2_Id = 2;
        //multilepton_JetClosestMw1_P4.SetPtEtaPhiE(vSelectedJets.at(ik1).pt(), vSelectedJets.at(ik1).eta(), vSelectedJets.at(ik1).phi(), vSelectedJets.at(ik1).E());
        //multilepton_JetClosestMw2_P4.SetPtEtaPhiE(vSelectedJets.at(ik2).pt(), vSelectedJets.at(ik2).eta(), vSelectedJets.at(ik2).phi(), vSelectedJets.at(ik2).E());
    }
    if (il1!=-1 && il2!=-1){
        Jet1.SetPtEtaPhiE(vSelectedJets.at(il1).pt(), vSelectedJets.at(il1).eta(), vSelectedJets.at(il1).phi(), vSelectedJets.at(il1).E());
	Jet2.SetPtEtaPhiE(vSelectedJets.at(il2).pt(), vSelectedJets.at(il2).eta(), vSelectedJets.at(il2).phi(), vSelectedJets.at(il2).E());
	FillJetInfoOutputTree(&multilepton_JetLowestMjj1_Id, 3, &multilepton_JetLowestMjj1_P4, Jet1, &multilepton_JetLowestMjj1_CSV, vSelectedJets.at(il1).CSVv2(), &multilepton_JetLowestMjj1_JEC_Up, &multilepton_JetLowestMjj1_JEC_Down, vSelectedJets.at(il1).JES_uncert(), &multilepton_JetLowestMjj1_JER_Up, &multilepton_JetLowestMjj1_JER_Down, vSelectedJets.at(il1).pt_JER(), vSelectedJets.at(il1).pt_JER_up(), vSelectedJets.at(il1).pt_JER_down());
	FillJetInfoOutputTree(&multilepton_JetLowestMjj2_Id, 3, &multilepton_JetLowestMjj2_P4, Jet2, &multilepton_JetLowestMjj2_CSV, vSelectedJets.at(il2).CSVv2(), &multilepton_JetLowestMjj2_JEC_Up, &multilepton_JetLowestMjj2_JEC_Down, vSelectedJets.at(il2).JES_uncert(), &multilepton_JetLowestMjj2_JER_Up, &multilepton_JetLowestMjj2_JER_Down, vSelectedJets.at(il2).pt_JER(), vSelectedJets.at(il2).pt_JER_up(), vSelectedJets.at(il2).pt_JER_down());
        //multilepton_JetLowestMjj1_Id = 3;
        //multilepton_JetLowestMjj2_Id = 3;
        //multilepton_JetLowestMjj1_P4.SetPtEtaPhiE(vSelectedJets.at(il1).pt(), vSelectedJets.at(il1).eta(), vSelectedJets.at(il1).phi(), vSelectedJets.at(il1).E());
        //multilepton_JetLowestMjj2_P4.SetPtEtaPhiE(vSelectedJets.at(il2).pt(), vSelectedJets.at(il2).eta(), vSelectedJets.at(il2).phi(), vSelectedJets.at(il2).E());
    }

    //2nd pair (first one: closest to Mw)
    if (is2lss && ij1!=-1 && ij2!=-1){
        if (im1!=-1 && im2==-1){
	    Jet1.SetPtEtaPhiE(vSelectedJets.at(im1).pt(), vSelectedJets.at(im1).eta(), vSelectedJets.at(im1).phi(), vSelectedJets.at(im1).E());
	    FillJetInfoOutputTree(&multilepton_JetHighestPt1_2ndPair_Id, 1, &multilepton_JetHighestPt1_2ndPair_P4, Jet1, &multilepton_JetHighestPt1_2ndPair_CSV, vSelectedJets.at(im1).CSVv2(), &multilepton_JetHighestPt1_2ndPair_JEC_Up, &multilepton_JetHighestPt1_2ndPair_JEC_Down, vSelectedJets.at(im1).JES_uncert(), &multilepton_JetHighestPt1_2ndPair_JER_Up, &multilepton_JetHighestPt1_2ndPair_JER_Down, vSelectedJets.at(im1).pt_JER(), vSelectedJets.at(im1).pt_JER_up(), vSelectedJets.at(im1).pt_JER_down());
            //multilepton_JetHighestPt1_2ndPair_Id = 1;
            //multilepton_JetHighestPt1_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(im1).pt(), vSelectedJets.at(im1).eta(), vSelectedJets.at(im1).phi(), vSelectedJets.at(im1).E());
        }
        if (im1!=-1 && im2!=-1){
            Jet1.SetPtEtaPhiE(vSelectedJets.at(im1).pt(), vSelectedJets.at(im1).eta(), vSelectedJets.at(im1).phi(), vSelectedJets.at(im1).E());
            Jet2.SetPtEtaPhiE(vSelectedJets.at(im2).pt(), vSelectedJets.at(im2).eta(), vSelectedJets.at(im2).phi(), vSelectedJets.at(im2).E());
	    FillJetInfoOutputTree(&multilepton_JetHighestPt1_2ndPair_Id, 1, &multilepton_JetHighestPt1_2ndPair_P4, Jet1, &multilepton_JetHighestPt1_2ndPair_CSV, vSelectedJets.at(im1).CSVv2(), &multilepton_JetHighestPt1_2ndPair_JEC_Up, &multilepton_JetHighestPt1_2ndPair_JEC_Down, vSelectedJets.at(im1).JES_uncert(), &multilepton_JetHighestPt1_2ndPair_JER_Up, &multilepton_JetHighestPt1_2ndPair_JER_Down, vSelectedJets.at(im1).pt_JER(), vSelectedJets.at(im1).pt_JER_up(), vSelectedJets.at(im1).pt_JER_down());
	    FillJetInfoOutputTree(&multilepton_JetHighestPt2_2ndPair_Id, 1, &multilepton_JetHighestPt2_2ndPair_P4, Jet2, &multilepton_JetHighestPt2_2ndPair_CSV, vSelectedJets.at(im2).CSVv2(), &multilepton_JetHighestPt2_2ndPair_JEC_Up, &multilepton_JetHighestPt2_2ndPair_JEC_Down, vSelectedJets.at(im2).JES_uncert(), &multilepton_JetHighestPt2_2ndPair_JER_Up, &multilepton_JetHighestPt2_2ndPair_JER_Down, vSelectedJets.at(im2).pt_JER(), vSelectedJets.at(im2).pt_JER_up(), vSelectedJets.at(im2).pt_JER_down());
            //multilepton_JetHighestPt1_2ndPair_Id = 1;
            //multilepton_JetHighestPt1_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(im1).pt(), vSelectedJets.at(im1).eta(), vSelectedJets.at(im1).phi(), vSelectedJets.at(im1).E());
            //multilepton_JetHighestPt2_2ndPair_Id = 1;
            //multilepton_JetHighestPt2_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(im2).pt(), vSelectedJets.at(im2).eta(), vSelectedJets.at(im2).phi(), vSelectedJets.at(im2).E());   
        }
        if (io1!=-1 && io2!=-1){
	    Jet1.SetPtEtaPhiE(vSelectedJets.at(ip1).pt(), vSelectedJets.at(ip1).eta(), vSelectedJets.at(ip1).phi(), vSelectedJets.at(ip1).E());
	    Jet2.SetPtEtaPhiE(vSelectedJets.at(io2).pt(), vSelectedJets.at(io2).eta(), vSelectedJets.at(io2).phi(), vSelectedJets.at(io2).E());
	    FillJetInfoOutputTree(&multilepton_JetClosestMw1_2ndPair_Id, 2, &multilepton_JetClosestMw1_2ndPair_P4, Jet1, &multilepton_JetClosestMw1_2ndPair_CSV, vSelectedJets.at(io1).CSVv2(), &multilepton_JetClosestMw1_2ndPair_JEC_Up, &multilepton_JetClosestMw1_2ndPair_JEC_Down, vSelectedJets.at(io1).JES_uncert(), &multilepton_JetClosestMw1_2ndPair_JER_Up, &multilepton_JetClosestMw1_2ndPair_JER_Down, vSelectedJets.at(io1).pt_JER(), vSelectedJets.at(io1).pt_JER_up(), vSelectedJets.at(io1).pt_JER_down());
            FillJetInfoOutputTree(&multilepton_JetClosestMw2_2ndPair_Id, 2, &multilepton_JetClosestMw2_2ndPair_P4, Jet2, &multilepton_JetClosestMw2_2ndPair_CSV, vSelectedJets.at(io2).CSVv2(), &multilepton_JetClosestMw2_2ndPair_JEC_Up, &multilepton_JetClosestMw2_2ndPair_JEC_Down, vSelectedJets.at(io2).JES_uncert(), &multilepton_JetClosestMw2_2ndPair_JER_Up, &multilepton_JetClosestMw2_2ndPair_JER_Down, vSelectedJets.at(io2).pt_JER(), vSelectedJets.at(io2).pt_JER_up(), vSelectedJets.at(io2).pt_JER_down());
            //multilepton_JetClosestMw1_2ndPair_Id = 2;
            //multilepton_JetClosestMw2_2ndPair_Id = 2;
            //multilepton_JetClosestMw1_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(io1).pt(), vSelectedJets.at(io1).eta(), vSelectedJets.at(io1).phi(), vSelectedJets.at(io1).E());
            //multilepton_JetClosestMw2_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(io2).pt(), vSelectedJets.at(io2).eta(), vSelectedJets.at(io2).phi(), vSelectedJets.at(io2).E());
        }
        if (ip1!=-1 && ip2!=-1){
	    Jet1.SetPtEtaPhiE(vSelectedJets.at(ip1).pt(), vSelectedJets.at(ip1).eta(), vSelectedJets.at(ip1).phi(), vSelectedJets.at(ip1).E());
	    Jet2.SetPtEtaPhiE(vSelectedJets.at(ip2).pt(), vSelectedJets.at(ip2).eta(), vSelectedJets.at(ip2).phi(), vSelectedJets.at(ip2).E());
	    FillJetInfoOutputTree(&multilepton_JetLowestMjj1_2ndPair_Id, 3, &multilepton_JetLowestMjj1_2ndPair_P4, Jet1, &multilepton_JetLowestMjj1_2ndPair_CSV, vSelectedJets.at(ip1).CSVv2(), &multilepton_JetLowestMjj1_2ndPair_JEC_Up, &multilepton_JetLowestMjj1_2ndPair_JEC_Down, vSelectedJets.at(ip1).JES_uncert(), &multilepton_JetLowestMjj1_2ndPair_JER_Up, &multilepton_JetLowestMjj1_2ndPair_JER_Down, vSelectedJets.at(ip1).pt_JER(), vSelectedJets.at(ip1).pt_JER_up(), vSelectedJets.at(ip1).pt_JER_down());
            FillJetInfoOutputTree(&multilepton_JetLowestMjj2_2ndPair_Id, 3, &multilepton_JetLowestMjj2_2ndPair_P4, Jet2, &multilepton_JetLowestMjj2_2ndPair_CSV, vSelectedJets.at(ip2).CSVv2(), &multilepton_JetLowestMjj2_2ndPair_JEC_Up, &multilepton_JetLowestMjj2_2ndPair_JEC_Down, vSelectedJets.at(ip2).JES_uncert(), &multilepton_JetLowestMjj2_2ndPair_JER_Up, &multilepton_JetLowestMjj2_2ndPair_JER_Down, vSelectedJets.at(ip2).pt_JER(), vSelectedJets.at(ip2).pt_JER_up(), vSelectedJets.at(ip2).pt_JER_down());
            //multilepton_JetLowestMjj1_2ndPair_Id = 3;
            //multilepton_JetLowestMjj2_2ndPair_Id = 3;
            //multilepton_JetLowestMjj1_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(ip1).pt(), vSelectedJets.at(ip1).eta(), vSelectedJets.at(ip1).phi(), vSelectedJets.at(ip1).E());
            //multilepton_JetLowestMjj2_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(ip2).pt(), vSelectedJets.at(ip2).eta(), vSelectedJets.at(ip2).phi(), vSelectedJets.at(ip2).E());
        }
    }

    multilepton_mET.SetPtEtaPhiE(vEvent->at(0).metpt(), 0, vEvent->at(0).metphi(), vEvent->at(0).metpt());
    multilepton_mETcov00 = vEvent->at(0).metcov00();
    multilepton_mETcov01 = vEvent->at(0).metcov01();
    multilepton_mETcov10 = vEvent->at(0).metcov10();
    multilepton_mETcov11 = vEvent->at(0).metcov11();
    multilepton_mHT = vEvent->at(0).metsumet();

    mc_ttZhypAllowed = 0;
/*
    if(vSelectedLeptons.size()==3) 
    {
        if ( vSelectedLeptons.at(0).charge()==vSelectedLeptons.at(1).charge() && vSelectedLeptons.at(1).charge()==vSelectedLeptons.at(2).charge() ) mc_ttZhypAllowed =-1;
        else if (  ( vSelectedLeptons.at(0).id() == -vSelectedLeptons.at(1).id() ) 
                || ( vSelectedLeptons.at(0).id() == -vSelectedLeptons.at(2).id() ) 
                || ( vSelectedLeptons.at(1).id() == -vSelectedLeptons.at(2).id() ))
            mc_ttZhypAllowed = 1; }
*/

    if (multilepton_Lepton1_Id!=-999 && multilepton_Lepton2_Id!=-999 && multilepton_Lepton3_Id!=-999){
	if (multilepton_Lepton1_Id*multilepton_Lepton2_Id>0 && multilepton_Lepton2_Id*multilepton_Lepton3_Id>0) mc_ttZhypAllowed =-1;
	else if ( (multilepton_Lepton1_Id==-multilepton_Lepton2_Id)
		|| (multilepton_Lepton1_Id==-multilepton_Lepton3_Id)
		|| (multilepton_Lepton2_Id==-multilepton_Lepton3_Id))
			mc_ttZhypAllowed = 1; }


        mc_nJets25 = vSelectedJets.size();
        mc_nBtagJets25 = vSelectedBTagJets.size();
        mc_nMediumBtagJets25 = vSelectedMediumBTagJets.size();
        mc_nNonBtagJets25 = vSelectedNonBTagJets.size();

        tOutput->Fill();

        //if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);

}

void TTbarHiggsMultileptonAnalysis::FillJetInfoOutputTree(int* tree_Id, int Id, TLorentzVector* tree_P4, TLorentzVector P4, float* tree_CSV, float CSV, float* tree_JEC_Up, float* tree_JEC_Down, float JEC_value, float* tree_JER_Up, float* tree_JER_Down, float JER, float JER_Up, float JER_Down){

  *tree_Id = Id;
  *tree_P4 = P4;
  
  *tree_CSV = CSV;

  *tree_JEC_Up = P4.E()*(1.+JEC_value);
  *tree_JEC_Down = P4.E()*(1.-JEC_value);

  *tree_JER_Up = P4.E()*JER_Up/JER;
  *tree_JER_Down = P4.E()*JER_Down/JER;

  return;
}

void TTbarHiggsMultileptonAnalysis::PrintLHCOforMadweight_MC(int evt)
{  

    //std::cout <<" _MC 1"<< std::endl;

    /*if( proc<-1 || proc > 6 )f
      {
      std::cout << "proc can only take following values: -1,1,2,3,4,5,6" << std::endl;
      std::cout << "3l final state specific" << std::endl;    
      std::cout << "1,2,3,4: specific to the ttH final state, cf patches provided to madweight" << std::endl;
      std::cout << "5: specific to the ttZ with l+l-l+ final state" << std::endl;
      std::cout << "6: specific to the ttZ with l+l-l- final state" << std::endl;
      std::cout << "-1: no selection on the final state applied" << std::endl;

      return;
      }*/


    fline0 = "0      " + std::string(Form("%d",evt)) + "     " + trig;

    int nobj = 1;	      

    int multilepton_Lepton1_Id_LHCO = -666;
    int multilepton_Lepton2_Id_LHCO = -666;
    int multilepton_Lepton3_Id_LHCO = -666;

    //
    // LHCO lepton ID convention
    //  
    if (abs(vTruth->at(0).Leptons_id().at(0))==11) multilepton_Lepton1_Id_LHCO = 1 ;
    else if (abs(vTruth->at(0).Leptons_id().at(0))==13) multilepton_Lepton1_Id_LHCO = 2 ;
    if (abs(vTruth->at(0).Leptons_id().at(1))==11) multilepton_Lepton2_Id_LHCO = 1 ;
    else if (abs(vTruth->at(0).Leptons_id().at(1))==13) multilepton_Lepton2_Id_LHCO = 2 ;
    if (abs(vTruth->at(0).Leptons_id().at(2))==11) multilepton_Lepton3_Id_LHCO = 1 ;
    else if (abs(vTruth->at(0).Leptons_id().at(2))==13) multilepton_Lepton3_Id_LHCO = 2 ;

    //
    // LHCO phi convention
    //
    float multilepton_Lepton1_phi       = Phi_0_2Pi(vTruth->at(0).Leptons_phi().at(0));
    float multilepton_Lepton2_phi       = Phi_0_2Pi(vTruth->at(0).Leptons_phi().at(1));
    float multilepton_Lepton3_phi       = Phi_0_2Pi(vTruth->at(0).Leptons_phi().at(2));
    float multilepton_Bjet1_phi         = Phi_0_2Pi(vTruth->at(0).Bjets_phi().at(0));
    float multilepton_Bjet2_phi         = Phi_0_2Pi(vTruth->at(0).Bjets_phi().at(1));	  
    float multilepton_JetHighestPt1_phi = Phi_0_2Pi(vTruth->at(0).JetsHighestPt_phi().at(0));
    float multilepton_JetHighestPt2_phi = Phi_0_2Pi(vTruth->at(0).JetsHighestPt_phi().at(1));
    float multilepton_JetClosestMw1_phi = Phi_0_2Pi(vTruth->at(0).JetsClosestMw_phi().at(0));
    float multilepton_JetClosestMw2_phi = Phi_0_2Pi(vTruth->at(0).JetsClosestMw_phi().at(1));
    float multilepton_JetLowestMjj1_phi = Phi_0_2Pi(vTruth->at(0).JetsLowestMjj_phi().at(0));
    float multilepton_JetLowestMjj2_phi = Phi_0_2Pi(vTruth->at(0).JetsLowestMjj_phi().at(1));

    //std::cout <<" _MC 22"<< std::endl;

    // l1
    std::string l1_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d", 
                nobj,multilepton_Lepton1_Id_LHCO,vTruth->at(0).Leptons_eta().at(0),multilepton_Lepton1_phi,vTruth->at(0).Leptons_pt().at(0),0.0,vTruth->at(0).Leptons_id().at(0)/abs(vTruth->at(0).Leptons_id().at(0)),0,0,0,0));
    nobj++; 
    //std::cout <<" _MC 23"<< std::endl;
    // l2
    std::string l2_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d",
                nobj,multilepton_Lepton2_Id_LHCO,vTruth->at(0).Leptons_eta().at(1),multilepton_Lepton2_phi,vTruth->at(1).Leptons_pt().at(1),0.0,vTruth->at(0).Leptons_id().at(1)/abs(vTruth->at(0).Leptons_id().at(1)),0,0,0,0));
    nobj++;	
    //std::cout <<" _MC 24"<< std::endl;
    // l3
    std::string l3_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d",
                nobj,multilepton_Lepton3_Id_LHCO,vTruth->at(0).Leptons_eta().at(2),multilepton_Lepton3_phi,vTruth->at(2).Leptons_pt().at(2),0.0,vTruth->at(0).Leptons_id().at(2)/abs(vTruth->at(0).Leptons_id().at(2)),0,0,0,0));
    nobj++;
    //std::cout <<" _MC 3"<< std::endl;

    //										    
    std::string j1_fline;
    std::string j2_fline;

    if ( _processLHCO_MC == 5 || _processLHCO_MC == 6 || _processLHCO_MC == 4 || _processLHCO_MC == 3 )
    {			     
        // j1
        j1_fline = std::string(Form("%d	  %d	 %.2f	    %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
                    nobj,4,vTruth->at(0).JetsClosestMw_eta().at(0),multilepton_JetClosestMw1_phi,vTruth->at(0).JetsClosestMw_pt().at(0),0.0,1,0,0,0,0));
        nobj++;

        // j2
        j2_fline = std::string(Form("%d	  %d	  %.2f      %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
                    nobj,4,vTruth->at(0).JetsClosestMw_eta().at(1),multilepton_JetClosestMw2_phi,vTruth->at(0).JetsClosestMw_pt().at(1),0.0,1,0,0,0,0));
        nobj++;
    }
    else if ( _processLHCO_MC == 1 || _processLHCO_MC == 2)
    {
        // j1
        j1_fline = std::string(Form("%d	  %d	  %.2f      %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
                    nobj,4,vTruth->at(0).JetsLowestMjj_eta().at(0),multilepton_JetLowestMjj1_phi,vTruth->at(0).JetsLowestMjj_pt().at(0),0.0,1,0,0,0,0));
        nobj++;

        // j2
        j2_fline = std::string(Form("%d	  %d	   %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
                    nobj,4,vTruth->at(0).JetsLowestMjj_eta().at(1),multilepton_JetLowestMjj2_phi,vTruth->at(0).JetsLowestMjj_pt().at(1),0.0,1,0,0,0,0));
        nobj++;
    }

    //
    TLorentzVector BJet1;
    BJet1.SetPtEtaPhiE(vTruth->at(0).Bjets_pt().at(0), vTruth->at(0).Bjets_eta().at(0), vTruth->at(0).Bjets_phi().at(0), vTruth->at(0).Bjets_E().at(0));

    TLorentzVector BJet2;
    BJet2.SetPtEtaPhiE(vTruth->at(0).Bjets_pt().at(1), vTruth->at(0).Bjets_eta().at(1), vTruth->at(0).Bjets_phi().at(1), vTruth->at(0).Bjets_E().at(1));


    // bj1
    std::string b1_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d",
                nobj,4,vTruth->at(0).Bjets_eta().at(0),multilepton_Bjet1_phi,vTruth->at(0).Bjets_pt().at(0),BJet1.M(),1,2,0,0,0));
    nobj++;


    // bj2 
    std::string b2_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d",
                nobj,4,vTruth->at(0).Bjets_eta().at(1),multilepton_Bjet2_phi,vTruth->at(0).Bjets_pt().at(1),BJet2.M(),1,2,0,0,0));
    nobj++;

    // met
    std::string met_fline = std::string(Form("%d     %d	  %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d	  %d",
                nobj,6,0.,Phi_0_2Pi(vTruth->at(0).metGen_phi()),vTruth->at(0).metGen_pt(),0.,0,0,0,0,0));
    nobj++;


    //    
    fout_MC << fline00   << std::endl;
    fout_MC << fline0    << std::endl;
    fout_MC << l1_fline  << std::endl;
    fout_MC << l2_fline  << std::endl;
    fout_MC << l3_fline  << std::endl;	  
    if (!(_processLHCO_MC == 7 || _processLHCO_MC == 8)) fout_MC << j1_fline  << std::endl;// don't print jets for ttW hypothesis
    if (!(_processLHCO_MC == 7 || _processLHCO_MC == 8)) fout_MC << j2_fline  << std::endl;// don't print jets for ttW hypothesis
    fout_MC << b1_fline  << std::endl;
    fout_MC << b2_fline  << std::endl;
    fout_MC << met_fline << std::endl;      

}

void TTbarHiggsMultileptonAnalysis::PrintLHCOforMadweight_RECO(int evt)
{

    if (vSelectedLeptons.size()!=3 || vSelectedBTagJets.size()<2 || vSelectedJets.size()<4 ) return; 

    fline0 = "0      " + std::string(Form("%d",evt)) + "     " + trig;

    int nobj = 1;	      

    int multilepton_Lepton1_Id_LHCO = -666;
    int multilepton_Lepton2_Id_LHCO = -666;
    int multilepton_Lepton3_Id_LHCO = -666;

    //
    // LHCO lepton ID convention
    //  
    if(abs(multilepton_Lepton1_Id)==11) multilepton_Lepton1_Id_LHCO = 1 ;
    else if(abs(multilepton_Lepton1_Id)==13) multilepton_Lepton1_Id_LHCO = 2 ;
    if(abs(multilepton_Lepton2_Id)==11) multilepton_Lepton2_Id_LHCO = 1 ;
    else if(abs(multilepton_Lepton2_Id)==13) multilepton_Lepton2_Id_LHCO = 2 ;
    if(abs(multilepton_Lepton3_Id)==11) multilepton_Lepton3_Id_LHCO = 1 ;
    else if(abs(multilepton_Lepton3_Id)==13) multilepton_Lepton3_Id_LHCO = 2 ;

    //
    // LHCO phi convention
    //

    float multilepton_Lepton1_phi       = Phi_0_2Pi(multilepton_Lepton1_P4.Phi());
    float multilepton_Lepton2_phi       = Phi_0_2Pi(multilepton_Lepton1_P4.Phi());
    float multilepton_Lepton3_phi       = Phi_0_2Pi(multilepton_Lepton1_P4.Phi());
    float multilepton_Bjet1_phi         = Phi_0_2Pi(multilepton_Bjet1_P4.Phi());
    float multilepton_Bjet2_phi         = Phi_0_2Pi(multilepton_Bjet2_P4.Phi());	 
    float multilepton_JetHighestPt1_phi = Phi_0_2Pi(multilepton_JetHighestPt1_P4.Phi());
    float multilepton_JetHighestPt2_phi = Phi_0_2Pi(multilepton_JetHighestPt1_P4.Phi());
    float multilepton_JetClosestMw1_phi = Phi_0_2Pi(multilepton_JetClosestMw1_P4.Phi());
    float multilepton_JetClosestMw2_phi = Phi_0_2Pi(multilepton_JetClosestMw2_P4.Phi());
    float multilepton_JetLowestMjj1_phi = Phi_0_2Pi(multilepton_JetLowestMjj1_P4.Phi());
    float multilepton_JetLowestMjj2_phi = Phi_0_2Pi(multilepton_JetLowestMjj2_P4.Phi());


    // l1
    std::string l1_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d", nobj,multilepton_Lepton1_Id_LHCO,multilepton_Lepton1_P4.Eta(),multilepton_Lepton1_phi,multilepton_Lepton1_P4.Pt(),0.0,multilepton_Lepton1_Id/abs(multilepton_Lepton1_Id),0,0,0,0));
    nobj++; 

    // l2
    std::string l2_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d", nobj,multilepton_Lepton2_Id_LHCO,multilepton_Lepton2_P4.Eta(),multilepton_Lepton2_phi,multilepton_Lepton2_P4.Pt(),0.0,multilepton_Lepton2_Id/abs(multilepton_Lepton2_Id),0,0,0,0));
    nobj++;	

    // l3
    std::string l3_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d", nobj,multilepton_Lepton3_Id_LHCO,multilepton_Lepton3_P4.Eta(),multilepton_Lepton3_phi,multilepton_Lepton3_P4.Pt(),0.0,multilepton_Lepton3_Id/abs(multilepton_Lepton3_Id),0,0,0,0));
    nobj++;

    //										    
    std::string j1_fline;
    std::string j2_fline;

    if ( _processLHCO_RECO == 5 || _processLHCO_RECO == 6 || _processLHCO_RECO == 4 || _processLHCO_RECO == 3 )
    {			     
        // j1
        j1_fline = std::string(Form("%d	  %d	 %.2f	    %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
                    nobj,4,multilepton_JetClosestMw1_P4.Eta(),multilepton_JetClosestMw1_phi,multilepton_JetClosestMw1_P4.Pt(),0.0,1,0,0,0,0));
        nobj++;

        // j2
        j2_fline = std::string(Form("%d	  %d	  %.2f      %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
                    nobj,4,multilepton_JetClosestMw2_P4.Eta(),multilepton_JetClosestMw2_phi,multilepton_JetClosestMw2_P4.Pt(),0.0,1,0,0,0,0));
        nobj++;
    }
    else if ( _processLHCO_RECO == 1 || _processLHCO_RECO == 2 ) 
    {
        // j1
        j1_fline = std::string(Form("%d	  %d	  %.2f      %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
                    nobj,4,multilepton_JetLowestMjj1_P4.Eta(),multilepton_JetLowestMjj1_phi,multilepton_JetLowestMjj1_P4.Pt(),0.0,1,0,0,0,0));
        nobj++;

        // j2
        j2_fline = std::string(Form("%d	  %d	   %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
                    nobj,4,multilepton_JetLowestMjj2_P4.Eta(),multilepton_JetLowestMjj2_phi,multilepton_JetLowestMjj2_P4.Pt(),0.0,1,0,0,0,0));
        nobj++;
    }

    // for B-jet mass
    TLorentzVector BJet1;
    BJet1.SetPtEtaPhiE(multilepton_Bjet1_P4.Pt(), multilepton_Bjet1_P4.Eta(), multilepton_Bjet1_P4.Phi(), multilepton_Bjet1_P4.E());

    TLorentzVector BJet2;
    BJet2.SetPtEtaPhiE(multilepton_Bjet2_P4.Pt(), multilepton_Bjet2_P4.Eta(), multilepton_Bjet2_P4.Phi(), multilepton_Bjet2_P4.E());


    // bj1
    std::string b1_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d",
                nobj,4,multilepton_Bjet1_P4.Eta(),multilepton_Bjet1_phi,multilepton_Bjet1_P4.Pt(),BJet1.M(),1,2,0,0,0));
    nobj++;


    // bj2 
    std::string b2_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d",
                nobj,4,multilepton_Bjet2_P4.Eta(),multilepton_Bjet2_phi,multilepton_Bjet2_P4.Pt(),BJet2.M(),1,2,0,0,0));
    nobj++;

    // met
    std::string met_fline = std::string(Form("%d     %d	  %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d	  %d",
                nobj,6,0.,Phi_0_2Pi(multilepton_mET.Phi()),multilepton_mET.Pt(),0.,0,0,0,0,0));
    nobj++;


    //    
    fout_RECO << fline00	<< std::endl;
    fout_RECO << fline0	<< std::endl;
    fout_RECO << l1_fline  << std::endl;
    fout_RECO << l2_fline  << std::endl;
    fout_RECO << l3_fline  << std::endl;	
    if (!(_processLHCO_RECO == 7 || _processLHCO_RECO == 8)) fout_RECO << j1_fline  << std::endl;// don't print jets for ttW hypothesis
    if (!(_processLHCO_RECO == 7 || _processLHCO_RECO == 8)) fout_RECO << j2_fline  << std::endl;// don't print jets for ttW hypothesis
    fout_RECO << b1_fline  << std::endl;
    fout_RECO << b2_fline  << std::endl;
    fout_RECO << met_fline << std::endl;      

}

/*void TTbarHiggsMultileptonAnalysis::ProducePUweight()
  {

// TTbarHiggsMultileptonAnalysis::TTbarHiggsMultileptonAnalysis(
// TString inputFileName, TChain *tree, TString the_sampleName, TString treeName, TString outputFile, bool isdata, float xsec, float lumi, int nowe, int nmax)

std::cout << "Initializing PU from MC..." << std::endl;

const char *fname_str1  = "input_WZJets_MC.txt";
const char *stream_str1 = "Nt";
TString inputFileName1 = *fname_str1;
TChain *tree1 = 0;
TString treeName1 = *stream_str1;
int nmax1 = -1;

std::ifstream infile1;
infile1.open(inputFileName1);
std::string ifile1 = "";
while( getline(infile1, ifile1) )
{
std::string fnameStr1 = std::string(ifile1);
tree1->Add(fnameStr1.c_str());
std::cout << "file: " << fnameStr1 << std::endl;
}
infile1.close();

Init(tree1);

TH1F * PU_MC = new TH1F("PU_MC","PU_MC",100,0,99);
//TH1F PU_MC("PU_MC","PU_MC",100,0,99); 

Long64_t nentries1 = fChain->GetEntries();
int nentries_max1 = nentries1;

for (Long64_t jentry=0; jentry<nentries_max1;jentry++)
{
std::cout << "n[" << jentry << "] / " << nentries_max1 << std::endl;
int PU_M = 0;
PU_M = vEvent->at(0).pv_n();
std::cout << "Number of PV in MC: " << PU_M << std::endl;
PU_MC->Fill(PU_M,1); 
}

std::cout << "Initializing PU from Data..." << std::endl;

const char *fname_str2  = "input_DoubleEG_DATA.txt";
const char *stream_str2 = "Nt";
TString inputFileName2 = *fname_str2;
TChain *tree2 = 0;
TString treeName2 = *stream_str2;
int nmax2 = -1;

std::ifstream infile2;
infile2.open(inputFileName2);
std::string ifile2 = "";
while( getline(infile2, ifile2) )
{
std::string fnameStr2 = std::string(ifile2);
tree2->Add(fnameStr2.c_str());
std::cout << "file: " << fnameStr2 << std::endl;
}
infile2.close();

Init(tree2);

TH1F * PU_DATA = new TH1F("PU_DATA","PU_DATA",100,0,99);
//TH1F PU_DATA("PU_DATA","PU_DATA",100,0,99); 

Long64_t nentries2 = fChain->GetEntries();
int nentries_max2 = nentries2;

for (Long64_t jentry=0; jentry<nentries_max2;jentry++)
{
    if (fChain == 0) break;
    std::cout << "n[" << jentry << "] / " << nentries_max2 << std::endl;
    int PU_D = 0;
    //PU_D = vEvent->at(0).pv_n();
    std::cout << "Number of PV in DATA: " << PU_D << std::endl;
    PU_DATA->Fill(PU_D,1);
}

TH1F * PU_weight = new TH1F("PU_weight","PU_weight",100,0,99);
PU_DATA->Divide(PU_MC);
PU_weight = PU_DATA;

std::cout << "PU weight produced..." << std::endl;
}


float TTbarHiggsMultileptonAnalysis::PUweight()
{
    float weight = 1;
    return weight;
}*/

void TTbarHiggsMultileptonAnalysis::selectBjets(std::string BjetSel, int* ibsel1, int* ibsel2, bool doSelectOnlyBjets){

    //Selects the two highest b-tag jets. If only one b-tag select just this one.
    int ib1=-1, ib2=-1;

    if (BjetSel=="HighestBtagDiscrim"){
        float btag_max=-999, btag_max2=-999;
        for (unsigned int ib=0; ib<vSelectedJets.size(); ib++){
            if (doSelectOnlyBjets && (vSelectedJets.at(ib).CSVv2()<0.423)) continue;
            if (vSelectedJets.at(ib).CSVv2()>btag_max){
                btag_max2 = btag_max;
                ib2 = ib1;
                btag_max = vSelectedJets.at(ib).CSVv2();
                ib1 = ib;
            }
            if (vSelectedJets.at(ib).CSVv2()<btag_max && vSelectedJets.at(ib).CSVv2()>btag_max2){
                btag_max2 = vSelectedJets.at(ib).CSVv2();
                ib2 = ib;
            }
        }
    }
    if (BjetSel=="BtagHighestPt"){
        float pt_max=0, pt_max2=0;
        for (unsigned int ib=0; ib<vSelectedJets.size(); ib++){
//            if (vSelectedJets.at(ib).CSVv2()<0.423) continue;
            if (vSelectedJets.at(ib).pt()>pt_max){
                pt_max2 = pt_max;
                ib2 = ib1;
                pt_max = vSelectedJets.at(ib).pt();
                ib1 = ib;
            }
            if (vSelectedJets.at(ib).pt()<pt_max && vSelectedJets.at(ib).pt()>pt_max2){
                pt_max2 = vSelectedJets.at(ib).pt();
                ib2 = ib;
            }
        }
    }

    *ibsel1 = ib1;
    *ibsel2 = ib2;

}

float TTbarHiggsMultileptonAnalysis::Phi_0_2Pi(float phi)
{
    float phi_0_2pi = phi;
    if (phi>= TMath::Pi()) phi_0_2pi  -= 2.*TMath::Pi();
    if (phi<0.)            phi_0_2pi  += 2.*TMath::Pi();
    return phi_0_2pi;
}

float TTbarHiggsMultileptonAnalysis::GetDeltaR(float eta1,float phi1,float eta2,float phi2)
{
    float DeltaPhi = TMath::Abs(phi2 - phi1);
    if (DeltaPhi > 3.141593 ) DeltaPhi -= 2.*3.141593;
    return TMath::Sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}
