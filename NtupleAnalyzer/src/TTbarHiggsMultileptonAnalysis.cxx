#include "../include/TTbarHiggsMultileptonAnalysis.h"
#include "TSystem.h"
#include "Helper.cxx"

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

	
	if ( !_isdata )
	{
          weight = _lumi*_xsec/_nowe;
          mc_weight = vEvent->at(0).mc_weight();
          //weight_PV = _h_PV->GetBinContent(pvn);
          weight = weight * mc_weight; //*weight_PV;
	}
	else 
	{
	    weight    = 1.;
        mc_weight = 1.;
        weight_PV = 1.; 

        int tab[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        int a = vEvent->at(0).ev_trigger_pass_byname_1();
        int n = 0, size = 0;

        do{
            n = a%10;
            a = a/10;
            tab[size] = n;
            size = size + 1;
        }while(a!=0);

        bool E = false, M = false, EE = false, MM = false, EM = false;
        int result_trigger = 0;
        if (size > 2) {if ( tab[3] == 1                ) E  = true;}
        if (size > 1) {if ( tab[2] == 1 || tab[2] == 2 ) M  = true;}
        if (size > 0) {if ( tab[1] == 1                ) EE = true;}
        if (size > 1) {if ( tab[1] == 2 || tab[1] == 5 ) MM = true;}
        if ( tab[0] == 2 || tab[0] == 5 )                EM = true;

        bool emdataset = _sampleName.Contains("MuonEG");
        bool mmdataset = _sampleName.Contains("DoubleMuon");
        bool eedataset = _sampleName.Contains("DoubleEG");
        bool mdataset  = _sampleName.Contains("SingleMuon");
        bool edataset  = _sampleName.Contains("SingleElectron");

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


        //---------------------------
        // initialisation
        //---------------------------
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
        is_2lss_TT_CR     = false;

        is_3l_TTH_SR      = false;
        is_3l_WZ_CR       = false; 
        is_3l_WZrel_CR    = false;
        is_3l_TTZ_CR      = false;
        is_CR_Zl          = false;
        
	//---------------------------
        //trigger
        //---------------------------
        is_trigger = false;
        if ( vEvent->at(0).ev_trigger_pass_byname_1() >= 1 ) is_trigger = true;

        //---------------------------
        //muons
        //---------------------------
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

        //---------------------------
        // electrons
        //---------------------------
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

        //---------------------------
        // taus
        //---------------------------
        for(unsigned int itau=0; itau < vTau->size() ; itau++)
        {
            Lepton l; l.setLepton(&vTau->at(itau),itau,1);

            vSelectedTaus.push_back(vTau->at(itau));
            //vSelectedLeptons.push_back(l);

            //vLeptons.push_back(l);

            theHistoManager->fillHisto("TauPt",                             "noSel",        "",   "",  vTau->at(itau).pt(),               weight);
            theHistoManager->fillHisto("TauEta",                            "noSel",        "",   "",  vTau->at(itau).eta(),              weight);
        }

        std::sort(vSelectedLeptons.begin(), vSelectedLeptons.end(), SortingLeptonPt);

        //---------------------------
        // b-jets
        //---------------------------
        nLooseBJets  = 0;
        nMediumBJets = 0;

        for(unsigned int ijet=0; ijet < vJet->size() ; ijet++)
        {

            if( vJet->at(ijet).CSVv2() > 0.423 ) nLooseBJets++;
            if( vJet->at(ijet).CSVv2() > 0.814 ) nMediumBJets++;

            if(vJet->at(ijet).CSVv2() >= 0.423 ) vSelectedBTagJets.push_back(vJet->at(ijet));
            else                                 vSelectedNonBTagJets.push_back(vJet->at(ijet));
            if(vJet->at(ijet).CSVv2() >= 0.814 ) vSelectedMediumBTagJets.push_back(vJet->at(ijet));

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

        //---------------------------
        //Selection for signal and control regions
        //---------------------------
       
        TwoLeptonsSameSignSelection_TTH2l(jentry);
        //TwoLeptonsSameSignSelection_LepMVA_sideband(jentry);
        //TwoLeptonsSameSignSelection_JetMultiplicity_sideband(jentry);
        TwoLeptonsSameSignSelection_TT_CR(jentry);

        ThreeLeptonSelection_TTH3l(jentry);
        ThreeLeptonSelection_CR_WZ(jentry);
        ThreeLeptonSelection_CR_WZrelaxed(jentry);
        ThreeLeptonSelection_CR_Zl(jentry);
        ThreeLeptonSelection_TTZ(jentry);

        //std::cout <<is_CR_TTl<<" "<< is_CR_Zl <<" " << is_CR_WZ<<" " << is_TTH3l<< std::endl;
        //if (is_TTH3l==true ) std::cout <<"is_TTH3l" << std::endl;
        if ( is_3l_TTH_SR || is_3l_TTZ_CR ) fillOutputTree();
	//if (is_2lss_TTH_SR) fillOutputTree();

        //---------------------------
        //Madweight LHCO stuff
        //---------------------------
        if ( !_isdata && _printLHCO_MC && ThreeLeptonSelection_TTH3l_MC()) PrintLHCOforMadweight_MC(jentry);
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

    bool same_sign      = ( vSelectedLeptons.at(0).charge() == vSelectedLeptons.at(1).charge() );
    if(!same_sign)      return;

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

    is_2lss_TTH_SR = true;   

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

void TTbarHiggsMultileptonAnalysis::TwoLeptonsSameSignSelection_TT_CR(int evt)
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

    is_2lss_TT_CR = true;   

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

    nLep                  = ( vSelectedLeptons.size()				    >= 3 );
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

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        for(int j=0; j<vSelectedLeptons.size(); j++)
        {
            if (  ( i                           != j                                                          )
               && ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id()                               )
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

    for(int i=0; i<vLeptons.size(); i++)
    {
        for(int j=0; j<vLeptons.size(); j++)
        {
            if (  ( i                           != j                                                          )
               && ( vLeptons.at(i).id() == -vLeptons.at(j).id()                               )
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

    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel", "TTZ_CR",   "",    vLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel", "TTZ_CR",   "",    vLeptons.size(),   vSelectedBTagJets.size(),  weight);

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

    theHistoManager->fillHisto("CutFlow",                           "noSel", "TTZ_CR",   "",   2, weight);
    theHistoManager->fillHisto("CutFlow",                 "PassingTightMVA", "TTZ_CR",   "",   1, weight);

    // ##########
    // # Z peak #
    // ##########

    bool pass_Zpeak = false;
    float Delta = 10, ZM = -1., Zpt = -1., Lep1Z = -1, Lep2Z = -1;

    for(int i=0; i<vLeptons.size(); i++)
    {
        for(int j=0; j<vLeptons.size(); j++)
        {
            if (  ( i                           != j                                                          )
               && ( vLeptons.at(i).id() == -vLeptons.at(j).id()                               )
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

    for(int i=0; i<vLeptons.size(); i++)
    {
        lepton_px = lepton_px + vLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vLeptons.at(i).p4().Py();
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
    for(int i=0; i<vLeptons.size(); i++)
    {
        for(int j=0; j<vLeptons.size(); j++)
        {
            if (  ( i                   != j                            )
               && ( vLeptons.at(i).id() == -vLeptons.at(j).id() ) )
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

    theHistoManager->fillHisto("LeadingLeptonPt",                "FinalCut", "TTZ_CR",   "",   vLeptons.at(0).pt()   , weight);
    theHistoManager->fillHisto("SubLeadingLeptonPt",             "FinalCut", "TTZ_CR",   "",   vLeptons.at(1).pt()   , weight);    

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
        theHistoManager->fillHisto("LeadingLeptonPt",                "FinalCut4j", "TTZ_CR",   "",   vLeptons.at(0).pt()  , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",             "FinalCut4j", "TTZ_CR",   "",   vLeptons.at(1).pt()  , weight);

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
        is_CR_Zl = true;
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
    tOutput->Branch("mc_weight",&mc_weight,"mc_weight/F"); 
    tOutput->Branch("weight",&weight,"weight/F");
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

    tOutput->Branch("is_2lss_TT_CR",&is_2lss_TT_CR,"is_2lss_TT_CR/B");
    //tOutput->Branch("is_3l_WZ_CR",&is_3l_WZ_CR,"is_3l_WZ_CR/B");
    tOutput->Branch("is_3l_TTZ_CR",&is_3l_TTZ_CR,"is_3l_TTZ_CR/B");
    tOutput->Branch("is_CR_Zl",&is_CR_Zl,"is_CR_Zl/B");

    tOutput->Branch("is_trigger",&is_trigger,"is_trigger/B");

    tOutput->Branch("multilepton_Bjet1_Id",&multilepton_Bjet1_Id,"multilepton_Bjet1_Id/I");
    tOutput->Branch("multilepton_Bjet1_P4","TLorentzVector",&multilepton_Bjet1_P4);
    tOutput->Branch("multilepton_Bjet2_Id",&multilepton_Bjet2_Id,"multilepton_Bjet2_Id/I");
    tOutput->Branch("multilepton_Bjet2_P4","TLorentzVector",&multilepton_Bjet2_P4);
    tOutput->Branch("multilepton_Lepton1_Id",&multilepton_Lepton1_Id,"multilepton_Lepton1_Id/I");
    tOutput->Branch("multilepton_Lepton1_P4","TLorentzVector",&multilepton_Lepton1_P4);
    tOutput->Branch("multilepton_Lepton2_Id",&multilepton_Lepton2_Id,"multilepton_Lepton2_Id/I");
    tOutput->Branch("multilepton_Lepton2_P4","TLorentzVector",&multilepton_Lepton2_P4);
    tOutput->Branch("multilepton_Lepton3_Id",&multilepton_Lepton3_Id,"multilepton_Lepton3_Id/I");
    tOutput->Branch("multilepton_Lepton3_P4","TLorentzVector",&multilepton_Lepton3_P4);
    tOutput->Branch("multilepton_Lepton4_Id",&multilepton_Lepton4_Id,"multilepton_Lepton4_Id/I");
    tOutput->Branch("multilepton_Lepton4_P4","TLorentzVector",&multilepton_Lepton4_P4);
    tOutput->Branch("multilepton_JetHighestPt1_Id",&multilepton_JetHighestPt1_Id,"multilepton_JetHighestPt1_Id/I");
    tOutput->Branch("multilepton_JetHighestPt1_P4","TLorentzVector",&multilepton_JetHighestPt1_P4);
    tOutput->Branch("multilepton_JetHighestPt2_Id",&multilepton_JetHighestPt2_Id,"multilepton_JetHighestPt2_Id/I");
    tOutput->Branch("multilepton_JetHighestPt2_P4","TLorentzVector",&multilepton_JetHighestPt2_P4);
    tOutput->Branch("multilepton_JetClosestMw1_Id",&multilepton_JetClosestMw1_Id,"multilepton_JetClosestMw1_Id/I");
    tOutput->Branch("multilepton_JetClosestMw1_P4","TLorentzVector",&multilepton_JetClosestMw1_P4);
    tOutput->Branch("multilepton_JetClosestMw2_Id",&multilepton_JetClosestMw2_Id,"multilepton_JetClosestMw2_Id/I");
    tOutput->Branch("multilepton_JetClosestMw2_P4","TLorentzVector",&multilepton_JetClosestMw2_P4);
    tOutput->Branch("multilepton_JetLowestMjj1_Id",&multilepton_JetLowestMjj1_Id,"multilepton_JetLowestMjj1_Id/I");
    tOutput->Branch("multilepton_JetLowestMjj1_P4","TLorentzVector",&multilepton_JetLowestMjj1_P4);
    tOutput->Branch("multilepton_JetLowestMjj2_Id",&multilepton_JetLowestMjj2_Id,"multilepton_JetLowestMjj2_Id/I");
    tOutput->Branch("multilepton_JetLowestMjj2_P4","TLorentzVector",&multilepton_JetLowestMjj2_P4);

    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_Id",&multilepton_JetHighestPt1_2ndPair_Id,"multilepton_JetHighestPt1_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_P4","TLorentzVector",&multilepton_JetHighestPt1_2ndPair_P4);
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_Id",&multilepton_JetHighestPt2_2ndPair_Id,"multilepton_JetHighestPt2_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_P4","TLorentzVector",&multilepton_JetHighestPt2_2ndPair_P4);
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_Id",&multilepton_JetClosestMw1_2ndPair_Id,"multilepton_JetClosestMw1_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_P4","TLorentzVector",&multilepton_JetClosestMw1_2ndPair_P4);
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_Id",&multilepton_JetClosestMw2_2ndPair_Id,"multilepton_JetClosestMw2_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_P4","TLorentzVector",&multilepton_JetClosestMw2_2ndPair_P4);
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_Id",&multilepton_JetLowestMjj1_2ndPair_Id,"multilepton_JetLowestMjj1_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_P4","TLorentzVector",&multilepton_JetLowestMjj1_2ndPair_P4);
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_Id",&multilepton_JetLowestMjj2_2ndPair_Id,"multilepton_JetLowestMjj2_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_P4","TLorentzVector",&multilepton_JetLowestMjj2_2ndPair_P4);

    tOutput->Branch("multilepton_mET","TLorentzVector",&multilepton_mET);
    tOutput->Branch("multilepton_Ptot","TLorentzVector",&multilepton_Ptot);

    return;
}


void TTbarHiggsMultileptonAnalysis::fillOutputTree(){

    bool is2lss=false, is3l=false, is4l=false;

    if ( vSelectedLeptons.size()>=4 ) is4l = true;
    if ( vSelectedLeptons.size()==3 || (vSelectedLeptons.size() == 2 &&  vFakeLeptons.size() == 1 )) is3l = true;
    if ( vSelectedLeptons.size()==2 ) is2lss = true;
    if (!is2lss && !is3l && !is4l) return;

    if (vSelectedJets.size()<2) return;
    if (!(vSelectedBTagJets.size()>=2 || (vSelectedMediumBTagJets.size()==1))) return; 
    
    //if (vSelectedLeptons.size()<2) return; // 2lss only at the moment
    if (vSelectedLeptons.size()<4) return; // 4l only at the moment

    //std::cout << "lept="<<vSelectedLeptons.size()<<" fake="<<vFakeLeptons.size()<<std::endl;
    //std::cout << "btag="<<vSelectedBTagJets.size()<<" nonbtag="<<vSelectedNonBTagJets.size()<<std::endl;

    //2lss
    if (is2lss && vSelectedBTagJets.size()>=2 && vSelectedJets.size()-2>=4) catJets = kCat_2lss_2b_4j;
    else if (is2lss && vSelectedBTagJets.size()==1 && vSelectedJets.size()-1>=4) catJets = kCat_2lss_1b_4j;
    else if (is2lss && vSelectedBTagJets.size()>=2 && vSelectedJets.size()-2==3) catJets = kCat_2lss_2b_3j;
    else if (is2lss && vSelectedBTagJets.size()==1 && vSelectedJets.size()-1==3) catJets = kCat_2lss_1b_3j;
    else if (is2lss && vSelectedBTagJets.size()>=2 && vSelectedJets.size()-2==2) catJets = kCat_2lss_2b_2j;
    //4l 
    else if (is4l && vSelectedBTagJets.size()>=2) catJets = kCat_4l_2b;
    else if (is4l && vSelectedBTagJets.size()==1) catJets = kCat_4l_1b;
    //3l
    else if (is3l && vSelectedBTagJets.size()>=2 && vSelectedJets.size()-2>=2) catJets = kCat_3l_2b_2j;
    else if (is3l && vSelectedBTagJets.size()==1 && vSelectedJets.size()-1>=2) catJets = kCat_3l_1b_2j;
    else if (is3l && vSelectedBTagJets.size()>=2 && vSelectedJets.size()-2==1) catJets = kCat_3l_2b_1j;
    else if (is3l && vSelectedBTagJets.size()==1 && vSelectedJets.size()-1==1) catJets = kCat_3l_1b_1j;
    else if (is3l && vSelectedBTagJets.size()>=2 && vSelectedJets.size()-2==0) catJets = kCat_3l_2b_0j;
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
    else if (vSelectedLeptons.size()==2 && vFakeLeptons.size()==1)
    {
       multilepton_Lepton3_P4 = vFakeLeptons.at(0).p4();
       multilepton_Lepton3_Id = vFakeLeptons.at(0).id();
    }

    if (vSelectedLeptons.size()>=4)
    {
        multilepton_Lepton4_P4 = vSelectedLeptons.at(3).p4();
        multilepton_Lepton4_Id = vSelectedLeptons.at(3).id();
    }


            //Choosing 2 b-jets
            TLorentzVector Bjet1, Bjet2; 
            int ib1=-1, ib2=-1;
            selectBjets("HighestBtagDiscrim", &ib1, &ib2);
            if (ib1!=-1) Bjet1.SetPtEtaPhiE(vSelectedJets.at(ib1).pt(), vSelectedJets.at(ib1).eta(), vSelectedJets.at(ib1).phi(), vSelectedJets.at(ib1).E());
            if (ib2!=-1) Bjet2.SetPtEtaPhiE(vSelectedJets.at(ib2).pt(), vSelectedJets.at(ib2).eta(), vSelectedJets.at(ib2).phi(), vSelectedJets.at(ib2).E());

	    multilepton_Bjet1_Id = -999;
	    if (ib1!=-1){
              multilepton_Bjet1_P4 = Bjet1;
              multilepton_Bjet1_Id = 5;
	    }
	    multilepton_Bjet2_Id = -999;
	    if (ib2!=-1){
              multilepton_Bjet2_P4 = Bjet2;
              multilepton_Bjet2_Id = 5;
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

            if (ij1!=-1 && ij2==-1){
                multilepton_JetHighestPt1_Id = 1;
                multilepton_JetHighestPt1_P4.SetPtEtaPhiE(vSelectedJets.at(ij1).pt(), vSelectedJets.at(ij1).eta(), vSelectedJets.at(ij1).phi(), vSelectedJets.at(ij1).E());
            }   
            if (ij1!=-1 && ij2!=-1) {
                multilepton_JetHighestPt1_Id = 1;
                multilepton_JetHighestPt2_Id = 1;
                multilepton_JetHighestPt1_P4.SetPtEtaPhiE(vSelectedJets.at(ij1).pt(), vSelectedJets.at(ij1).eta(), vSelectedJets.at(ij1).phi(), vSelectedJets.at(ij1).E());
                multilepton_JetHighestPt2_P4.SetPtEtaPhiE(vSelectedJets.at(ij2).pt(), vSelectedJets.at(ij2).eta(), vSelectedJets.at(ij2).phi(), vSelectedJets.at(ij2).E());
            }
            if (ik1!=-1 && ik2!=-1){
                multilepton_JetClosestMw1_Id = 2;
                multilepton_JetClosestMw2_Id = 2;
                multilepton_JetClosestMw1_P4.SetPtEtaPhiE(vSelectedJets.at(ik1).pt(), vSelectedJets.at(ik1).eta(), vSelectedJets.at(ik1).phi(), vSelectedJets.at(ik1).E());
                multilepton_JetClosestMw2_P4.SetPtEtaPhiE(vSelectedJets.at(ik2).pt(), vSelectedJets.at(ik2).eta(), vSelectedJets.at(ik2).phi(), vSelectedJets.at(ik2).E());
            }
            if (il1!=-1 && il2!=-1){
                multilepton_JetLowestMjj1_Id = 3;
                multilepton_JetLowestMjj2_Id = 3;
                multilepton_JetLowestMjj1_P4.SetPtEtaPhiE(vSelectedJets.at(il1).pt(), vSelectedJets.at(il1).eta(), vSelectedJets.at(il1).phi(), vSelectedJets.at(il1).E());
                multilepton_JetLowestMjj2_P4.SetPtEtaPhiE(vSelectedJets.at(il2).pt(), vSelectedJets.at(il2).eta(), vSelectedJets.at(il2).phi(), vSelectedJets.at(il2).E());
            }

            //2nd pair (first one: closest to Mw)
            if (is2lss && ij1!=-1 && ij2!=-1){
              if (im1!=-1 && im2==-1){
                multilepton_JetHighestPt1_2ndPair_Id = 1;
                multilepton_JetHighestPt1_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(im1).pt(), vSelectedJets.at(im1).eta(), vSelectedJets.at(im1).phi(), vSelectedJets.at(im1).E());
              }
              if (im1!=-1 && im2!=-1){
                multilepton_JetHighestPt1_2ndPair_Id = 1;
                multilepton_JetHighestPt1_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(im1).pt(), vSelectedJets.at(im1).eta(), vSelectedJets.at(im1).phi(), vSelectedJets.at(im1).E());
                multilepton_JetHighestPt2_2ndPair_Id = 1;
                multilepton_JetHighestPt2_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(im2).pt(), vSelectedJets.at(im2).eta(), vSelectedJets.at(im2).phi(), vSelectedJets.at(im2).E());   
              }
              if (io1!=-1 && io2!=-1){
                multilepton_JetClosestMw1_2ndPair_Id = 2;
                multilepton_JetClosestMw2_2ndPair_Id = 2;
                multilepton_JetClosestMw1_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(io1).pt(), vSelectedJets.at(io1).eta(), vSelectedJets.at(io1).phi(), vSelectedJets.at(io1).E());
                multilepton_JetClosestMw2_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(io2).pt(), vSelectedJets.at(io2).eta(), vSelectedJets.at(io2).phi(), vSelectedJets.at(io2).E());
              }
              if (ip1!=-1 && ip2!=-1){
                multilepton_JetLowestMjj1_2ndPair_Id = 3;
                multilepton_JetLowestMjj2_2ndPair_Id = 3;
                multilepton_JetLowestMjj1_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(ip1).pt(), vSelectedJets.at(ip1).eta(), vSelectedJets.at(ip1).phi(), vSelectedJets.at(ip1).E());
                multilepton_JetLowestMjj2_2ndPair_P4.SetPtEtaPhiE(vSelectedJets.at(ip2).pt(), vSelectedJets.at(ip2).eta(), vSelectedJets.at(ip2).phi(), vSelectedJets.at(ip2).E());
              }
            }

            multilepton_mET.SetPtEtaPhiE(vEvent->at(0).metpt(), 0, vEvent->at(0).metphi(), vEvent->at(0).metpt());

            mc_ttZhypAllowed = 0;

            if(vSelectedLeptons.size()==3) 
            {
                if ( vSelectedLeptons.at(0).charge()==vSelectedLeptons.at(1).charge() && vSelectedLeptons.at(1).charge()==vSelectedLeptons.at(2).charge() ) mc_ttZhypAllowed =-1;
                else if (  ( vSelectedLeptons.at(0).id() == -vSelectedLeptons.at(1).id() ) 
                        || ( vSelectedLeptons.at(0).id() == -vSelectedLeptons.at(2).id() ) 
                        || ( vSelectedLeptons.at(1).id() == -vSelectedLeptons.at(2).id() ))
                    mc_ttZhypAllowed = 1; }


                mc_nJets25 = vSelectedJets.size();
                mc_nBtagJets25 = vSelectedBTagJets.size();
                mc_nMediumBtagJets25 = vSelectedMediumBTagJets.size();
                mc_nNonBtagJets25 = vSelectedNonBTagJets.size();

                tOutput->Fill();

                //if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);

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

void TTbarHiggsMultileptonAnalysis::selectBjets(std::string BjetSel, int* ibsel1, int* ibsel2){

    //Selects the two highest b-tag jets. If only one b-tag select just this one.
    int ib1=-1, ib2=-1;

    if (BjetSel=="HighestBtagDiscrim"){
        float btag_max=-999, btag_max2=-999;
        for (unsigned int ib=0; ib<vSelectedJets.size(); ib++){
            if (vSelectedJets.at(ib).CSVv2()<0.423) continue;
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
            if (vSelectedJets.at(ib).CSVv2()<0.423) continue;
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
