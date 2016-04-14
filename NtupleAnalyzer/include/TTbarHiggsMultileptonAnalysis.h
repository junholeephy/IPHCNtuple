#ifndef TTbarHiggsMultileptonAnalysis_H
#define TTbarHiggsMultileptonAnalysis_H

#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TChain.h>
#include <TFile.h>
#include "Event.h"
#include "Muon.h"
#include "Tau.h"
#include "Electron.h"
#include "Jet.h"
#include "Truth.h"
#include <TObject.h>
#include <TROOT.h>
#include <iostream>
#include <fstream>

#include "HistoManager.h"
#include "Lepton.h"

class TTbarHiggsMultileptonAnalysis
{

    public:

        TTbarHiggsMultileptonAnalysis();
        TTbarHiggsMultileptonAnalysis(TString inputFileName, TChain *tree, TString sampleName, TString treeName,
                TString outputFileName, bool isdata, float xsec, float lumi, int nowe, int nmax);

        ~TTbarHiggsMultileptonAnalysis();

        void createHistograms();
        void writeHistograms();

        void TwoLeptonsSameSignSelection_TTH2l(int evt);

        void TwoLeptonsSameSignSelection_LepMVA_sideband(int evt);
        void TwoLeptonsSameSignSelection_JetMultiplicity_sideband(int evt);
        void DiLeptonSelection_TT_CR(int evt);

        void ThreeLeptonSelection_TTH3l(int evt);

        void ThreeLeptonSelection_CR_WZ(int evt);
        void ThreeLeptonSelection_CR_WZrelaxed(int evt);
        void ThreeLeptonSelection_TTZ(int evt);

        void ThreeLeptonSelection_CR_Zl(int evt);

        bool ThreeLeptonSelection_TTH3l_MC();


        TChain *fChain;   //!pointer to the analyzed TTree or TChain

        std::vector<Event>    *vEvent     = new std::vector<Event>();
        std::vector<Electron> *vElectron  = new std::vector<Electron>();
        std::vector<Muon>     *vMuon      = new std::vector<Muon>();
        std::vector<Tau>      *vTau       = new std::vector<Tau>();
        std::vector<Jet>      *vJet       = new std::vector<Jet>();
        std::vector<Truth>    *vTruth     = new std::vector<Truth>();

        std::vector<Lepton>   vLeptons;
        std::vector<Muon>	  vSelectedMuons;
        std::vector<Electron> vSelectedElectrons;
        std::vector<Tau>      vSelectedTaus;
        std::vector<Lepton>   vSelectedLeptons;
        std::vector<Muon>	  vFakeMuons;     // inverted MVA cut
        std::vector<Electron> vFakeElectrons; // inverted MVA cut
        std::vector<Lepton>   vFakeLeptons;   // inverted MVA cut
        std::vector<Jet>	  vSelectedNonBTagJets;
        std::vector<Jet>	  vSelectedBTagJets;
        std::vector<Jet>      vSelectedMediumBTagJets;
        std::vector<Jet>	  vSelectedJets;

        int nLooseBJets;
        int nMediumBJets;

        bool is_2lss_TTH_SR;
        bool is_2lss_JM_SB;
        bool is_2lss_LepMVA_SB;
        bool is_emu_TT_CR;

        bool is_3l_TTH_SR;    // TTH 3l analysis
        bool is_3l_WZ_CR;     // WZ CR w/ 3l (selected) or more, no b-jets, Z peak
        bool is_3l_WZrel_CR;  // WZ CR w/ 3l (loose) or more, no medium b-jets, Z peak
        bool is_3l_TTZ_CR;    // TTZ 3l analysis (for the future..)

        bool is_Zl_CR;

        bool is_trigger;

        virtual void     Init(TChain *tree);
        virtual void     Loop();

        void initializeOutputTree();
        void selectBjets(std::string, int*, int*);
        void fillOutputTree();

        // needed to print info in LHCO text format (madweight)
        void InitLHCO(int process_MC, int process_RECO);
        void PrintLHCOforMadweight_MC(int evt);
        void PrintLHCOforMadweight_RECO(int evt);

        float Phi_0_2Pi(float phi);
        float GetDeltaR(float eta1,float phi1,float eta2,float phi2);

        TTree* tOutput;
        Int_t mc_event;
        Float_t weight;
        Float_t mc_weight;
        Float_t weight_PV; // PU reweighting from PV distribution
        Int_t mc_3l_category, mc_ttbar_decay, mc_boson_decay, mc_ttZhypAllowed, mc_nJets25, mc_nBtagJets25, mc_nMediumBtagJets25, mc_nNonBtagJets25;
        Int_t catJets;
        Int_t multilepton_Bjet1_Id, multilepton_Bjet2_Id;
        Int_t multilepton_Lepton1_Id, multilepton_Lepton2_Id, multilepton_Lepton3_Id, multilepton_Lepton4_Id;
        Int_t multilepton_JetHighestPt1_Id, multilepton_JetHighestPt2_Id, multilepton_JetClosestMw1_Id, multilepton_JetClosestMw2_Id, multilepton_JetLowestMjj1_Id, multilepton_JetLowestMjj2_Id;
        Int_t multilepton_JetHighestPt1_2ndPair_Id, multilepton_JetHighestPt2_2ndPair_Id, multilepton_JetClosestMw1_2ndPair_Id, multilepton_JetClosestMw2_2ndPair_Id, multilepton_JetLowestMjj1_2ndPair_Id, multilepton_JetLowestMjj2_2ndPair_Id;
        TLorentzVector multilepton_Bjet1_P4, multilepton_Bjet2_P4;
        TLorentzVector multilepton_Lepton1_P4, multilepton_Lepton2_P4, multilepton_Lepton3_P4, multilepton_Lepton4_P4;
        TLorentzVector multilepton_JetHighestPt1_P4, multilepton_JetHighestPt2_P4, multilepton_JetClosestMw1_P4, multilepton_JetClosestMw2_P4, multilepton_JetLowestMjj1_P4, multilepton_JetLowestMjj2_P4;
        TLorentzVector  multilepton_JetHighestPt1_2ndPair_P4, multilepton_JetHighestPt2_2ndPair_P4, multilepton_JetClosestMw1_2ndPair_P4, multilepton_JetClosestMw2_2ndPair_P4, multilepton_JetLowestMjj1_2ndPair_P4, multilepton_JetLowestMjj2_2ndPair_P4;
        TLorentzVector multilepton_mET, multilepton_Ptot;
	Float_t multilepton_mHT;
	Double_t multilepton_mETcov00;
        Double_t multilepton_mETcov01;
        Double_t multilepton_mETcov10;
        Double_t multilepton_mETcov11;

    private:

        HistoManager * theHistoManager;

        TFile * outputfile;

        TString _sampleName;
        TString _outputFileName;

        bool _isdata;
        float _xsec;
        float _lumi;
        int   _nowe; // number of weighted events
        int   _nmax; // max number of events to process

        TFile * _file_PVreweighting;
        TH1F* _h_PV;

        // needed to print info in LHCO text format (madweight)

        bool _printLHCO_MC;
        int _processLHCO_MC;
        std::ofstream fout_MC;

        bool _printLHCO_RECO;
        int _processLHCO_RECO;
        std::ofstream fout_RECO;

        std::string fline00 ;
        std::string del;
        std::string trig;
        std::string fline0;

};

#endif
