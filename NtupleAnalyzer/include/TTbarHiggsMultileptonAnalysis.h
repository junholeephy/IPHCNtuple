#ifndef TTbarHiggsMultileptonAnalysis_H
#define TTbarHiggsMultileptonAnalysis_H

#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TChain.h>
#include <TFile.h>
#include "Muon.h"
#include "Electron.h"
#include "Jet.h"
#include "Event.h"
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
   TTbarHiggsMultileptonAnalysis(TString inputFileName, TChain *tree, TString theSampleName, TString treeName, 
                                 TString outputPath, bool isdata, float xsec, float lumi, int nowe, int nmax);

   ~TTbarHiggsMultileptonAnalysis();

   void createHistograms();
   void writeHistograms();

   void PrintEventList(std::vector<Lepton> leptons,std::vector<Jet> jets);
   
   // needed to print info in LHCO text format (madweight)
   void InitLHCO(int process_MC, int process_RECO);
   void PrintLHCOforMadweight_MC(int evt);
   void PrintLHCOforMadweight_RECO(int evt);
   float Phi_0_2Pi(float phi);
   
   //
   void ThreeLeptonSelection(std::vector<Lepton> vSelectedLeptons,
                             std::vector<Jet>    vSelectedJets,
                             std::vector<Jet>    vSelectedBTagJets,
                             std::vector<Jet>    vSelectedNonBTagJets,
                             int                 nLoose,
                             int                 nMedium,
			     int evt);
   
   bool ThreeLeptonSelection_MC();
   
  
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain

   std::vector<Electron>  *vElectron  = new std::vector<Electron>();
   std::vector<Muon>      *vMuon      = new std::vector<Muon>();
   std::vector<Event>     *vEvent     = new std::vector<Event>();
   std::vector<Jet>       *vJet       = new std::vector<Jet>();
   std::vector<Truth>     *vTruth     = new std::vector<Truth>();

   std::vector<Muon>	 vSelectedMuons;
   std::vector<Electron> vSelectedElectrons;
   std::vector<Lepton>   vSelectedLeptons;
   std::vector<Jet>	 vSelectedNonBTagJets;
   std::vector<Jet>	 vSelectedBTagJets;
   std::vector<Jet>      vSelectedMediumBTagJets;
   std::vector<Jet>	 vSelectedJets;


   virtual void     Init(TChain *tree);

   virtual void     Loop();

   void initializeOutputTree();
   void selectBjets(std::string, int*, int*);
   void fillOutputTree();

   TTree* tOutput;
   Int_t mc_event;
   Float_t mc_weight;
   Int_t mc_3l_category, mc_ttbar_decay, mc_boson_decay, mc_ttZhypAllowed, mc_nJets25, mc_nBtagJets25, mc_nMediumBtagJets25, mc_nNonBtagJets25;
   Int_t multilepton_Bjet1_Id, multilepton_Bjet2_Id;
   Int_t multilepton_Lepton1_Id, multilepton_Lepton2_Id, multilepton_Lepton3_Id;
   Int_t multilepton_JetHighestPt1_Id, multilepton_JetHighestPt2_Id, multilepton_JetClosestMw1_Id, multilepton_JetClosestMw2_Id, multilepton_JetLowestMjj1_Id, multilepton_JetLowestMjj2_Id;
   TLorentzVector multilepton_Bjet1_P4, multilepton_Bjet2_P4;
   TLorentzVector multilepton_Lepton1_P4, multilepton_Lepton2_P4, multilepton_Lepton3_P4;
   TLorentzVector multilepton_JetHighestPt1_P4, multilepton_JetHighestPt2_P4, multilepton_JetClosestMw1_P4, multilepton_JetClosestMw2_P4, multilepton_JetLowestMjj1_P4, multilepton_JetLowestMjj2_P4;
   TLorentzVector multilepton_mET, multilepton_Ptot;





 private:

   HistoManager * theHistoManager;

   TFile * outputfile;

   FILE *fevc;

   TString sampleName;
   TString _outputFile;
   
   bool _isdata;
   float _xsec;
   float _lumi;
   int _nowe;    // number of weighted events
   int _nmax;    // max number of events to process 
 
      
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
