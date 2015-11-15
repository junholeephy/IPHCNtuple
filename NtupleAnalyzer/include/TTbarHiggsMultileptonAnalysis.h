#ifndef TTbarHiggsMultileptonAnalysis_H
#define TTbarHiggsMultileptonAnalysis_H

#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TFile.h>
#include "Muon.h"
#include "Electron.h"
#include "Jet.h"
#include "Event.h"
#include "Truth.h"
#include <TObject.h>
#include <TROOT.h>
#include <iostream>

#include "HistoManager.h"
#include "Lepton.h"

class TTbarHiggsMultileptonAnalysis
{

 public:
   TTbarHiggsMultileptonAnalysis();
   TTbarHiggsMultileptonAnalysis(TString inputfilename, TTree *tree, TString theSampleName);

   ~TTbarHiggsMultileptonAnalysis();

   void createHistograms();
   void writeHistograms();

   void PrintEventList(std::vector<Lepton> leptons,std::vector<Jet> jets);

   void ThreeLeptonSelection(std::vector<Lepton> vSelectedLeptons,
                             std::vector<Jet>    vSelectedJets,
                             std::vector<Jet>    vSelectedBTagJets,
                             std::vector<Jet>    vSelectedNonBTagJets,
                             int                 nLoose,
                             int                 nMedium);

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   std::vector<Electron>             *vElectron             = new std::vector<Electron>();
   std::vector<Muon>                 *vMuon                 = new std::vector<Muon>();
   std::vector<Event>                *vEvent                = new std::vector<Event>();
   std::vector<Jet>                  *vJet                  = new std::vector<Jet>();
   std::vector<Truth>                *vTruth                = new std::vector<Truth>();

        std::vector<Muon>     vSelectedMuons;
        std::vector<Electron> vSelectedElectrons;
        std::vector<Lepton>   vSelectedLeptons;
        std::vector<Jet>      vSelectedNonBTagJets;
        std::vector<Jet>      vSelectedBTagJets;
        std::vector<Jet>      vSelectedJets;


   virtual void     Init(TTree *tree);

   virtual void     Loop();

   void initializeOutputTree();
   void selectBjets(std::string, int*, int*);
   void fillOutputTree();

   TTree* tOutput;
   Int_t mc_event;
   Float_t mc_weight;
   Int_t mc_3l_category, mc_ttbar_decay, mc_boson_decay, mc_ttZhypAllowed, mc_nJets25, mc_nBtagJets25, mc_nNonBtagJets25;
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
};

#endif
