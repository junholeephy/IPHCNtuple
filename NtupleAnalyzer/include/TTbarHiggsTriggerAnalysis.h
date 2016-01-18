#ifndef TTbarHiggsTriggerAnalysis_H
#define TTbarHiggsTriggerAnalysis_H

#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TChain.h>
#include <TFile.h>
#include "Muon.h"
#include "Electron.h"
#include "Jet.h"
#include "TriggerObj.h"
#include "Event.h"
#include "Truth.h"
#include <TObject.h>
#include <TROOT.h>
#include <iostream>
#include <fstream>
#include "Base.h"

#include "HistoManager.h"
#include "Lepton.h"

class TTbarHiggsTriggerAnalysis
{

 public:
   TTbarHiggsTriggerAnalysis();
   TTbarHiggsTriggerAnalysis(TString inputFileName, TChain *tree, TString sampleName, TString treeName, TString outputFileName,
                             bool isdata, float xsec, float lumi, int nowe, int nmax);

   ~TTbarHiggsTriggerAnalysis();

   void createHistograms();
   void writeHistograms();
  
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain

   virtual void     Init(TChain *tree);

   virtual void     Loop();
  
   std::vector<Electron>   *vElectron   = new std::vector<Electron>();
   std::vector<Muon>       *vMuon       = new std::vector<Muon>();
   std::vector<Event>      *vEvent      = new std::vector<Event>();
   std::vector<Jet>        *vJet        = new std::vector<Jet>();
   std::vector<Truth>      *vTruth      = new std::vector<Truth>();
   std::vector<TriggerObj> *vTriggerObj = new std::vector<TriggerObj>();
   
   std::vector<Muon>	   vSelectedMuons;
   std::vector<Electron>   vSelectedElectrons;
   std::vector<Lepton>     vSelectedLeptons;
   
   std::vector<TriggerObj> vSelectedTriggerObj_IsoMu20;
   std::vector<int>        vSelectedTriggerObj_IsoMu20_recoMatched;
   std::vector<TriggerObj> vSelectedTriggerObj_Ele23;
   std::vector<int>        vSelectedTriggerObj_Ele23_recoMatched;
   
   std::vector<TriggerObj> vSelectedTriggerObj_Ele17Ele12Leg1;
   std::vector<TriggerObj> vSelectedTriggerObj_Ele17Ele12Leg2;
   
  
 private:

   HistoManager * theHistoManager;

   TFile * _outputFile;
   TString _outputFileName;

   TString _sampleName;
  
   bool  _isdata;
   float _xsec;
   float _lumi;
   int   _nowe; // number of weighted events
   int   _nmax; // max number of events to process
     
};

#endif
