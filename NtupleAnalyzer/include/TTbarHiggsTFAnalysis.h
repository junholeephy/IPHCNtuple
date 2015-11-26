#ifndef TTbarHiggsTFAnalysis_H
#define TTbarHiggsTFAnalysis_H

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

class TTbarHiggsTFAnalysis
{

 public:
   TTbarHiggsTFAnalysis();
   TTbarHiggsTFAnalysis(TString inputFileName, TChain *tree, TString theSampleName, TString treeName);

   ~TTbarHiggsTFAnalysis();

   void createHistograms();
   void writeHistograms();
  
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain

   virtual void     Init(TChain *tree);

   virtual void     Loop();
  
   std::vector<Electron>  *vElectron  = new std::vector<Electron>();
   std::vector<Muon>      *vMuon      = new std::vector<Muon>();
   std::vector<Event>     *vEvent     = new std::vector<Event>();
   std::vector<Jet>       *vJet       = new std::vector<Jet>();
   std::vector<Truth>     *vTruth     = new std::vector<Truth>();


 private:

   HistoManager * theHistoManager;

   TFile * outputfile;


   TString sampleName;
   
};

#endif
