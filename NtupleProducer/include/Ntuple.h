#ifndef NTUPLE_H
#define NTUPLE_H

#include "Electron.h"
#include "Muon.h"
#include "Event.h"
#include "Jet.h"
#include "Truth.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

class Ntuple
{
 public:
   Ntuple();
   virtual ~Ntuple();
   
   std::vector<Electron>* NtElectron;
   std::vector<Muon>* NtMuon;
   std::vector<Event>* NtEvent;
   std::vector<Jet>* NtJet;
   std::vector<Truth>* NtTruth;
   
   void Init();
   
   void setBranchAddress();
   void createVar();
   void clearVar();
   void fill();
   
   TFile*  m_file;
   
 private:
   TTree*  m_tree;
   TChain* m_chain;
};

#endif
