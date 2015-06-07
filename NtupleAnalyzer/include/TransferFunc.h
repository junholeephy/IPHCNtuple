#ifndef TRANSFERFUNC_H
#define TRANSFERFUNC_H

#include "Analyzer.h"

class TransferFunc
{
   
 public:
   
   TransferFunc();
   virtual ~TransferFunc();
   
   void setElectron(std::vector<Electron> *v)                  {_v_Electron = v;};
   void setMuon(std::vector<Muon> *v)                      {_v_Muon = v;};
   void setJet(std::vector<Jet> *v)                      {_v_Jet = v;};
   void setEvent(std::vector<Event> *v)                      {_v_Event = v;};
   void setTruth(std::vector<Truth> *v)                      {_v_Truth = v;};

   float GetDeltaR(float eta1,float phi1,float eta2,float phi2);
   
   void init();
   void run();
   void close();
   
 protected:

   std::vector<Electron>             *_v_Electron;
   std::vector<Muon>             *_v_Muon;
   std::vector<Electron>             *_v_ElectronTight;
   std::vector<Muon>             *_v_MuonTight;
   std::vector<Truth>             *_v_Truth;
   
   std::vector<Event>             *_v_Event;
   std::vector<Jet>             *_v_Jet;
   std::vector<Jet>             *_v_JetTight;
   
 private:

   int m_elTruth_n;
   float m_elTruth_pt[1000];
   float m_elTruth_E[1000];
   float m_elTruth_eta[1000];
   int m_elTruth_label[1000];
   float m_elRec_pt[1000];
   float m_elRec_E[1000];
   float m_elRec_eta[1000];

   int m_muTruth_n;
   float m_muTruth_pt[1000];
   float m_muTruth_E[1000];
   float m_muTruth_eta[1000];
   int m_muTruth_label[1000];
   float m_muRec_pt[1000];
   float m_muRec_E[1000];
   float m_muRec_eta[1000];

   int m_qTruth_n;
   float m_qTruth_pt[1000];
   float m_qTruth_E[1000];
   float m_qTruth_eta[1000];
   int m_qTruth_label[1000];
   float m_qRec_pt[1000];
   float m_qRec_E[1000];
   float m_qRec_eta[1000];

   int m_bTruth_n;
   float m_bTruth_pt[1000];
   float m_bTruth_E[1000];
   float m_bTruth_eta[1000];
   int m_bTruth_label[1000];
   float m_bRec_pt[1000];
   float m_bRec_E[1000];
   float m_bRec_eta[1000];
   
   TRandom3 *rnd;
   
   TFile *_fout;
   TTree *_tr;
};

#endif
