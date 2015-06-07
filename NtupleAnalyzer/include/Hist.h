#ifndef HIST_H
#define HIST_H

#include "Analyzer.h"

//#include "FakeWeight.h"

class Hist
{
   
 public:
   
   Hist();
   virtual ~Hist();
   
   void setElectron(std::vector<Electron> *v)                  {_v_Electron = v;};
   void setMuon(std::vector<Muon> *v)                      {_v_Muon = v;};
   void setJet(std::vector<Jet> *v)                      {_v_Jet = v;};
   void setEvent(std::vector<Event> *v)                      {_v_Event = v;};

//   void setFakeWeight(SKYPLOT::FakeWeight fakeWeight)      {_fakeWeight = fakeWeight;}; 
	
   void init();
   void fill();
   void close();
   bool printout(bool doPrint);

   std::vector<int> filterPt(double Pt1, double Pt2, double Pt3,
			     std::vector<Electron>* ntElectron,
			     std::vector<Muon>* ntMuon);

   unsigned int passSel;
   
   std::vector<TLorentzVector> lep_tlv(std::string sys,int llc,bool isTIGHT);
   std::vector<TLorentzVector> jet_tlv(std::string sys); 
   
//   float getVar(std::string sys,int ijet,std::string varName,int ibin);
//   std::vector<float> getVarVec(std::string sys,int ijet,std::string varName,int ibin);
//   std::pair<float,float> getVar2d(std::string sys,int ijet,std::string varName,int ibin);
//   std::vector<float> getVar3d(std::string sys,int ijet,std::string varName,int ibin);
//   float getPt(std::string sys);
//   double DeltaR(double eta1,double phi1,double eta2,double phi2);

//   double getSF(std::vector<std::pair<double,double> > vBin,
//		std::vector<std::pair<double,double> > vSf,
//		double var);
//   double get2DSF(std::vector<std::pair<double,double> > vBinX,
//		  std::vector<std::pair<double,double> > vBinY,
//		  std::vector<std::pair<double,double> > vSf,
//		  double varX,
//		  double varY);
	
   std::string dilchan[100];                                                                                           
   std::string dilchar[100];                                                                                           
   std::string jets[100];                                                                                              
   std::string sys[100];                                                                                               
   std::string sys_low[100];                                                                                           
   std::string sys_up[100];                                                                                            
   std::string histname_dilep[100];                                                                                    
   std::string histname_lep[100];                                                                                      
   std::string histname_jet[100];                                                                                      
   std::string type[100];                                                                                              
   std::string sel[100];                                                                                               
   
   std::string histname_lep_2d[100];                                                                                   
   
   std::vector<std::string> hname;
	
   // var / ptbin / bins / x1,x2
//   std::vector<std::vector<std::vector<std::pair<double,double> > > > rwBin;
//   std::vector<std::vector<std::vector<std::pair<double,double> > > > rwSf;
//   std::vector<std::vector<std::vector<std::pair<double,double> > > > rwBin_btag;
//   std::vector<std::vector<std::vector<std::pair<double,double> > > > rwSf_btag;
   
   // var[0] / ptbin / bins / v1,v2
//   std::vector<std::vector<std::vector<std::pair<double,double> > > > rw3DBinX;
//   std::vector<std::vector<std::vector<std::pair<double,double> > > > rw3DBinY;
//   std::vector<std::vector<std::vector<std::pair<double,double> > > > rw3DBinZ;
//   std::vector<std::vector<std::vector<std::pair<double,double> > > > rw3DSf;
   
   // var[0] / ptbin / bins / v1,v2
//   std::vector<std::vector<std::vector<std::pair<double,double> > > > rw2DBinX;
//   std::vector<std::vector<std::vector<std::pair<double,double> > > > rw2DBinY;
//   std::vector<std::vector<std::vector<std::pair<double,double> > > > rw2DSf;
	
   int sel_n;
   int type_n;
   int sys_n;                                                                                                          
   int sys_low_n;
   int sys_up_n;
   int n_jets;
	
   bool fillThis;

   // Tree variables
   double m_weight;
   int m_channel;
   int m_sel;
   
   int m_l1_type;
   int m_l1_charge;
   double m_l1_pt;
   double m_l1_eta;
   double m_l1_phi;
   double m_l1_m;
   
   int m_l2_type;
   int m_l2_charge;
   double m_l2_pt;
   double m_l2_eta;
   double m_l2_phi;
   double m_l2_m;

   int m_l3_type;
   int m_l3_charge;
   double m_l3_pt;
   double m_l3_eta;
   double m_l3_phi;
   double m_l3_m;

   int m_l4_type;
   int m_l4_charge;
   double m_l4_pt;
   double m_l4_eta;
   double m_l4_phi;
   double m_l4_m;
   
   int m_j1_ntrk;
   double m_j1_btag;
   double m_j1_pt;
   double m_j1_eta;
   double m_j1_phi;
   double m_j1_m;

   int m_j2_ntrk;
   double m_j2_btag;
   double m_j2_pt;
   double m_j2_eta;
   double m_j2_phi;
   double m_j2_m;

   int m_j3_ntrk;
   double m_j3_btag;
   double m_j3_pt;
   double m_j3_eta;
   double m_j3_phi;
   double m_j3_m;

   int m_j4_ntrk;
   double m_j4_btag;
   double m_j4_pt;
   double m_j4_eta;
   double m_j4_phi;
   double m_j4_m;

   int m_j5_ntrk;
   double m_j5_btag;
   double m_j5_pt;
   double m_j5_eta;
   double m_j5_phi;
   double m_j5_m;

   int m_j6_ntrk;
   double m_j6_btag;
   double m_j6_pt;
   double m_j6_eta;
   double m_j6_phi;
   double m_j6_m;
   
   double m_met;
   double m_met_phi;
   
 protected:

   std::vector<Electron>             *_v_Electron;
   std::vector<Muon>             *_v_Muon;
   std::vector<Electron>             *_v_ElectronTight;
   std::vector<Muon>             *_v_MuonTight;
   
   std::vector<Event>             *_v_Event;
   std::vector<Jet>             *_v_Jet;
   std::vector<Jet>             *_v_JetTight;
   
//   SKYPLOT::FakeWeight           _fakeWeight;
  
   std::map<std::string, TH1D*> *_m1d_Event;
   std::map<std::string, TH1D*> *_m1d_Dilepton;
   std::map<std::string, TH1D*> *_m1d_Lepton;
   std::map<std::string, TH1D*> *_m1d_Jet;

   std::vector<std::pair<std::vector<std::string>,double*> > *_s_Dilepton;
	
 private:

   // [FLAVOUR][PT][ETA][BTAG][ADDSEL][VAR][2*(NSYS-1)+1]
   std::string histNAMES[5][19][1][10][1][10][25];
   std::string histNAMES_2d[5][19][1][10][1][10][25];
   std::string histNAMES_3d[5][19][1][10][1][10][25];
   
//   JetCorrectionUncertainty *jesTotal;
   
   double cJER[5];
   double cJER_down[5];
   double cJER_up[5];
	
//	TLorentzVector *v_jet;
//	TLorentzVector *v_jet_sys_jesTotalUp;
//	TLorentzVector *v_jet_sys_jesTotalLow;
//	TLorentzVector *v_jet_sys_jerTotalUp;
//	TLorentzVector *v_jet_sys_jerTotalLow;
//	TLorentzVector *v_mu;
	
   FILE *_fevc;
   std::ofstream _fevcVal;
	
   TRandom3 *rnd;
   
   TFile *_fout;
   TTree *_trout;
};

#endif
