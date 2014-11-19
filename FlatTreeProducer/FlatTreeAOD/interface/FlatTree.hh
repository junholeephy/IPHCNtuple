#ifndef FLATTREE_H
#define FLATTREE_H

#include <TTree.h>
#include <TLorentzVector.h>
#include <string>
#include <iostream>
#include <vector>

#define DEFVAL -666

class FlatTree
{
 public:
   
   FlatTree(TTree* _tree) {tree = _tree;};
   TTree* tree;

   std::map<std::string,bool> conf;
   
   void Init();
   void CreateBranches();
   bool doWrite(std::string name);
   
   int ev_run;
   int ev_id;
   int ev_lumi;
   
   int mc_id;
   int mc_f1;
   int mc_f2;
   float mc_x1;
   float mc_x2;
   float mc_scale;
   float mc_ptHat;
   
   int mc_pu_intime_NumInt;
   int mc_pu_trueNumInt;
   int mc_pu_before_npu;
   int mc_pu_after_npu;

   int mc_pu_Npvi;
   std::vector<int> mc_pu_Nzpositions;
   std::vector<int> mc_pu_BunchCrossing;
   std::vector<std::vector<float> > mc_pu_zpositions;
   std::vector<std::vector<float> > mc_pu_sumpT_lowpT;
   std::vector<std::vector<float> > mc_pu_sumpT_highpT;
   std::vector<std::vector<int> > mc_pu_ntrks_lowpT;
   std::vector<std::vector<int> > mc_pu_ntrks_highpT;

   // Electrons
   
   int el_n;
   std::vector<float> el_pt;
   std::vector<float> el_eta;
   std::vector<float> el_phi;
   std::vector<float> el_m;
   std::vector<float> el_E;
   std::vector<int> el_id;
   std::vector<int> el_charge;
	
   std::vector<float> el_scleta;
   std::vector<int> el_passConversionVeto;
   std::vector<int> el_isGsfCtfScPixChargeConsistent;
   std::vector<float> el_dB3D;
   std::vector<float> el_edB3D;
   
   std::vector<float> el_neutralHadronIso;
   std::vector<float> el_chargedHadronIso;
   std::vector<float> el_puChargedHadronIso;
   std::vector<float> el_ecalIso;
   std::vector<float> el_hcalIso;
   std::vector<float> el_particleIso;
   std::vector<float> el_photonIso;
   std::vector<float> el_trackIso;

   std::vector<float> el_pfIso_sumChargedHadronPt;
   std::vector<float> el_pfIso_sumNeutralHadronEt;
   std::vector<float> el_pfIso_sumPhotonEt;
   std::vector<float> el_pfIso_sumPUPt;
   
   std::vector<int> el_isLoose;
   std::vector<int> el_isTight;
   std::vector<int> el_isRobustLoose;
   std::vector<int> el_isRobustTight;
   std::vector<int> el_isRobustHighEnergy;
   
   std::vector<float> el_vx;
   std::vector<float> el_vy;
   std::vector<float> el_vz;
   
   std::vector<bool> el_isGsf;
   std::vector<float> el_dxy;
   std::vector<float> el_dz;
   std::vector<float> el_dxyError;
   std::vector<float> el_dzError;
   
   std::vector<int> el_numberOfHits;
   
   std::vector<float> el_sigmaIetaIeta;
   std::vector<float> el_hadronicOverEm;
   std::vector<float> el_dr03TkSumPt;
   std::vector<float> el_dr03EcalRecHitSumEt;
   std::vector<float> el_dr03HcalTowerSumEt;
   std::vector<int> el_numberOfLostHits;
   
   // Muons

   int mu_n;
   std::vector<float> mu_pt;
   std::vector<float> mu_eta;
   std::vector<float> mu_phi;
   std::vector<float> mu_m;
   std::vector<float> mu_E;   
   std::vector<int> mu_id;
   std::vector<int> mu_charge;
   
   std::vector<float> mu_dB3D;
   std::vector<float> mu_edB3D;
	
   std::vector<float> mu_neutralHadronIso;
   std::vector<float> mu_chargedHadronIso;
   std::vector<float> mu_puChargedHadronIso;
   std::vector<float> mu_ecalIso;
   std::vector<float> mu_hcalIso;
   std::vector<float> mu_photonIso;
   std::vector<float> mu_trackIso;

   std::vector<float> mu_pfIso03_sumChargedHadronPt;
   std::vector<float> mu_pfIso03_sumNeutralHadronEt;
   std::vector<float> mu_pfIso03_sumPhotonEt;
   std::vector<float> mu_pfIso03_sumPUPt;
   
   std::vector<int> mu_isGlobalMuon;
   std::vector<int> mu_isTrackerMuon;
   std::vector<int> mu_isStandAloneMuon;
   std::vector<int> mu_isCaloMuon;
   std::vector<int> mu_isPFMuon;
   
   std::vector<float> mu_vx;
   std::vector<float> mu_vy;
   std::vector<float> mu_vz;
   
   std::vector<float> mu_globalTrack_dxy;
   std::vector<float> mu_globalTrack_dz;
   std::vector<float> mu_globalTrack_dxyError;
   std::vector<float> mu_globalTrack_dzError;

   std::vector<float> mu_innerTrack_dxy;
   std::vector<float> mu_innerTrack_dz;
   std::vector<float> mu_innerTrack_dxyError;
   std::vector<float> mu_innerTrack_dzError;
   
   std::vector<float> mu_innerTrack_pt;
   std::vector<float> mu_innerTrack_ptError;
   
   std::vector<int> mu_numberOfMatches;
   std::vector<int> mu_numberOfValidMuonHits;

   // Jets

   int jet_n;
   std::vector<float> jet_pt;
   std::vector<float> jet_eta;
   std::vector<float> jet_phi;
   std::vector<float> jet_m;
   std::vector<float> jet_E;
   
   std::vector<float> jet_CSV;
   std::vector<int> jet_flavour;

   std::vector<float> jet_neutralHadronEnergy;
   std::vector<float> jet_neutralEmEnergy;
   std::vector<float> jet_chargedHadronEnergy;
   std::vector<float> jet_chargedEmEnergy;
   std::vector<float> jet_electronEnergy;
   std::vector<float> jet_muonEnergy;
   std::vector<float> jet_photonEnergy;
   
   std::vector<float> jet_pileupJetId;
   
   std::vector<float> jet_gen_pt;
   std::vector<float> jet_gen_eta;
   std::vector<float> jet_gen_phi;
   std::vector<float> jet_gen_m;
   std::vector<float> jet_gen_E;
   
   std::vector<int> jet_gen_status;
   std::vector<int> jet_gen_id;
   
   // ttH
   int mc_truth_tth_channel;
   
   TLorentzVector mc_truth_tth_h0_p4;

   TLorentzVector mc_truth_tth_h0W1_p4;
   TLorentzVector mc_truth_tth_h0W2_p4;
   TLorentzVector mc_truth_tth_h0Wl1_p4;
   TLorentzVector mc_truth_tth_h0Wnu1_p4;
   TLorentzVector mc_truth_tth_h0Wtau1_p4;
   TLorentzVector mc_truth_tth_h0Wnutau1_p4;
   TLorentzVector mc_truth_tth_h0Wtaul1_p4;
   TLorentzVector mc_truth_tth_h0Wtaunu1_p4;
   TLorentzVector mc_truth_tth_h0Wtaunutau1_p4;
   TLorentzVector mc_truth_tth_h0Wl2_p4;
   TLorentzVector mc_truth_tth_h0Wnu2_p4;
   TLorentzVector mc_truth_tth_h0Wtau2_p4;
   TLorentzVector mc_truth_tth_h0Wnutau2_p4;
   TLorentzVector mc_truth_tth_h0Wtaul2_p4;
   TLorentzVector mc_truth_tth_h0Wtaunu2_p4;
   TLorentzVector mc_truth_tth_h0Wtaunutau2_p4;
   TLorentzVector mc_truth_tth_h0Wq11_p4;
   TLorentzVector mc_truth_tth_h0Wq21_p4;
   TLorentzVector mc_truth_tth_h0Wq12_p4;
   TLorentzVector mc_truth_tth_h0Wq22_p4;
   
   TLorentzVector mc_truth_tth_h0Z1_p4;
   TLorentzVector mc_truth_tth_h0Z2_p4;
   TLorentzVector mc_truth_tth_h0Zl11_p4;
   TLorentzVector mc_truth_tth_h0Zl21_p4;
   TLorentzVector mc_truth_tth_h0Ztau11_p4;
   TLorentzVector mc_truth_tth_h0Ztau21_p4;
   TLorentzVector mc_truth_tth_h0Ztaul11_p4;
   TLorentzVector mc_truth_tth_h0Ztaul21_p4;
   TLorentzVector mc_truth_tth_h0Ztaunu11_p4;
   TLorentzVector mc_truth_tth_h0Ztaunu21_p4;
   TLorentzVector mc_truth_tth_h0Ztaunutau11_p4;
   TLorentzVector mc_truth_tth_h0Ztaunutau21_p4;
   TLorentzVector mc_truth_tth_h0Zq11_p4;
   TLorentzVector mc_truth_tth_h0Zq21_p4;
   TLorentzVector mc_truth_tth_h0Zl12_p4;
   TLorentzVector mc_truth_tth_h0Zl22_p4;
   TLorentzVector mc_truth_tth_h0Ztau12_p4;
   TLorentzVector mc_truth_tth_h0Ztau22_p4;
   TLorentzVector mc_truth_tth_h0Ztaul12_p4;
   TLorentzVector mc_truth_tth_h0Ztaul22_p4;
   TLorentzVector mc_truth_tth_h0Ztaunu12_p4;
   TLorentzVector mc_truth_tth_h0Ztaunu22_p4;
   TLorentzVector mc_truth_tth_h0Ztaunutau12_p4;
   TLorentzVector mc_truth_tth_h0Ztaunutau22_p4;
   TLorentzVector mc_truth_tth_h0Zq12_p4;
   TLorentzVector mc_truth_tth_h0Zq22_p4;
   TLorentzVector mc_truth_tth_h0Znu11_p4;
   TLorentzVector mc_truth_tth_h0Znu21_p4;
   TLorentzVector mc_truth_tth_h0Znu12_p4;
   TLorentzVector mc_truth_tth_h0Znu22_p4;
   
   TLorentzVector mc_truth_tth_h0tau1_p4;
   TLorentzVector mc_truth_tth_h0tau2_p4;
   TLorentzVector mc_truth_tth_h0taul1_p4;
   TLorentzVector mc_truth_tth_h0taunutau1_p4;
   TLorentzVector mc_truth_tth_h0taunu1_p4;
   TLorentzVector mc_truth_tth_h0taul2_p4;
   TLorentzVector mc_truth_tth_h0taunutau2_p4;
   TLorentzVector mc_truth_tth_h0taunu2_p4;
   
   TLorentzVector mc_truth_tth_t1_p4;
   TLorentzVector mc_truth_tth_t2_p4;
   TLorentzVector mc_truth_tth_tb1_p4;
   TLorentzVector mc_truth_tth_tb2_p4;
   
   TLorentzVector mc_truth_tth_tW1_p4;
   TLorentzVector mc_truth_tth_tWnu1_p4;
   TLorentzVector mc_truth_tth_tWnutau1_p4;
   TLorentzVector mc_truth_tth_tWl1_p4;
   TLorentzVector mc_truth_tth_tWtau1_p4;
   TLorentzVector mc_truth_tth_tWtaunu1_p4;
   TLorentzVector mc_truth_tth_tWtaunutau1_p4;
   TLorentzVector mc_truth_tth_tWtaul1_p4;
   TLorentzVector mc_truth_tth_tWq11_p4;
   TLorentzVector mc_truth_tth_tWq21_p4;

   TLorentzVector mc_truth_tth_tW2_p4;
   TLorentzVector mc_truth_tth_tWnu2_p4;
   TLorentzVector mc_truth_tth_tWnutau2_p4;
   TLorentzVector mc_truth_tth_tWl2_p4;
   TLorentzVector mc_truth_tth_tWtau2_p4;
   TLorentzVector mc_truth_tth_tWtaunu2_p4;
   TLorentzVector mc_truth_tth_tWtaunutau2_p4;
   TLorentzVector mc_truth_tth_tWtaul2_p4;
   TLorentzVector mc_truth_tth_tWq12_p4;
   TLorentzVector mc_truth_tth_tWq22_p4;

   TLorentzVector mc_truth_tth_j1_p4;
   TLorentzVector mc_truth_tth_j2_p4;
   TLorentzVector mc_truth_tth_j3_p4;

   // pdgId

   int mc_truth_tth_h0_id;

   int mc_truth_tth_h0W1_id;
   int mc_truth_tth_h0W2_id;
   int mc_truth_tth_h0Wl1_id;
   int mc_truth_tth_h0Wnu1_id;
   int mc_truth_tth_h0Wtau1_id;
   int mc_truth_tth_h0Wnutau1_id;
   int mc_truth_tth_h0Wtaul1_id;
   int mc_truth_tth_h0Wtaunu1_id;
   int mc_truth_tth_h0Wtaunutau1_id;
   int mc_truth_tth_h0Wl2_id;
   int mc_truth_tth_h0Wnu2_id;
   int mc_truth_tth_h0Wtau2_id;
   int mc_truth_tth_h0Wnutau2_id;
   int mc_truth_tth_h0Wtaul2_id;
   int mc_truth_tth_h0Wtaunu2_id;
   int mc_truth_tth_h0Wtaunutau2_id;
   int mc_truth_tth_h0Wq11_id;
   int mc_truth_tth_h0Wq21_id;
   int mc_truth_tth_h0Wq12_id;
   int mc_truth_tth_h0Wq22_id;
   
   int mc_truth_tth_h0Z1_id;
   int mc_truth_tth_h0Z2_id;
   int mc_truth_tth_h0Zl11_id;
   int mc_truth_tth_h0Zl21_id;
   int mc_truth_tth_h0Ztau11_id;
   int mc_truth_tth_h0Ztau21_id;
   int mc_truth_tth_h0Ztaul11_id;
   int mc_truth_tth_h0Ztaul21_id;
   int mc_truth_tth_h0Ztaunu11_id;
   int mc_truth_tth_h0Ztaunu21_id;
   int mc_truth_tth_h0Ztaunutau11_id;
   int mc_truth_tth_h0Ztaunutau21_id;
   int mc_truth_tth_h0Zq11_id;
   int mc_truth_tth_h0Zq21_id;
   int mc_truth_tth_h0Zl12_id;
   int mc_truth_tth_h0Zl22_id;
   int mc_truth_tth_h0Ztau12_id;
   int mc_truth_tth_h0Ztau22_id;
   int mc_truth_tth_h0Ztaul12_id;
   int mc_truth_tth_h0Ztaul22_id;
   int mc_truth_tth_h0Ztaunu12_id;
   int mc_truth_tth_h0Ztaunu22_id;
   int mc_truth_tth_h0Ztaunutau12_id;
   int mc_truth_tth_h0Ztaunutau22_id;
   int mc_truth_tth_h0Zq12_id;
   int mc_truth_tth_h0Zq22_id;
   int mc_truth_tth_h0Znu11_id;
   int mc_truth_tth_h0Znu21_id;
   int mc_truth_tth_h0Znu12_id;
   int mc_truth_tth_h0Znu22_id;
   
   int mc_truth_tth_h0tau1_id;
   int mc_truth_tth_h0tau2_id;
   int mc_truth_tth_h0taul1_id;
   int mc_truth_tth_h0taunutau1_id;
   int mc_truth_tth_h0taunu1_id;
   int mc_truth_tth_h0taul2_id;
   int mc_truth_tth_h0taunutau2_id;
   int mc_truth_tth_h0taunu2_id;
   
   int mc_truth_tth_t1_id;
   int mc_truth_tth_t2_id;
   int mc_truth_tth_tb1_id;
   int mc_truth_tth_tb2_id;
   
   int mc_truth_tth_tW1_id;
   int mc_truth_tth_tWnu1_id;
   int mc_truth_tth_tWnutau1_id;
   int mc_truth_tth_tWl1_id;
   int mc_truth_tth_tWtau1_id;
   int mc_truth_tth_tWtaunu1_id;
   int mc_truth_tth_tWtaunutau1_id;
   int mc_truth_tth_tWtaul1_id;
   int mc_truth_tth_tWq11_id;
   int mc_truth_tth_tWq21_id;

   int mc_truth_tth_tW2_id;
   int mc_truth_tth_tWnu2_id;
   int mc_truth_tth_tWnutau2_id;
   int mc_truth_tth_tWl2_id;
   int mc_truth_tth_tWtau2_id;
   int mc_truth_tth_tWtaunu2_id;
   int mc_truth_tth_tWtaunutau2_id;
   int mc_truth_tth_tWtaul2_id;
   int mc_truth_tth_tWq12_id;
   int mc_truth_tth_tWq22_id;

   int mc_truth_tth_j1_id;
   int mc_truth_tth_j2_id;
   int mc_truth_tth_j3_id;

   // status

   int mc_truth_tth_h0_status;

   int mc_truth_tth_h0W1_status;
   int mc_truth_tth_h0W2_status;
   int mc_truth_tth_h0Wl1_status;
   int mc_truth_tth_h0Wnu1_status;
   int mc_truth_tth_h0Wtau1_status;
   int mc_truth_tth_h0Wnutau1_status;
   int mc_truth_tth_h0Wtaul1_status;
   int mc_truth_tth_h0Wtaunu1_status;
   int mc_truth_tth_h0Wtaunutau1_status;
   int mc_truth_tth_h0Wl2_status;
   int mc_truth_tth_h0Wnu2_status;
   int mc_truth_tth_h0Wtau2_status;
   int mc_truth_tth_h0Wnutau2_status;
   int mc_truth_tth_h0Wtaul2_status;
   int mc_truth_tth_h0Wtaunu2_status;
   int mc_truth_tth_h0Wtaunutau2_status;
   int mc_truth_tth_h0Wq11_status;
   int mc_truth_tth_h0Wq21_status;
   int mc_truth_tth_h0Wq12_status;
   int mc_truth_tth_h0Wq22_status;
   
   int mc_truth_tth_h0Z1_status;
   int mc_truth_tth_h0Z2_status;
   int mc_truth_tth_h0Zl11_status;
   int mc_truth_tth_h0Zl21_status;
   int mc_truth_tth_h0Ztau11_status;
   int mc_truth_tth_h0Ztau21_status;
   int mc_truth_tth_h0Ztaul11_status;
   int mc_truth_tth_h0Ztaul21_status;
   int mc_truth_tth_h0Ztaunu11_status;
   int mc_truth_tth_h0Ztaunu21_status;
   int mc_truth_tth_h0Ztaunutau11_status;
   int mc_truth_tth_h0Ztaunutau21_status;
   int mc_truth_tth_h0Zq11_status;
   int mc_truth_tth_h0Zq21_status;
   int mc_truth_tth_h0Zl12_status;
   int mc_truth_tth_h0Zl22_status;
   int mc_truth_tth_h0Ztau12_status;
   int mc_truth_tth_h0Ztau22_status;
   int mc_truth_tth_h0Ztaul12_status;
   int mc_truth_tth_h0Ztaul22_status;
   int mc_truth_tth_h0Ztaunu12_status;
   int mc_truth_tth_h0Ztaunu22_status;
   int mc_truth_tth_h0Ztaunutau12_status;
   int mc_truth_tth_h0Ztaunutau22_status;
   int mc_truth_tth_h0Zq12_status;
   int mc_truth_tth_h0Zq22_status;
   int mc_truth_tth_h0Znu11_status;
   int mc_truth_tth_h0Znu21_status;
   int mc_truth_tth_h0Znu12_status;
   int mc_truth_tth_h0Znu22_status;
   
   int mc_truth_tth_h0tau1_status;
   int mc_truth_tth_h0tau2_status;
   int mc_truth_tth_h0taul1_status;
   int mc_truth_tth_h0taunutau1_status;
   int mc_truth_tth_h0taunu1_status;
   int mc_truth_tth_h0taul2_status;
   int mc_truth_tth_h0taunutau2_status;
   int mc_truth_tth_h0taunu2_status;
   
   int mc_truth_tth_t1_status;
   int mc_truth_tth_t2_status;
   int mc_truth_tth_tb1_status;
   int mc_truth_tth_tb2_status;
   
   int mc_truth_tth_tW1_status;
   int mc_truth_tth_tWnu1_status;
   int mc_truth_tth_tWnutau1_status;
   int mc_truth_tth_tWl1_status;
   int mc_truth_tth_tWtau1_status;
   int mc_truth_tth_tWtaunu1_status;
   int mc_truth_tth_tWtaunutau1_status;
   int mc_truth_tth_tWtaul1_status;
   int mc_truth_tth_tWq11_status;
   int mc_truth_tth_tWq21_status;

   int mc_truth_tth_tW2_status;
   int mc_truth_tth_tWnu2_status;
   int mc_truth_tth_tWnutau2_status;
   int mc_truth_tth_tWl2_status;
   int mc_truth_tth_tWtau2_status;
   int mc_truth_tth_tWtaunu2_status;
   int mc_truth_tth_tWtaunutau2_status;
   int mc_truth_tth_tWtaul2_status;
   int mc_truth_tth_tWq12_status;
   int mc_truth_tth_tWq22_status;

   int mc_truth_tth_j1_status;
   int mc_truth_tth_j2_status;
   int mc_truth_tth_j3_status;   

   // tZq
   int mc_truth_tzq_channel;

   // TLV

   TLorentzVector mc_truth_tzq_Z_p4;
   TLorentzVector mc_truth_tzq_Zl1_p4;
   TLorentzVector mc_truth_tzq_Zl2_p4;
   TLorentzVector mc_truth_tzq_Ztau1_p4;
   TLorentzVector mc_truth_tzq_Ztau2_p4;
   TLorentzVector mc_truth_tzq_Ztaul1_p4;
   TLorentzVector mc_truth_tzq_Ztaul2_p4;
   TLorentzVector mc_truth_tzq_Ztaunu1_p4;
   TLorentzVector mc_truth_tzq_Ztaunu2_p4;
   TLorentzVector mc_truth_tzq_Ztaunutau1_p4;
   TLorentzVector mc_truth_tzq_Ztaunutau2_p4;
   
   TLorentzVector mc_truth_tzq_t_p4;
   TLorentzVector mc_truth_tzq_tb_p4;
   TLorentzVector mc_truth_tzq_tW_p4;
   TLorentzVector mc_truth_tzq_tWnu_p4;
   TLorentzVector mc_truth_tzq_tWnutau_p4;
   TLorentzVector mc_truth_tzq_tWl_p4;
   TLorentzVector mc_truth_tzq_tWtau_p4;
   TLorentzVector mc_truth_tzq_tWtaunu_p4;
   TLorentzVector mc_truth_tzq_tWtaunutau_p4;
   TLorentzVector mc_truth_tzq_tWtaul_p4;
   TLorentzVector mc_truth_tzq_tWq1_p4;
   TLorentzVector mc_truth_tzq_tWq2_p4;

   TLorentzVector mc_truth_tzq_j1_p4;
   TLorentzVector mc_truth_tzq_j2_p4;
   TLorentzVector mc_truth_tzq_j3_p4;

   // pdgId

   int mc_truth_tzq_Z_id;
   int mc_truth_tzq_Zl1_id;
   int mc_truth_tzq_Zl2_id;
   int mc_truth_tzq_Ztau1_id;
   int mc_truth_tzq_Ztau2_id;
   int mc_truth_tzq_Ztaul1_id;
   int mc_truth_tzq_Ztaul2_id;
   int mc_truth_tzq_Ztaunu1_id;
   int mc_truth_tzq_Ztaunu2_id;
   int mc_truth_tzq_Ztaunutau1_id;
   int mc_truth_tzq_Ztaunutau2_id;
   
   int mc_truth_tzq_t_id;
   int mc_truth_tzq_tb_id;
   int mc_truth_tzq_tW_id;
   int mc_truth_tzq_tWnu_id;
   int mc_truth_tzq_tWnutau_id;
   int mc_truth_tzq_tWl_id;
   int mc_truth_tzq_tWtau_id;
   int mc_truth_tzq_tWtaunu_id;
   int mc_truth_tzq_tWtaunutau_id;
   int mc_truth_tzq_tWtaul_id;
   int mc_truth_tzq_tWq1_id;
   int mc_truth_tzq_tWq2_id;
   
   int mc_truth_tzq_j1_id;
   int mc_truth_tzq_j2_id;
   int mc_truth_tzq_j3_id;

   // status
   
   int mc_truth_tzq_Z_status;
   int mc_truth_tzq_Zl1_status;
   int mc_truth_tzq_Zl2_status;
   int mc_truth_tzq_Ztau1_status;
   int mc_truth_tzq_Ztau2_status;
   int mc_truth_tzq_Ztaul1_status;
   int mc_truth_tzq_Ztaul2_status;
   int mc_truth_tzq_Ztaunu1_status;
   int mc_truth_tzq_Ztaunu2_status;
   int mc_truth_tzq_Ztaunutau1_status;
   int mc_truth_tzq_Ztaunutau2_status;
   
   int mc_truth_tzq_t_status;
   int mc_truth_tzq_tb_status;
   int mc_truth_tzq_tW_status;
   int mc_truth_tzq_tWnu_status;
   int mc_truth_tzq_tWnutau_status;
   int mc_truth_tzq_tWl_status;
   int mc_truth_tzq_tWtau_status;
   int mc_truth_tzq_tWtaunu_status;
   int mc_truth_tzq_tWtaunutau_status;
   int mc_truth_tzq_tWtaul_status;
   int mc_truth_tzq_tWq1_status;
   int mc_truth_tzq_tWq2_status;
   
   int mc_truth_tzq_j1_status;
   int mc_truth_tzq_j2_status;
   int mc_truth_tzq_j3_status;
};

#endif
