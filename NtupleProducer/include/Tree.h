#ifndef Tree_h
#define Tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>

#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

class Tree {
public :
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           ev_run;
   Int_t           ev_id;
   Int_t           ev_lumi;
   Float_t         ev_rho;
   Float_t         met_pt;
   Float_t         met_phi;
   Float_t         met_sumet;
   Float_t         pv_x;
   Float_t         pv_y;
   Float_t         pv_z;
   Int_t           mc_id;
   Int_t           mc_f1;
   Int_t           mc_f2;
   Float_t         mc_x1;
   Float_t         mc_x2;
   Float_t         mc_scale;
   Float_t         mc_ptHat;
   Float_t         mc_weight;
   Int_t           mc_pu_intime_NumInt;
   Int_t           mc_pu_trueNumInt;
   Int_t           mc_pu_before_npu;
   Int_t           mc_pu_after_npu;
   Int_t           mc_pu_Npvi;
   vector<int>     *mc_pu_Nzpositions;
   vector<int>     *mc_pu_BunchCrossing;
   vector<vector<float> > *mc_pu_zpositions;
   vector<vector<float> > *mc_pu_sumpT_lowpT;
   vector<vector<float> > *mc_pu_sumpT_highpT;
   vector<vector<int> > *mc_pu_ntrks_lowpT;
   vector<vector<int> > *mc_pu_ntrks_highpT;
   Int_t           el_n;
   vector<float>   *el_pt;
   vector<float>   *el_eta;
   vector<float>   *el_phi;
   vector<float>   *el_m;
   vector<float>   *el_E;
   vector<int>     *el_id;
   vector<int>     *el_charge;
   vector<float>   *el_neutralHadronIso;
   vector<float>   *el_chargedHadronIso;
   vector<float>   *el_puChargedHadronIso;
   vector<float>   *el_ecalIso;
   vector<float>   *el_hcalIso;
   vector<float>   *el_particleIso;
   vector<float>   *el_photonIso;
   vector<float>   *el_trackIso;
   vector<int>     *el_isLoose;
   vector<int>     *el_isTight;
   vector<int>     *el_isRobustLoose;
   vector<int>     *el_isRobustTight;
   vector<int>     *el_isRobustHighEnergy;
   vector<float>   *el_vx;
   vector<float>   *el_vy;
   vector<float>   *el_vz;
   vector<bool>    *el_isGsf;
   vector<float>   *el_dxy;
   vector<float>   *el_dz;
   vector<float>   *el_dxyError;
   vector<float>   *el_dzError;
   vector<float>   *el_mvaNonTrigV0;
   vector<int>     *el_numberOfHits;
   vector<float>   *el_pfIso_sumChargedHadronPt;
   vector<float>   *el_pfIso_sumNeutralHadronEt;
   vector<float>   *el_pfIso_sumPhotonEt;
   vector<float>   *el_pfIso_sumPUPt;
   vector<float>   *el_lepMVA;
   vector<float>   *el_lepMVA_neuRelIso;
   vector<float>   *el_lepMVA_chRelIso;
   vector<float>   *el_lepMVA_jetDR;
   vector<float>   *el_lepMVA_jetPtRatio;
   vector<float>   *el_lepMVA_jetBTagCSV;
   vector<float>   *el_lepMVA_sip3d;
   vector<float>   *el_lepMVA_dxy;
   vector<float>   *el_lepMVA_dz;
   vector<float>   *el_lepMVA_mvaId;
   vector<int>     *el_isGsfCtfScPixChargeConsistent;
   vector<int>     *el_passConversionVeto;
   vector<float>   *el_deltaEtaSuperClusterTrackAtVtx;
   vector<float>   *el_deltaPhiSuperClusterTrackAtVtx;
   vector<float>   *el_see;
   vector<float>   *el_hadronicOverEm;
   vector<float>   *el_scleta;
   vector<float>   *el_dB3D;
   vector<float>   *el_edB3D;
   vector<bool>    *el_hasMatchedConversion;
   Int_t           mu_n;
   vector<float>   *mu_pt;
   vector<float>   *mu_eta;
   vector<float>   *mu_phi;
   vector<float>   *mu_m;
   vector<float>   *mu_E;
   vector<int>     *mu_id;
   vector<int>     *mu_charge;
   vector<float>   *mu_neutralHadronIso;
   vector<float>   *mu_chargedHadronIso;
   vector<float>   *mu_ecalIso;
   vector<float>   *mu_hcalIso;
   vector<float>   *mu_photonIso;
   vector<float>   *mu_trackIso;
   vector<int>     *mu_isGlobalMuon;
   vector<int>     *mu_isTrackerMuon;
   vector<int>     *mu_isStandAloneMuon;
   vector<int>     *mu_isCaloMuon;
   vector<int>     *mu_isPFMuon;
   vector<bool>    *mu_isTightMuon;
   vector<float>   *mu_vx;
   vector<float>   *mu_vy;
   vector<float>   *mu_vz;
   vector<float>   *mu_segmentCompatibility;
   vector<int>     *mu_hasGlobalTrack;
   vector<float>   *mu_globalTrack_dxy;
   vector<float>   *mu_globalTrack_dz;
   vector<float>   *mu_globalTrack_dxyError;
   vector<float>   *mu_globalTrack_dzError;
   vector<float>   *mu_globalTrack_normalizedChi2;
   vector<float>   *mu_combinedQuality_chi2LocalPosition;
   vector<float>   *mu_combinedQuality_trkKink;
   vector<int>     *mu_hasInnerTrack;
   vector<float>   *mu_innerTrack_dxy;
   vector<float>   *mu_innerTrack_dz;
   vector<float>   *mu_innerTrack_dxyError;
   vector<float>   *mu_innerTrack_dzError;
   vector<float>   *mu_innerTrack_validFraction;
   vector<float>   *mu_bestTrack_dxy;
   vector<float>   *mu_bestTrack_dz;
   vector<float>   *mu_bestTrack_dxyError;
   vector<float>   *mu_bestTrack_dzError;
   vector<int>     *mu_numberOfMatches;
   vector<int>     *mu_numberOfValidMuonHits;
   vector<float>   *mu_pfIso03_sumChargedHadronPt;
   vector<float>   *mu_pfIso03_sumNeutralHadronEt;
   vector<float>   *mu_pfIso03_sumPhotonEt;
   vector<float>   *mu_pfIso03_sumPUPt;
   vector<float>   *mu_lepMVA;
   vector<float>   *mu_lepMVA_neuRelIso;
   vector<float>   *mu_lepMVA_chRelIso;
   vector<float>   *mu_lepMVA_jetDR;
   vector<float>   *mu_lepMVA_jetPtRatio;
   vector<float>   *mu_lepMVA_jetBTagCSV;
   vector<float>   *mu_lepMVA_sip3d;
   vector<float>   *mu_lepMVA_dxy;
   vector<float>   *mu_lepMVA_dz;
   vector<float>   *mu_lepMVA_mvaId;
   vector<float>   *mu_innerTrack_pt;
   vector<float>   *mu_innerTrack_ptError;
   vector<float>   *mu_dB3D;
   vector<float>   *mu_edB3D;
   Int_t           jet_n;
   vector<float>   *jet_pt;
   vector<float>   *jet_eta;
   vector<float>   *jet_phi;
   vector<float>   *jet_m;
   vector<float>   *jet_E;
   vector<int>     *jet_ntrk;
   vector<float>   *jet_CSV;
   vector<float>   *jet_CSVv2;
   vector<int>     *jet_partonFlavour;
   vector<int>     *jet_hadronFlavour;
   vector<float>   *jet_neutralHadronEnergy;
   vector<float>   *jet_neutralEmEnergy;
   vector<float>   *jet_chargedHadronEnergy;
   vector<float>   *jet_chargedEmEnergy;
   vector<float>   *jet_electronEnergy;
   vector<float>   *jet_muonEnergy;
   vector<float>   *jet_photonEnergy;
   vector<float>   *jet_genJet_pt;
   vector<float>   *jet_genJet_eta;
   vector<float>   *jet_genJet_phi;
   vector<float>   *jet_genJet_m;
   vector<float>   *jet_genJet_E;
   vector<int>     *jet_genJet_status;
   vector<int>     *jet_genJet_id;
   vector<float>   *jet_pileupJetId;

   Int_t gen_n;
   std::vector<float> *gen_pt;
   std::vector<float> *gen_eta;
   std::vector<float> *gen_phi;
   std::vector<float> *gen_m;
   std::vector<int> *gen_id;
   std::vector<int> *gen_status;
   std::vector<int> *gen_index;
   std::vector<int> *gen_mother_index;
   
   Int_t           mc_truth_h0_id;
   Int_t           mc_truth_h0W1_id;
   Int_t           mc_truth_h0Wl1_id;
   Int_t           mc_truth_h0Wq11_id;
   Int_t           mc_truth_h0Wq21_id;
   Int_t           mc_truth_h0W2_id;
   Int_t           mc_truth_h0Wl2_id;
   Int_t           mc_truth_h0Wq12_id;
   Int_t           mc_truth_h0Wq22_id;
   Int_t           mc_truth_h0Z1_id;
   Int_t           mc_truth_h0Zl11_id;
   Int_t           mc_truth_h0Zl21_id;
   Int_t           mc_truth_h0Zq11_id;
   Int_t           mc_truth_h0Zq21_id;
   Int_t           mc_truth_h0Z2_id;
   Int_t           mc_truth_h0Zl12_id;
   Int_t           mc_truth_h0Zl22_id;
   Int_t           mc_truth_h0Zq12_id;
   Int_t           mc_truth_h0Zq22_id;
   Int_t           mc_truth_h0tau1_id;
   Int_t           mc_truth_h0tau2_id;
   Int_t           mc_truth_h0b1_id;
   Int_t           mc_truth_h0b2_id;
   
   Int_t           mc_truth_t1_id;
   Int_t           mc_truth_t2_id;
   Int_t           mc_truth_tb1_id;
   Int_t           mc_truth_tb2_id;
   Int_t           mc_truth_tW1_id;
   Int_t           mc_truth_tW2_id;
   Int_t           mc_truth_tWl1_id;
   Int_t           mc_truth_tWl2_id;
   Int_t           mc_truth_tWq11_id;
   Int_t           mc_truth_tWq21_id;
   Int_t           mc_truth_tWq12_id;
   Int_t           mc_truth_tWq22_id;

   Int_t           mc_truth_t_id;
   Int_t           mc_truth_tb_id;
   Int_t           mc_truth_tW_id;
   Int_t           mc_truth_tWl_id;
   Int_t           mc_truth_tWq1_id;
   Int_t           mc_truth_tWq2_id;
   
   TLorentzVector           *mc_truth_t1_p4;
   TLorentzVector           *mc_truth_t2_p4;
   TLorentzVector           *mc_truth_tb1_p4;
   TLorentzVector           *mc_truth_tb2_p4;
   TLorentzVector           *mc_truth_tW1_p4;
   TLorentzVector           *mc_truth_tW2_p4;
   TLorentzVector           *mc_truth_tWl1_p4;
   TLorentzVector           *mc_truth_tWl2_p4;
   TLorentzVector           *mc_truth_tWq11_p4;
   TLorentzVector           *mc_truth_tWq21_p4;
   TLorentzVector           *mc_truth_tWq12_p4;
   TLorentzVector           *mc_truth_tWq22_p4;

   TLorentzVector           *mc_truth_t_p4;
   TLorentzVector           *mc_truth_tb_p4;
   TLorentzVector           *mc_truth_tW_p4;
   TLorentzVector           *mc_truth_tWl_p4;
   TLorentzVector           *mc_truth_tWq1_p4;
   TLorentzVector           *mc_truth_tWq2_p4;
   
   TLorentzVector           *mc_truth_h0_p4;
   TLorentzVector           *mc_truth_h0W1_p4;
   TLorentzVector           *mc_truth_h0Wl1_p4;
   TLorentzVector           *mc_truth_h0Wq11_p4;
   TLorentzVector           *mc_truth_h0Wq21_p4;
   TLorentzVector           *mc_truth_h0W2_p4;
   TLorentzVector           *mc_truth_h0Wl2_p4;
   TLorentzVector           *mc_truth_h0Wq12_p4;
   TLorentzVector           *mc_truth_h0Wq22_p4;
   TLorentzVector           *mc_truth_h0Z1_p4;
   TLorentzVector           *mc_truth_h0Zl11_p4;
   TLorentzVector           *mc_truth_h0Zl21_p4;
   TLorentzVector           *mc_truth_h0Zq11_p4;
   TLorentzVector           *mc_truth_h0Zq21_p4;
   TLorentzVector           *mc_truth_h0Z2_p4;
   TLorentzVector           *mc_truth_h0Zl12_p4;
   TLorentzVector           *mc_truth_h0Zl22_p4;
   TLorentzVector           *mc_truth_h0Zq12_p4;
   TLorentzVector           *mc_truth_h0Zq22_p4;
   TLorentzVector           *mc_truth_h0tau1_p4;
   TLorentzVector           *mc_truth_h0tau2_p4;
   TLorentzVector           *mc_truth_h0b1_p4;
   TLorentzVector           *mc_truth_h0b2_p4;
   
   // List of branches
   TBranch        *b_ev_run;   //!
   TBranch        *b_ev_id;   //!
   TBranch        *b_ev_lumi;   //!
   TBranch        *b_ev_rho;   //!
   TBranch        *b_met_pt;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_sumet;   //!
   TBranch        *b_pv_x;   //!
   TBranch        *b_pv_y;   //!
   TBranch        *b_pv_z;   //!
   TBranch        *b_mc_weight;   //!
   TBranch        *b_mc_id;   //!
   TBranch        *b_mc_f1;   //!
   TBranch        *b_mc_f2;   //!
   TBranch        *b_mc_x1;   //!
   TBranch        *b_mc_x2;   //!
   TBranch        *b_mc_scale;   //!
   TBranch        *b_mc_ptHat;   //!
   TBranch        *b_mc_pu_intime_NumInt;   //!
   TBranch        *b_mc_pu_trueNumInt;   //!
   TBranch        *b_mc_pu_before_npu;   //!
   TBranch        *b_mc_pu_after_npu;   //!
   TBranch        *b_mc_pu_Npvi;   //!
   TBranch        *b_mc_pu_Nzpositions;   //!
   TBranch        *b_mc_pu_BunchCrossing;   //!
   TBranch        *b_mc_pu_zpositions;   //!
   TBranch        *b_mc_pu_sumpT_lowpT;   //!
   TBranch        *b_mc_pu_sumpT_highpT;   //!
   TBranch        *b_mc_pu_ntrks_lowpT;   //!
   TBranch        *b_mc_pu_ntrks_highpT;   //!
   TBranch        *b_el_n;   //!
   TBranch        *b_el_pt;   //!
   TBranch        *b_el_eta;   //!
   TBranch        *b_el_phi;   //!
   TBranch        *b_el_m;   //!
   TBranch        *b_el_E;   //!
   TBranch        *b_el_id;   //!
   TBranch        *b_el_charge;   //!
   TBranch        *b_el_neutralHadronIso;   //!
   TBranch        *b_el_chargedHadronIso;   //!
   TBranch        *b_el_puChargedHadronIso;   //!
   TBranch        *b_el_ecalIso;   //!
   TBranch        *b_el_hcalIso;   //!
   TBranch        *b_el_particleIso;   //!
   TBranch        *b_el_photonIso;   //!
   TBranch        *b_el_trackIso;   //!
   TBranch        *b_el_isLoose;   //!
   TBranch        *b_el_isTight;   //!
   TBranch        *b_el_isRobustLoose;   //!
   TBranch        *b_el_isRobustTight;   //!
   TBranch        *b_el_isRobustHighEnergy;   //!
   TBranch        *b_el_vx;   //!
   TBranch        *b_el_vy;   //!
   TBranch        *b_el_vz;   //!
   TBranch        *b_el_isGsf;   //!
   TBranch        *b_el_dxy;   //!
   TBranch        *b_el_dz;   //!
   TBranch        *b_el_dxyError;   //!
   TBranch        *b_el_dzError;   //!
   TBranch        *b_el_mvaNonTrigV0;   //!
   TBranch        *b_el_numberOfHits;   //!
   TBranch        *b_el_pfIso_sumChargedHadronPt;   //!
   TBranch        *b_el_pfIso_sumNeutralHadronEt;   //!
   TBranch        *b_el_pfIso_sumPhotonEt;   //!
   TBranch        *b_el_pfIso_sumPUPt;   //!
   TBranch        *b_el_lepMVA;   //!
   TBranch        *b_el_lepMVA_neuRelIso;
   TBranch        *b_el_lepMVA_chRelIso;
   TBranch        *b_el_lepMVA_jetDR;
   TBranch        *b_el_lepMVA_jetPtRatio;
   TBranch        *b_el_lepMVA_jetBTagCSV;
   TBranch        *b_el_lepMVA_sip3d;
   TBranch        *b_el_lepMVA_dxy;
   TBranch        *b_el_lepMVA_dz;
   TBranch        *b_el_lepMVA_mvaId;
   TBranch        *b_el_isGsfCtfScPixChargeConsistent;   //!
   TBranch        *b_el_passConversionVeto;   //!
   TBranch        *b_el_deltaEtaSuperClusterTrackAtVtx;
   TBranch        *b_el_deltaPhiSuperClusterTrackAtVtx;
   TBranch        *b_el_see;
   TBranch        *b_el_hadronicOverEm;
   TBranch        *b_el_scleta;
   TBranch        *b_el_dB3D;
   TBranch        *b_el_edB3D;
   TBranch        *b_el_hasMatchedConversion;
   TBranch        *b_mu_n;   //!
   TBranch        *b_mu_pt;   //!
   TBranch        *b_mu_eta;   //!
   TBranch        *b_mu_phi;   //!
   TBranch        *b_mu_m;   //!
   TBranch        *b_mu_E;   //!
   TBranch        *b_mu_id;   //!
   TBranch        *b_mu_charge;   //!
   TBranch        *b_mu_neutralHadronIso;   //!
   TBranch        *b_mu_chargedHadronIso;   //!
   TBranch        *b_mu_ecalIso;   //!
   TBranch        *b_mu_hcalIso;   //!
   TBranch        *b_mu_photonIso;   //!
   TBranch        *b_mu_trackIso;   //!
   TBranch        *b_mu_isGlobalMuon;   //!
   TBranch        *b_mu_isTrackerMuon;   //!
   TBranch        *b_mu_isStandAloneMuon;   //!
   TBranch        *b_mu_isCaloMuon;   //!
   TBranch        *b_mu_isPFMuon;   //!
   TBranch        *b_mu_isTightMuon;   //!
   TBranch        *b_mu_vx;   //!
   TBranch        *b_mu_vy;   //!
   TBranch        *b_mu_vz;   //!
   TBranch        *b_mu_segmentCompatibility;   //!
   TBranch        *b_mu_hasGlobalTrack;   //!
   TBranch        *b_mu_globalTrack_dxy;   //!
   TBranch        *b_mu_globalTrack_dz;   //!
   TBranch        *b_mu_globalTrack_dxyError;   //!
   TBranch        *b_mu_globalTrack_dzError;   //!
   TBranch        *b_mu_globalTrack_normalizedChi2;   //!
   TBranch        *b_mu_combinedQuality_chi2LocalPosition;   //!
   TBranch        *b_mu_combinedQuality_trkKink;   //!
   TBranch        *b_mu_hasInnerTrack;   //!
   TBranch        *b_mu_innerTrack_dxy;   //!
   TBranch        *b_mu_innerTrack_dz;   //!
   TBranch        *b_mu_innerTrack_dxyError;   //!
   TBranch        *b_mu_innerTrack_dzError;   //!
   TBranch        *b_mu_innerTrack_validFraction;   //!
   TBranch        *b_mu_bestTrack_dxy;   //!
   TBranch        *b_mu_bestTrack_dz;   //!
   TBranch        *b_mu_bestTrack_dxyError;   //!
   TBranch        *b_mu_bestTrack_dzError;   //!
   TBranch        *b_mu_numberOfMatches;   //!
   TBranch        *b_mu_numberOfValidMuonHits;   //!
   TBranch        *b_mu_pfIso03_sumChargedHadronPt;   //!
   TBranch        *b_mu_pfIso03_sumNeutralHadronEt;   //!
   TBranch        *b_mu_pfIso03_sumPhotonEt;   //!
   TBranch        *b_mu_pfIso03_sumPUPt;   //!
   TBranch        *b_mu_lepMVA;   //!
   TBranch        *b_mu_lepMVA_neuRelIso;
   TBranch        *b_mu_lepMVA_chRelIso;
   TBranch        *b_mu_lepMVA_jetDR;
   TBranch        *b_mu_lepMVA_jetPtRatio;
   TBranch        *b_mu_lepMVA_jetBTagCSV;
   TBranch        *b_mu_lepMVA_sip3d;
   TBranch        *b_mu_lepMVA_dxy;
   TBranch        *b_mu_lepMVA_dz;
   TBranch        *b_mu_lepMVA_mvaId;
   TBranch        *b_mu_innerTrack_pt;   //!
   TBranch        *b_mu_innerTrack_ptError;   //!
   TBranch        *b_mu_dB3D;   //!
   TBranch        *b_mu_edB3D;   //!
   TBranch        *b_jet_n;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_m;   //!
   TBranch        *b_jet_E;   //!
   TBranch        *b_jet_ntrk;   //!
   TBranch        *b_jet_CSV;   //!
   TBranch        *b_jet_CSVv2;   //!
   TBranch        *b_jet_partonFlavour;   //!
   TBranch        *b_jet_hadronFlavour;   //!
   TBranch        *b_jet_neutralHadronEnergy;   //!
   TBranch        *b_jet_neutralEmEnergy;   //!
   TBranch        *b_jet_chargedHadronEnergy;   //!
   TBranch        *b_jet_chargedEmEnergy;   //!
   TBranch        *b_jet_electronEnergy;   //!
   TBranch        *b_jet_muonEnergy;   //!
   TBranch        *b_jet_photonEnergy;   //!
   TBranch        *b_jet_genJet_pt;   //!
   TBranch        *b_jet_genJet_eta;   //!
   TBranch        *b_jet_genJet_phi;   //!
   TBranch        *b_jet_genJet_m;   //!
   TBranch        *b_jet_genJet_E;   //!
   TBranch        *b_jet_genJet_status;   //!
   TBranch        *b_jet_genJet_id;   //!
   TBranch        *b_jet_pileupJetId;   //!

   TBranch        *b_gen_n;
   TBranch        *b_gen_pt;
   TBranch        *b_gen_eta;
   TBranch        *b_gen_phi;
   TBranch        *b_gen_m;
   TBranch        *b_gen_id;
   TBranch        *b_gen_status;
   TBranch        *b_gen_index;
   TBranch        *b_gen_mother_index;
   
   TBranch        *b_mc_truth_h0_id;
   TBranch        *b_mc_truth_h0W1_id;
   TBranch        *b_mc_truth_h0Wl1_id;
   TBranch        *b_mc_truth_h0Wq11_id;
   TBranch        *b_mc_truth_h0Wq21_id;
   TBranch        *b_mc_truth_h0W2_id;
   TBranch        *b_mc_truth_h0Wl2_id;
   TBranch        *b_mc_truth_h0Wq12_id;
   TBranch        *b_mc_truth_h0Wq22_id;
   TBranch        *b_mc_truth_h0Z1_id;
   TBranch        *b_mc_truth_h0Zl11_id;
   TBranch        *b_mc_truth_h0Zl21_id;
   TBranch        *b_mc_truth_h0Zq11_id;
   TBranch        *b_mc_truth_h0Zq21_id;
   TBranch        *b_mc_truth_h0Z2_id;
   TBranch        *b_mc_truth_h0Zl12_id;
   TBranch        *b_mc_truth_h0Zl22_id;
   TBranch        *b_mc_truth_h0Zq12_id;
   TBranch        *b_mc_truth_h0Zq22_id;
   TBranch        *b_mc_truth_h0tau1_id;
   TBranch        *b_mc_truth_h0tau2_id;
   TBranch        *b_mc_truth_h0b1_id;
   TBranch        *b_mc_truth_h0b2_id;
   
   TBranch           *b_mc_truth_t1_id;
   TBranch           *b_mc_truth_t2_id;
   TBranch           *b_mc_truth_tb1_id;
   TBranch           *b_mc_truth_tb2_id;
   TBranch           *b_mc_truth_tW1_id;
   TBranch           *b_mc_truth_tW2_id;
   TBranch           *b_mc_truth_tWl1_id;
   TBranch           *b_mc_truth_tWl2_id;
   TBranch           *b_mc_truth_tWq11_id;
   TBranch           *b_mc_truth_tWq21_id;
   TBranch           *b_mc_truth_tWq12_id;
   TBranch           *b_mc_truth_tWq22_id;

   TBranch           *b_mc_truth_t_id;
   TBranch           *b_mc_truth_tb_id;
   TBranch           *b_mc_truth_tW_id;
   TBranch           *b_mc_truth_tWl_id;
   TBranch           *b_mc_truth_tWq1_id;
   TBranch           *b_mc_truth_tWq2_id;
   
   Tree(TChain *tree=0,char* fname="output.root",std::string treename="FlatTree/tree");
   virtual ~Tree();
   virtual Int_t    GetEntry(Long64_t entry);
//   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TChain *ch);
   virtual void     registerInputBranches(TChain *ch);
};

#endif
