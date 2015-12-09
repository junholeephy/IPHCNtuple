
#ifndef READGENFLATTREE_H
#define READGENFLATTREE_H


class ReadGenFlatTree {

  public:

  string InputFileName;
  TFile* fInput;
  TTree* tInput;
 
//  string OutputFileName;
  TFile* fOutput;
  TTree* tOutput;

  ReadGenFlatTree();
  ~ReadGenFlatTree();
  void InitializeDryRun(string);
  void InitializeMEMRun(string);
  void FillGenMultilepton(Long64_t, MultiLepton*);  
  int ApplyGenSelection(Long64_t, MultiLepton*);
  void WriteMultilepton(MultiLepton*);
  void ReadMultilepton(Long64_t, MultiLepton*);

  Float_t mc_weight;
  Int_t mc_truth_h0_id;
  TLorentzVector* mc_truth_h0_p4;
  Int_t mc_truth_h0Wl1_id;
  TLorentzVector* mc_truth_h0Wl1_p4;
  Int_t mc_truth_h0Wl2_id;
  TLorentzVector* mc_truth_h0Wl2_p4;
  TLorentzVector* mc_truth_h0Wnu1_p4;
  TLorentzVector* mc_truth_h0Wnu2_p4;
  Int_t mc_truth_h0Wq11_id;
  TLorentzVector* mc_truth_h0Wq11_p4;
  Int_t mc_truth_h0Wq12_id;
  TLorentzVector* mc_truth_h0Wq12_p4;
  Int_t mc_truth_h0Wq21_id;
  TLorentzVector* mc_truth_h0Wq21_p4;
  Int_t mc_truth_h0Wq22_id;
  TLorentzVector* mc_truth_h0Wq22_p4;
  Int_t mc_truth_t1_id;
  TLorentzVector* mc_truth_t1_p4;
  Int_t mc_truth_t2_id;
  TLorentzVector* mc_truth_t2_p4;
  Int_t mc_truth_tb1_id;
  TLorentzVector* mc_truth_tb1_p4;
  Int_t mc_truth_tb2_id;
  TLorentzVector* mc_truth_tb2_p4;
  Int_t mc_truth_tWl1_id;
  TLorentzVector* mc_truth_tWl1_p4;
  TLorentzVector* mc_truth_tWnu1_p4;
  Int_t mc_truth_tWl2_id;
  TLorentzVector* mc_truth_tWl2_p4;
  TLorentzVector* mc_truth_tWnu2_p4;
  Int_t mc_truth_tWq11_id;
  TLorentzVector* mc_truth_tWq11_p4;
  Int_t mc_truth_tWq21_id;
  TLorentzVector* mc_truth_tWq21_p4;
  Int_t mc_truth_tWq12_id;
  TLorentzVector* mc_truth_tWq12_p4;
  Int_t mc_truth_tWq22_id;
  TLorentzVector* mc_truth_tWq22_p4;
  Int_t mc_truth_Z_id;
  TLorentzVector* mc_truth_Z_p4;
  Int_t mc_truth_Zl1_id;
  TLorentzVector* mc_truth_Zl1_p4;
  Int_t mc_truth_Zl2_id;
  TLorentzVector* mc_truth_Zl2_p4;
  Int_t mc_truth_W_id;
  TLorentzVector*  mc_truth_W_p4;
  Int_t mc_truth_Wl_id;
  TLorentzVector*  mc_truth_Wl_p4;
  Int_t genJet_n;
  Int_t mc_truth_gammal1_id;
  TLorentzVector*  mc_truth_gammal1_p4;
  Int_t mc_truth_gammal2_id;
  TLorentzVector*  mc_truth_gammal2_p4;
  std::vector<float>* genJet_pt;
  std::vector<float>* genJet_eta;
  std::vector<float>* genJet_phi;
  std::vector<float>* genJet_E;

  TBranch* b_mc_event;
  TBranch* b_mc_weight;
  TBranch* b_mc_truth_h0_id;
  TBranch* b_mc_truth_h0_p4;
  TBranch* b_mc_truth_h0Wl1_id;
  TBranch* b_mc_truth_h0Wl1_p4;
  TBranch* b_mc_truth_h0Wl2_id;
  TBranch* b_mc_truth_h0Wl2_p4;
  TBranch* b_mc_truth_h0Wnu1_p4;
  TBranch* b_mc_truth_h0Wnu2_p4;
  TBranch* b_mc_truth_h0Wq11_id;
  TBranch* b_mc_truth_h0Wq11_p4;
  TBranch* b_mc_truth_h0Wq12_id;
  TBranch* b_mc_truth_h0Wq12_p4;
  TBranch* b_mc_truth_h0Wq21_id;
  TBranch* b_mc_truth_h0Wq21_p4;
  TBranch* b_mc_truth_h0Wq22_id;
  TBranch* b_mc_truth_h0Wq22_p4;
  TBranch* b_mc_truth_t1_id;
  TBranch* b_mc_truth_t1_p4;
  TBranch* b_mc_truth_t2_id;
  TBranch* b_mc_truth_t2_p4;
  TBranch* b_mc_truth_tb1_id;
  TBranch* b_mc_truth_tb1_p4;
  TBranch* b_mc_truth_tb2_id;
  TBranch* b_mc_truth_tb2_p4;
  TBranch* b_mc_truth_tWl1_id;
  TBranch* b_mc_truth_tWl1_p4;
  TBranch* b_mc_truth_tWl2_id;
  TBranch* b_mc_truth_tWl2_p4;
  TBranch* b_mc_truth_tWnu1_p4;
  TBranch* b_mc_truth_tWnu2_p4;
  TBranch* b_mc_truth_tWq11_id;
  TBranch* b_mc_truth_tWq11_p4;
  TBranch* b_mc_truth_tWq21_id;
  TBranch* b_mc_truth_tWq21_p4;
  TBranch* b_mc_truth_tWq12_id;
  TBranch* b_mc_truth_tWq12_p4;
  TBranch* b_mc_truth_tWq22_id;
  TBranch* b_mc_truth_tWq22_p4;
  TBranch* b_mc_truth_Z_id;
  TBranch* b_mc_truth_Z_p4;
  TBranch* b_mc_truth_Zl1_id;
  TBranch* b_mc_truth_Zl1_p4;
  TBranch* b_mc_truth_Zl2_id;
  TBranch* b_mc_truth_Zl2_p4;
  TBranch* b_mc_truth_W_id;
  TBranch* b_mc_truth_W_p4;
  TBranch* b_mc_truth_Wl_id;
  TBranch* b_mc_truth_Wl_p4;
  TBranch* b_mc_truth_gammal1_id;
  TBranch* b_mc_truth_gammal1_p4;
  TBranch* b_mc_truth_gammal2_id;
  TBranch* b_mc_truth_gammal2_p4;
  TBranch* b_genJet_n;
  TBranch* b_genJet_pt;
  TBranch* b_genJet_eta;
  TBranch* b_genJet_phi;
  TBranch* b_genJet_E;

  Int_t mc_event;
  Float_t mc_totp4_px;
  Float_t mc_totp4_py;
  Float_t mc_totp4_pt;
  Float_t mc_thad_pt;
  Float_t mc_thad_b_pt;
  Float_t mc_thad_b_eta;
  Float_t mc_thad_j1_pt;
  Float_t mc_thad_j1_eta;
  Float_t mc_thad_j2_pt;
  Float_t mc_thad_j2_eta;
  Float_t mc_tlep_pt;
  Float_t mc_tlep_b_pt;
  Float_t mc_tlep_b_eta;
  Float_t mc_tlep_l_pt;
  Float_t mc_tlep_l_eta;
  Float_t mc_tlep2_pt;
  Float_t mc_tlep2_b_pt;
  Float_t mc_tlep2_b_eta;
  Float_t mc_tlep2_l_pt;
  Float_t mc_tlep2_l_eta;
  Int_t   mc_ttbar_decay;
  Int_t   mc_boson_decay;
  Float_t mc_boson_pt;
  Float_t mc_boson_l1_pt;
  Float_t mc_boson_l1_eta;
  Float_t mc_boson_l2_pt;
  Float_t mc_boson_l2_eta;
  Float_t mc_boson_ll_mass;
  Float_t mc_boson_ll_pt;
  Float_t mc_boson_ll_dphi;
  Float_t mc_boson_j1_pt;
  Float_t mc_boson_j1_eta;
  Float_t mc_boson_j2_pt;
  Float_t mc_boson_j2_eta;
  Float_t mc_boson_jj_mass;
  Float_t mc_boson_jj_pt;
  Float_t mc_boson_jj_dphi;
  Float_t mc_met;
  Int_t mc_njets25;

  Int_t mc_3l_category;

  Int_t mc_ttZhypAllowed;
  Int_t mc_hasLLcombZpeak;
  Int_t mc_passMllGt12;
  Int_t mc_passLepPresel;
  Int_t mc_passJetPresel25;
  Int_t mc_passBjetPresel25;

  Int_t multilepton_Bjet1_Id;
  TLorentzVector multilepton_Bjet1_P4;
  Int_t multilepton_Bjet2_Id;
  TLorentzVector multilepton_Bjet2_P4;
  Int_t multilepton_Lepton1_Id;
  TLorentzVector multilepton_Lepton1_P4;
  Int_t multilepton_Lepton2_Id;
  TLorentzVector multilepton_Lepton2_P4;
  Int_t multilepton_Lepton3_Id;
  TLorentzVector multilepton_Lepton3_P4;
  Int_t multilepton_JetHighestPt1_Id;
  TLorentzVector multilepton_JetHighestPt1_P4;
  Int_t multilepton_JetHighestPt2_Id;
  TLorentzVector multilepton_JetHighestPt2_P4;
  Float_t multilepton_JetHighestPt_Mjj;
  Int_t multilepton_JetClosestMw1_Id;
  TLorentzVector multilepton_JetClosestMw1_P4;
  Int_t multilepton_JetClosestMw2_Id;
  TLorentzVector multilepton_JetClosestMw2_P4;
  Float_t multilepton_JetClosestMw_Mjj;
  Int_t multilepton_JetLowestMjj1_Id;
  TLorentzVector multilepton_JetLowestMjj1_P4;
  Int_t multilepton_JetLowestMjj2_Id;
  TLorentzVector multilepton_JetLowestMjj2_P4;
  Float_t multilepton_JetLowestMjj_Mjj;
  TLorentzVector multilepton_mET;
  TLorentzVector multilepton_Ptot;

  TLorentzVector* multilepton_Bjet1_P4_ptr;
  TLorentzVector* multilepton_Bjet2_P4_ptr;
  TLorentzVector* multilepton_Lepton1_P4_ptr;
  TLorentzVector* multilepton_Lepton2_P4_ptr;
  TLorentzVector* multilepton_Lepton3_P4_ptr;
  TLorentzVector* multilepton_JetHighestPt1_P4_ptr;
  TLorentzVector* multilepton_JetHighestPt2_P4_ptr;
  TLorentzVector* multilepton_JetClosestMw1_P4_ptr;
  TLorentzVector* multilepton_JetClosestMw2_P4_ptr;
  TLorentzVector* multilepton_JetLowestMjj1_P4_ptr;
  TLorentzVector* multilepton_JetLowestMjj2_P4_ptr;

  TLorentzVector* multilepton_mET_ptr;
  TLorentzVector* multilepton_Ptot_ptr;

  Double_t mc_mem_tthfl_weight;
  Double_t mc_mem_tthfl_weight_log;
  Double_t mc_mem_tthfl_weight_err;
  Float_t mc_mem_tthfl_weight_chi2;
  Float_t mc_mem_tthfl_weight_time;
  Double_t mc_mem_tthfl_weight_max;
  Double_t mc_mem_tthfl_weight_avg;

  Double_t mc_mem_tthsl_weight;
  Double_t mc_mem_tthsl_weight_log;
  Double_t mc_mem_tthsl_weight_err;
  Float_t mc_mem_tthsl_weight_chi2;
  Float_t mc_mem_tthsl_weight_time;
  Double_t mc_mem_tthsl_weight_max;
  Double_t mc_mem_tthsl_weight_avg;

  Double_t mc_mem_tth_weight;
  Double_t mc_mem_tth_weight_log;
  Double_t mc_mem_tth_weight_err;
  Float_t mc_mem_tth_weight_chi2;
  Float_t mc_mem_tth_weight_time;
  Double_t mc_mem_tth_weight_max;
  Double_t mc_mem_tth_weight_avg;

  Double_t mc_mem_ttz_weight;
  Double_t mc_mem_ttz_weight_log;
  Double_t mc_mem_ttz_weight_err;
  Float_t mc_mem_ttz_weight_chi2;
  Float_t mc_mem_ttz_weight_time;
  Double_t mc_mem_ttz_weight_max;
  Double_t mc_mem_ttz_weight_avg;

  Double_t mc_mem_ttw_weight;
  Double_t mc_mem_ttw_weight_log;
  Double_t mc_mem_ttw_weight_err;
  Float_t mc_mem_ttw_weight_chi2;
  Float_t mc_mem_ttw_weight_time;
  Double_t mc_mem_ttw_weight_max;
  Double_t mc_mem_ttw_weight_avg;

  Double_t mc_mem_ttwjj_weight;
  Double_t mc_mem_ttwjj_weight_log;
  Double_t mc_mem_ttwjj_weight_err;
  Float_t mc_mem_ttwjj_weight_chi2;
  Float_t mc_mem_ttwjj_weight_time;
  Double_t mc_mem_ttwjj_weight_max;
  Double_t mc_mem_ttwjj_weight_avg;

  //Double_t mc_mem_ttz_tthfl_likelihood;
  //Double_t mc_mem_ttz_tthsl_likelihood;  
  //Double_t mc_mem_ttw_tthfl_likelihood;
  //Double_t mc_mem_ttw_tthsl_likelihood;
  Double_t mc_mem_ttz_tth_likelihood;
  Double_t mc_mem_ttz_tth_likelihood_nlog;
  Double_t mc_mem_ttz_tth_likelihood_max;
  Double_t mc_mem_ttz_tth_likelihood_avg;

  Double_t mc_mem_ttw_tth_likelihood;
  Double_t mc_mem_ttw_tth_likelihood_nlog;
  Double_t mc_mem_ttw_tth_likelihood_max;
  Double_t mc_mem_ttw_tth_likelihood_avg;

  Double_t mc_mem_ttwjj_tth_likelihood;
  Double_t mc_mem_ttwjj_tth_likelihood_nlog;
  Double_t mc_mem_ttwjj_tth_likelihood_max;
  Double_t mc_mem_ttwjj_tth_likelihood_avg;

  Double_t mc_mem_ttv_tth_likelihood;
  Double_t mc_mem_ttv_tth_likelihood_nlog;
  Double_t mc_mem_ttv_tth_likelihood_max;
  Double_t mc_mem_ttv_tth_likelihood_avg;

  Double_t mc_mem_ttvjj_tth_likelihood;
  Double_t mc_mem_ttvjj_tth_likelihood_nlog;
  Double_t mc_mem_ttvjj_tth_likelihood_max;
  Double_t mc_mem_ttvjj_tth_likelihood_avg;

  TBranch* b_mc_3l_category;
  TBranch* b_mc_ttbar_decay;
  TBranch* b_mc_boson_decay;
  TBranch* b_mc_ttZhypAllowed;
  TBranch* b_multilepton_Bjet1_Id;
  TBranch* b_multilepton_Bjet1_P4;
  TBranch* b_multilepton_Bjet2_Id;
  TBranch* b_multilepton_Bjet2_P4;
  TBranch* b_multilepton_Lepton1_Id;
  TBranch* b_multilepton_Lepton1_P4;
  TBranch* b_multilepton_Lepton2_Id;
  TBranch* b_multilepton_Lepton2_P4;
  TBranch* b_multilepton_Lepton3_Id;
  TBranch* b_multilepton_Lepton3_P4;
  TBranch* b_multilepton_JetHighestPt1_Id;
  TBranch* b_multilepton_JetHighestPt1_P4;
  TBranch* b_multilepton_JetHighestPt2_Id;
  TBranch* b_multilepton_JetHighestPt2_P4;
  TBranch* b_multilepton_JetClosestMw1_Id;
  TBranch* b_multilepton_JetClosestMw1_P4;
  TBranch* b_multilepton_JetClosestMw2_Id;
  TBranch* b_multilepton_JetClosestMw2_P4;
  TBranch* b_multilepton_JetLowestMjj1_Id;
  TBranch* b_multilepton_JetLowestMjj1_P4;
  TBranch* b_multilepton_JetLowestMjj2_Id;
  TBranch* b_multilepton_JetLowestMjj2_P4;
  TBranch* b_multilepton_JetHighestPt_Mjj;
  TBranch* b_multilepton_JetClosestMw_Mjj;
  TBranch* b_multilepton_JetLowestMjj_Mjj;
  TBranch* b_multilepton_mET;
  TBranch* b_multilepton_Ptot;


  TH1F** hMEPhaseSpace_Error;
  TH1F** hMEPhaseSpace_ErrorTot;
  TH1F** hMEPhaseSpace_ErrorTot_Fail;
  TH1F** hMEPhaseSpace_ErrorTot_Pass;

  private:
};

ReadGenFlatTree::ReadGenFlatTree(){

}

ReadGenFlatTree::~ReadGenFlatTree(){

}

void ReadGenFlatTree::InitializeDryRun(string InputFileName){

  cout << "Opening input tree"<<endl;

  fInput = TFile::Open(InputFileName.c_str(),"READ");
  //tInput = (TTree*)fInput->Get("FlatTree/tree");
  tInput = (TTree*)fInput->Get("tree");

  mc_truth_h0_p4 = 0;
  mc_truth_h0Wl1_p4 = 0;
  mc_truth_h0Wl2_p4 = 0;
  mc_truth_h0Wnu1_p4 = 0;
  mc_truth_h0Wnu2_p4 = 0;
  mc_truth_h0Wq11_p4 = 0;
  mc_truth_h0Wq12_p4 = 0;
  mc_truth_h0Wq21_p4 = 0;
  mc_truth_h0Wq22_p4 = 0;
  mc_truth_t1_p4 = 0;
  mc_truth_t2_p4 = 0;
  mc_truth_tb1_p4 = 0;
  mc_truth_tb2_p4 = 0;
  mc_truth_tWl1_p4 = 0;
  mc_truth_tWnu1_p4 = 0;
  mc_truth_tWl2_p4 = 0;
  mc_truth_tWnu2_p4 = 0;
  mc_truth_tWq11_p4 = 0;
  mc_truth_tWq21_p4 = 0;
  mc_truth_tWq12_p4 = 0;
  mc_truth_tWq22_p4 = 0;
  mc_truth_Z_p4 = 0;
  mc_truth_Zl1_p4 = 0;
  mc_truth_Zl2_p4 = 0;
  mc_truth_W_p4 = 0;
  mc_truth_Wl_p4 = 0;
  mc_truth_gammal1_p4 = 0;
  mc_truth_gammal2_p4 = 0;
  genJet_pt = 0;
  genJet_eta = 0;
  genJet_phi = 0;
  genJet_E = 0;


  tInput->SetBranchAddress("mc_weight",&mc_weight,&b_mc_weight);
  tInput->SetBranchAddress("mc_truth_h0_id",&mc_truth_h0_id,&b_mc_truth_h0_id);
  tInput->SetBranchAddress("mc_truth_h0_p4",&mc_truth_h0_p4,&b_mc_truth_h0_p4);
  tInput->SetBranchAddress("mc_truth_h0Wl1_id",&mc_truth_h0Wl1_id,&b_mc_truth_h0Wl1_id);
  tInput->SetBranchAddress("mc_truth_h0Wl1_p4",&mc_truth_h0Wl1_p4,&b_mc_truth_h0Wl1_p4);
  tInput->SetBranchAddress("mc_truth_h0Wl2_id",&mc_truth_h0Wl2_id,&b_mc_truth_h0Wl2_id);
  tInput->SetBranchAddress("mc_truth_h0Wl2_p4",&mc_truth_h0Wl2_p4,&b_mc_truth_h0Wl2_p4);
  tInput->SetBranchAddress("mc_truth_h0Wnu1_p4",&mc_truth_h0Wnu1_p4,&b_mc_truth_h0Wnu1_p4);
  tInput->SetBranchAddress("mc_truth_h0Wnu2_p4",&mc_truth_h0Wnu2_p4,&b_mc_truth_h0Wnu2_p4);
  tInput->SetBranchAddress("mc_truth_h0Wq11_id",&mc_truth_h0Wq11_id,&b_mc_truth_h0Wq11_id);
  tInput->SetBranchAddress("mc_truth_h0Wq11_p4",&mc_truth_h0Wq11_p4,&b_mc_truth_h0Wq11_p4);
  tInput->SetBranchAddress("mc_truth_h0Wq12_id",&mc_truth_h0Wq12_id,&b_mc_truth_h0Wq12_id);
  tInput->SetBranchAddress("mc_truth_h0Wq12_p4",&mc_truth_h0Wq12_p4,&b_mc_truth_h0Wq12_p4);
  tInput->SetBranchAddress("mc_truth_h0Wq21_id",&mc_truth_h0Wq21_id,&b_mc_truth_h0Wq21_id);
  tInput->SetBranchAddress("mc_truth_h0Wq21_p4",&mc_truth_h0Wq21_p4,&b_mc_truth_h0Wq21_p4);
  tInput->SetBranchAddress("mc_truth_h0Wq22_id",&mc_truth_h0Wq22_id,&b_mc_truth_h0Wq22_id);
  tInput->SetBranchAddress("mc_truth_h0Wq22_p4",&mc_truth_h0Wq22_p4,&b_mc_truth_h0Wq22_p4);
  tInput->SetBranchAddress("mc_truth_t1_id",&mc_truth_t1_id,&b_mc_truth_t1_id);
  tInput->SetBranchAddress("mc_truth_t1_p4",&mc_truth_t1_p4,&b_mc_truth_t1_p4);
  tInput->SetBranchAddress("mc_truth_t2_id",&mc_truth_t2_id,&b_mc_truth_t2_id);
  tInput->SetBranchAddress("mc_truth_t2_p4",&mc_truth_t2_p4,&b_mc_truth_t2_p4);
  tInput->SetBranchAddress("mc_truth_tb1_id",&mc_truth_tb1_id,&b_mc_truth_tb1_id);
  tInput->SetBranchAddress("mc_truth_tb1_p4",&mc_truth_tb1_p4,&b_mc_truth_tb1_p4);
  tInput->SetBranchAddress("mc_truth_tb2_id",&mc_truth_tb2_id,&b_mc_truth_tb2_id);
  tInput->SetBranchAddress("mc_truth_tb2_p4",&mc_truth_tb2_p4,&b_mc_truth_tb2_p4);
  tInput->SetBranchAddress("mc_truth_tWl1_id",&mc_truth_tWl1_id,&b_mc_truth_tWl1_id);
  tInput->SetBranchAddress("mc_truth_tWl1_p4",&mc_truth_tWl1_p4,&b_mc_truth_tWl1_p4);
  tInput->SetBranchAddress("mc_truth_tWl2_id",&mc_truth_tWl2_id,&b_mc_truth_tWl2_id);
  tInput->SetBranchAddress("mc_truth_tWl2_p4",&mc_truth_tWl2_p4,&b_mc_truth_tWl2_p4);
  tInput->SetBranchAddress("mc_truth_tWnu1_p4",&mc_truth_tWnu1_p4,&b_mc_truth_tWnu1_p4);
  tInput->SetBranchAddress("mc_truth_tWnu2_p4",&mc_truth_tWnu2_p4,&b_mc_truth_tWnu2_p4);
  tInput->SetBranchAddress("mc_truth_tWq11_id",&mc_truth_tWq11_id,&b_mc_truth_tWq11_id);
  tInput->SetBranchAddress("mc_truth_tWq11_p4",&mc_truth_tWq11_p4,&b_mc_truth_tWq11_p4);
  tInput->SetBranchAddress("mc_truth_tWq21_id",&mc_truth_tWq21_id,&b_mc_truth_tWq21_id);
  tInput->SetBranchAddress("mc_truth_tWq21_p4",&mc_truth_tWq21_p4,&b_mc_truth_tWq21_p4);
  tInput->SetBranchAddress("mc_truth_tWq12_id",&mc_truth_tWq12_id,&b_mc_truth_tWq12_id);
  tInput->SetBranchAddress("mc_truth_tWq12_p4",&mc_truth_tWq12_p4,&b_mc_truth_tWq12_p4);
  tInput->SetBranchAddress("mc_truth_tWq22_id",&mc_truth_tWq22_id,&b_mc_truth_tWq22_id);
  tInput->SetBranchAddress("mc_truth_tWq22_p4",&mc_truth_tWq22_p4,&b_mc_truth_tWq22_p4);
  tInput->SetBranchAddress("mc_truth_Z_id",&mc_truth_Z_id,&b_mc_truth_Z_id);
  tInput->SetBranchAddress("mc_truth_Z_p4",&mc_truth_Z_p4,&b_mc_truth_Z_p4);
  tInput->SetBranchAddress("mc_truth_Zl1_id",&mc_truth_Zl1_id,&b_mc_truth_Zl1_id);
  tInput->SetBranchAddress("mc_truth_Zl1_p4",&mc_truth_Zl1_p4,&b_mc_truth_Zl1_p4);
  tInput->SetBranchAddress("mc_truth_Zl2_id",&mc_truth_Zl2_id,&b_mc_truth_Zl2_id);
  tInput->SetBranchAddress("mc_truth_Zl2_p4",&mc_truth_Zl2_p4,&b_mc_truth_Zl2_p4);
  tInput->SetBranchAddress("mc_truth_W_id",&mc_truth_W_id,&b_mc_truth_W_id);
  tInput->SetBranchAddress("mc_truth_W_p4",&mc_truth_W_p4,&b_mc_truth_W_p4);
  tInput->SetBranchAddress("mc_truth_Wl_id",&mc_truth_Wl_id,&b_mc_truth_Wl_id);
  tInput->SetBranchAddress("mc_truth_Wl_p4",&mc_truth_Wl_p4,&b_mc_truth_Wl_p4);
  tInput->SetBranchAddress("mc_truth_gammal1_id",&mc_truth_gammal1_id,&b_mc_truth_gammal1_id);
  tInput->SetBranchAddress("mc_truth_gammal1_p4",&mc_truth_gammal1_p4,&b_mc_truth_gammal1_p4);
  tInput->SetBranchAddress("mc_truth_gammal2_id",&mc_truth_gammal2_id,&b_mc_truth_gammal2_id);
  tInput->SetBranchAddress("mc_truth_gammal2_p4",&mc_truth_gammal2_p4,&b_mc_truth_gammal2_p4);
  tInput->SetBranchAddress("genJet_n",&genJet_n,&b_genJet_n);
  tInput->SetBranchAddress("genJet_pt",&genJet_pt,&b_genJet_pt);
  tInput->SetBranchAddress("genJet_eta",&genJet_eta,&b_genJet_eta);
  tInput->SetBranchAddress("genJet_phi",&genJet_phi,&b_genJet_phi);
  tInput->SetBranchAddress("genJet_E",&genJet_E,&b_genJet_E);

  cout << "Creating output tree"<<endl;

  fOutput = new TFile("output.root", "RECREATE");
  tOutput = new TTree("Tree", "Tree");

  tOutput->Branch("mc_event",&mc_event,"mc_event/I");
  tOutput->Branch("mc_weight",&mc_weight,"mc_weight/F");
  tOutput->Branch("mc_totp4_px",&mc_totp4_px,"mc_totp4_px/F");
  tOutput->Branch("mc_totp4_py",&mc_totp4_py,"mc_totp4_py/F");
  tOutput->Branch("mc_totp4_pt",&mc_totp4_pt,"mc_totp4_pt/F");
  tOutput->Branch("mc_thad_pt",&mc_thad_pt,"mc_thad_pt/F");
  tOutput->Branch("mc_thad_b_pt",&mc_thad_b_pt,"mc_thad_b_pt/F");
  tOutput->Branch("mc_thad_b_eta",&mc_thad_b_eta,"mc_thad_b_eta/F");
  tOutput->Branch("mc_thad_j1_pt",&mc_thad_j1_pt,"mc_thad_j1_pt/F");
  tOutput->Branch("mc_thad_j1_eta",&mc_thad_j1_eta,"mc_thad_j1_eta/F");
  tOutput->Branch("mc_thad_j2_pt",&mc_thad_j2_pt,"mc_thad_j2_pt/F");
  tOutput->Branch("mc_thad_j2_eta",&mc_thad_j2_eta,"mc_thad_j2_eta/F");
  tOutput->Branch("mc_tlep_pt",&mc_tlep_pt,"mc_tlep_pt/F");
  tOutput->Branch("mc_tlep_b_pt",&mc_tlep_b_pt,"mc_tlep_b_pt/F");
  tOutput->Branch("mc_tlep_b_eta",&mc_tlep_b_eta,"mc_tlep_b_eta/F");
  tOutput->Branch("mc_tlep_l_pt",&mc_tlep_l_pt,"mc_tlep_l_pt/F");
  tOutput->Branch("mc_tlep_l_eta",&mc_tlep_l_eta,"mc_tlep_l_eta/F");
  tOutput->Branch("mc_tlep2_pt",&mc_tlep2_pt,"mc_tlep2_pt/F");
  tOutput->Branch("mc_tlep2_b_pt",&mc_tlep2_b_pt,"mc_tlep2_b_pt/F");
  tOutput->Branch("mc_tlep2_b_eta",&mc_tlep2_b_eta,"mc_tlep2_b_eta/F");
  tOutput->Branch("mc_tlep2_l_pt",&mc_tlep2_l_pt,"mc_tlep2_l_pt/F");
  tOutput->Branch("mc_tlep2_l_eta",&mc_tlep2_l_eta,"mc_tlep2_l_eta/F");
  tOutput->Branch("mc_ttbar_decay",&mc_ttbar_decay,"mc_ttbar_decay/I");
  tOutput->Branch("mc_boson_decay",&mc_boson_decay,"mc_boson_decay/I");
  tOutput->Branch("mc_boson_pt",&mc_boson_pt,"mc_boson_pt/F");
  tOutput->Branch("mc_boson_l1_pt",&mc_boson_l1_pt,"mc_boson_l1_pt/F");
  tOutput->Branch("mc_boson_l1_eta",&mc_boson_l1_eta,"mc_boson_l1_eta/F");
  tOutput->Branch("mc_boson_l2_pt",&mc_boson_l2_pt,"mc_boson_l2_pt/F");
  tOutput->Branch("mc_boson_l2_eta",&mc_boson_l2_eta,"mc_boson_l2_eta/F");
  tOutput->Branch("mc_boson_ll_mass",&mc_boson_ll_mass,"mc_boson_ll_mass/F");
  tOutput->Branch("mc_boson_ll_pt",&mc_boson_ll_pt,"mc_boson_ll_pt/F");
  tOutput->Branch("mc_boson_ll_dphi",&mc_boson_ll_dphi,"mc_boson_ll_dphi/F");
  tOutput->Branch("mc_boson_j1_pt",&mc_boson_j1_pt,"mc_boson_j1_pt/F");
  tOutput->Branch("mc_boson_j1_eta",&mc_boson_j1_eta,"mc_boson_j1_eta/F");
  tOutput->Branch("mc_boson_j2_pt",&mc_boson_j2_pt,"mc_boson_j2_pt/F");
  tOutput->Branch("mc_boson_j2_eta",&mc_boson_j2_eta,"mc_boson_j2_eta/F");
  tOutput->Branch("mc_boson_jj_mass",&mc_boson_jj_mass,"mc_boson_jj_mass/F");
  tOutput->Branch("mc_boson_jj_pt",&mc_boson_jj_pt,"mc_boson_jj_pt/F");
  tOutput->Branch("mc_boson_jj_dphi",&mc_boson_jj_dphi,"mc_boson_jj_dphi/F");
  tOutput->Branch("mc_met",&mc_met,"mc_met/F");
  tOutput->Branch("mc_njets25",&mc_njets25,"mc_njets25/I");
  tOutput->Branch("mc_3l_category",&mc_3l_category,"mc_3l_category/I");
  tOutput->Branch("mc_ttZhypAllowed",&mc_ttZhypAllowed,"mc_ttZhypAllowed/I");
  tOutput->Branch("mc_hasLLcombZpeak",&mc_hasLLcombZpeak,"mc_hasLLcombZpeak/I");
  tOutput->Branch("mc_passMllGt12",&mc_passMllGt12,"mc_passMllGt12/I");
  tOutput->Branch("mc_passLepPresel",&mc_passLepPresel,"mc_passLepPresel/I");
  tOutput->Branch("mc_passJetPresel25",&mc_passJetPresel25,"mc_passJetPresel25/I");
  tOutput->Branch("mc_passBjetPresel25",&mc_passBjetPresel25,"mc_passBjetPresel25/I");

  cout << "multilepton variables"<<endl;
  tOutput->Branch("multilepton_Bjet1_Id",&multilepton_Bjet1_Id,"multilepton_Bjet1_Id/I");
  tOutput->Branch("multilepton_Bjet1_P4","TLorentzVector",&multilepton_Bjet1_P4);
  tOutput->Branch("multilepton_Bjet2_Id",&multilepton_Bjet2_Id,"multilepton_Bjet2_Id/I");
  tOutput->Branch("multilepton_Bjet2_P4","TLorentzVector",&multilepton_Bjet2_P4);
  tOutput->Branch("multilepton_Lepton1_Id",&multilepton_Lepton1_Id,"multilepton_Lepton1_Id/I");
  tOutput->Branch("multilepton_Lepton1_P4","TLorentzVector",&multilepton_Lepton1_P4);
  tOutput->Branch("multilepton_Lepton2_Id",&multilepton_Lepton2_Id,"multilepton_Lepton2_Id/I");
  tOutput->Branch("multilepton_Lepton2_P4","TLorentzVector",&multilepton_Lepton2_P4);
  tOutput->Branch("multilepton_Lepton3_Id",&multilepton_Lepton3_Id,"multilepton_Lepton3_Id/I");
  tOutput->Branch("multilepton_Lepton3_P4","TLorentzVector",&multilepton_Lepton3_P4);
  tOutput->Branch("multilepton_JetHighestPt1_Id",&multilepton_JetHighestPt1_Id,"multilepton_JetHighestPt1_Id/I");
  tOutput->Branch("multilepton_JetHighestPt1_P4","TLorentzVector",&multilepton_JetHighestPt1_P4);
  tOutput->Branch("multilepton_JetHighestPt2_Id",&multilepton_JetHighestPt2_Id,"multilepton_JetHighestPt2_Id/I");
  tOutput->Branch("multilepton_JetHighestPt2_P4","TLorentzVector",&multilepton_JetHighestPt2_P4);
  tOutput->Branch("multilepton_JetHighestPt_Mjj",&multilepton_JetHighestPt_Mjj,"multilepton_JetHighestPt_Mjj/F");
  tOutput->Branch("multilepton_JetClosestMw1_Id",&multilepton_JetClosestMw1_Id,"multilepton_JetClosestMw1_Id/I");
  tOutput->Branch("multilepton_JetClosestMw1_P4","TLorentzVector",&multilepton_JetClosestMw1_P4);
  tOutput->Branch("multilepton_JetClosestMw2_Id",&multilepton_JetClosestMw2_Id,"multilepton_JetClosestMw2_Id/I");
  tOutput->Branch("multilepton_JetClosestMw2_P4","TLorentzVector",&multilepton_JetClosestMw2_P4);
  tOutput->Branch("multilepton_JetClosestMw_Mjj",&multilepton_JetClosestMw_Mjj,"multilepton_JetClosestMw_Mjj/F");
  tOutput->Branch("multilepton_JetLowestMjj1_Id",&multilepton_JetLowestMjj1_Id,"multilepton_JetLowestMjj1_Id/I");
  tOutput->Branch("multilepton_JetLowestMjj1_P4","TLorentzVector",&multilepton_JetLowestMjj1_P4);
  tOutput->Branch("multilepton_JetLowestMjj2_Id",&multilepton_JetLowestMjj2_Id,"multilepton_JetLowestMjj2_Id/I");
  tOutput->Branch("multilepton_JetLowestMjj2_P4","TLorentzVector",&multilepton_JetLowestMjj2_P4);
  tOutput->Branch("multilepton_JetLowestMjj_Mjj",&multilepton_JetLowestMjj_Mjj,"multilepton_JetLowestMjj_Mjj/F");
  tOutput->Branch("multilepton_mET","TLorentzVector",&multilepton_mET);
  tOutput->Branch("multilepton_Ptot","TLorentzVector",&multilepton_Ptot);

  cout << "Tree created"<<endl;

}

void ReadGenFlatTree::InitializeMEMRun(string InputFileName){

  cout << "Opening input tree"<<endl;

  fInput = TFile::Open(InputFileName.c_str(),"READ");
  tInput = (TTree*)fInput->Get("Tree");

   multilepton_Bjet1_P4_ptr = 0;
   multilepton_Bjet2_P4_ptr = 0;
   multilepton_Lepton1_P4_ptr = 0;
   multilepton_Lepton2_P4_ptr = 0;
   multilepton_Lepton3_P4_ptr = 0;
   multilepton_JetHighestPt1_P4_ptr = 0;
   multilepton_JetHighestPt2_P4_ptr = 0;
   multilepton_JetClosestMw1_P4_ptr = 0;
   multilepton_JetClosestMw2_P4_ptr = 0;
   multilepton_JetLowestMjj1_P4_ptr = 0;
   multilepton_JetLowestMjj2_P4_ptr = 0;
   multilepton_mET_ptr = 0;
   multilepton_Ptot_ptr = 0;

  tInput->SetBranchAddress("mc_event",&mc_event,&b_mc_event);
  tInput->SetBranchAddress("mc_weight",&mc_weight,&b_mc_weight);
  tInput->SetBranchAddress("mc_3l_category",&mc_3l_category,&b_mc_3l_category);
  tInput->SetBranchAddress("mc_ttbar_decay",&mc_ttbar_decay,&b_mc_ttbar_decay);
  tInput->SetBranchAddress("mc_boson_decay",&mc_boson_decay,&b_mc_boson_decay);
  tInput->SetBranchAddress("mc_ttZhypAllowed",&mc_ttZhypAllowed,&b_mc_ttZhypAllowed);

  tInput->SetBranchAddress("multilepton_Bjet1_Id",&multilepton_Bjet1_Id,&b_multilepton_Bjet1_Id);
  tInput->SetBranchAddress("multilepton_Bjet1_P4",&multilepton_Bjet1_P4_ptr,&b_multilepton_Bjet1_P4);
  tInput->SetBranchAddress("multilepton_Bjet2_Id",&multilepton_Bjet2_Id,&b_multilepton_Bjet2_Id);
  tInput->SetBranchAddress("multilepton_Bjet2_P4",&multilepton_Bjet2_P4_ptr,&b_multilepton_Bjet2_P4);
  tInput->SetBranchAddress("multilepton_Lepton1_Id",&multilepton_Lepton1_Id,&b_multilepton_Lepton1_Id);
  tInput->SetBranchAddress("multilepton_Lepton1_P4",&multilepton_Lepton1_P4_ptr,&b_multilepton_Lepton1_P4);
  tInput->SetBranchAddress("multilepton_Lepton2_Id",&multilepton_Lepton2_Id,&b_multilepton_Lepton2_Id);
  tInput->SetBranchAddress("multilepton_Lepton2_P4",&multilepton_Lepton2_P4_ptr,&b_multilepton_Lepton2_P4);
  tInput->SetBranchAddress("multilepton_Lepton3_Id",&multilepton_Lepton3_Id,&b_multilepton_Lepton3_Id);
  tInput->SetBranchAddress("multilepton_Lepton3_P4",&multilepton_Lepton3_P4_ptr,&b_multilepton_Lepton3_P4);
  tInput->SetBranchAddress("multilepton_JetHighestPt1_Id",&multilepton_JetHighestPt1_Id,&b_multilepton_JetHighestPt1_Id);
  tInput->SetBranchAddress("multilepton_JetHighestPt1_P4",&multilepton_JetHighestPt1_P4_ptr,&b_multilepton_JetHighestPt1_P4);
  tInput->SetBranchAddress("multilepton_JetHighestPt2_Id",&multilepton_JetHighestPt2_Id,&b_multilepton_JetHighestPt2_Id);
  tInput->SetBranchAddress("multilepton_JetHighestPt2_P4",&multilepton_JetHighestPt2_P4_ptr,&b_multilepton_JetHighestPt2_P4);
  tInput->SetBranchAddress("multilepton_JetClosestMw1_Id",&multilepton_JetClosestMw1_Id,&b_multilepton_JetClosestMw1_Id);
  tInput->SetBranchAddress("multilepton_JetClosestMw1_P4",&multilepton_JetClosestMw1_P4_ptr,&b_multilepton_JetClosestMw1_P4);
  tInput->SetBranchAddress("multilepton_JetClosestMw2_Id",&multilepton_JetClosestMw2_Id,&b_multilepton_JetClosestMw2_Id);
  tInput->SetBranchAddress("multilepton_JetClosestMw2_P4",&multilepton_JetClosestMw2_P4_ptr,&b_multilepton_JetClosestMw2_P4);
  tInput->SetBranchAddress("multilepton_JetLowestMjj1_Id",&multilepton_JetLowestMjj1_Id,&b_multilepton_JetLowestMjj1_Id);
  tInput->SetBranchAddress("multilepton_JetLowestMjj1_P4",&multilepton_JetLowestMjj1_P4_ptr,&b_multilepton_JetLowestMjj1_P4);
  tInput->SetBranchAddress("multilepton_JetLowestMjj2_Id",&multilepton_JetLowestMjj2_Id,&b_multilepton_JetLowestMjj2_Id);
  tInput->SetBranchAddress("multilepton_JetLowestMjj2_P4",&multilepton_JetLowestMjj2_P4_ptr,&b_multilepton_JetLowestMjj2_P4);
  tInput->SetBranchAddress("multilepton_mET",&multilepton_mET_ptr,&b_multilepton_mET);
  tInput->SetBranchAddress("multilepton_Ptot",&multilepton_Ptot_ptr,&b_multilepton_Ptot);

  cout << "Creating output tree"<<endl;

  fOutput = new TFile("output.root", "RECREATE");
  tOutput = new TTree("Tree", "Tree");

  tOutput->Branch("mc_event",&mc_event,"mc_event/I");
  tOutput->Branch("mc_weight",&mc_weight,"mc_weight/F");
  tOutput->Branch("mc_3l_category",&mc_3l_category,"mc_3l_category/I");
  tOutput->Branch("mc_ttbar_decay",&mc_ttbar_decay,"mc_ttbar_decay/I");
  tOutput->Branch("mc_boson_decay",&mc_boson_decay,"mc_boson_decay/I");
  tOutput->Branch("mc_ttZhypAllowed",&mc_ttZhypAllowed,"mc_ttZhypAllowed/I");

  tOutput->Branch("multilepton_Bjet1_Id",&multilepton_Bjet1_Id,"multilepton_Bjet1_Id/I");
  tOutput->Branch("multilepton_Bjet1_P4","TLorentzVector",&multilepton_Bjet1_P4);
  tOutput->Branch("multilepton_Bjet2_Id",&multilepton_Bjet2_Id,"multilepton_Bjet2_Id/I");
  tOutput->Branch("multilepton_Bjet2_P4","TLorentzVector",&multilepton_Bjet2_P4);
  tOutput->Branch("multilepton_Lepton1_Id",&multilepton_Lepton1_Id,"multilepton_Lepton1_Id/I");
  tOutput->Branch("multilepton_Lepton1_P4","TLorentzVector",&multilepton_Lepton1_P4);
  tOutput->Branch("multilepton_Lepton2_Id",&multilepton_Lepton2_Id,"multilepton_Lepton2_Id/I");
  tOutput->Branch("multilepton_Lepton2_P4","TLorentzVector",&multilepton_Lepton2_P4);
  tOutput->Branch("multilepton_Lepton3_Id",&multilepton_Lepton3_Id,"multilepton_Lepton3_Id/I");
  tOutput->Branch("multilepton_Lepton3_P4","TLorentzVector",&multilepton_Lepton3_P4);
  tOutput->Branch("multilepton_JetHighestPt1_Id",&multilepton_JetHighestPt1_Id,"multilepton_JetHighestPt1_Id/I");
  tOutput->Branch("multilepton_JetHighestPt1_P4","TLorentzVector",&multilepton_JetHighestPt1_P4);
  tOutput->Branch("multilepton_JetHighestPt2_Id",&multilepton_JetHighestPt2_Id,"multilepton_JetHighestPt2_Id/I");
  tOutput->Branch("multilepton_JetHighestPt2_P4","TLorentzVector",&multilepton_JetHighestPt2_P4);
  tOutput->Branch("multilepton_JetClosestMw1_Id",&multilepton_JetClosestMw1_Id,"multilepton_JetClosestMw1_Id/I");
  tOutput->Branch("multilepton_JetClosestMw1_P4","TLorentzVector",&multilepton_JetClosestMw1_P4);
  tOutput->Branch("multilepton_JetClosestMw2_Id",&multilepton_JetClosestMw2_Id,"multilepton_JetClosestMw2_Id/I");
  tOutput->Branch("multilepton_JetClosestMw2_P4","TLorentzVector",&multilepton_JetClosestMw2_P4);
  tOutput->Branch("multilepton_JetLowestMjj1_Id",&multilepton_JetLowestMjj1_Id,"multilepton_JetLowestMjj1_Id/I");
  tOutput->Branch("multilepton_JetLowestMjj1_P4","TLorentzVector",&multilepton_JetLowestMjj1_P4);
  tOutput->Branch("multilepton_JetLowestMjj2_Id",&multilepton_JetLowestMjj2_Id,"multilepton_JetLowestMjj2_Id/I");
  tOutput->Branch("multilepton_JetLowestMjj2_P4","TLorentzVector",&multilepton_JetLowestMjj2_P4);
  tOutput->Branch("multilepton_mET","TLorentzVector",&multilepton_mET);
  tOutput->Branch("multilepton_Ptot","TLorentzVector",&multilepton_Ptot);

  tOutput->Branch("mc_mem_tthfl_weight",&mc_mem_tthfl_weight,"mc_mem_tthfl_weight/D");
  tOutput->Branch("mc_mem_tthfl_weight_log",&mc_mem_tthfl_weight_log,"mc_mem_tthfl_weight_log/D");
  tOutput->Branch("mc_mem_tthfl_weight_err",&mc_mem_tthfl_weight_err,"mc_mem_tthfl_weight_err/D");
  tOutput->Branch("mc_mem_tthfl_weight_chi2",&mc_mem_tthfl_weight_chi2,"mc_mem_tthfl_weight_chi2/F");
  tOutput->Branch("mc_mem_tthfl_weight_time",&mc_mem_tthfl_weight_time,"mc_mem_tthfl_weight_time/F");
  tOutput->Branch("mc_mem_tthfl_weight_max",&mc_mem_tthfl_weight_max,"mc_mem_tthfl_weight_max/D");
  tOutput->Branch("mc_mem_tthfl_weight_avg",&mc_mem_tthfl_weight_avg,"mc_mem_tthfl_weight_avg/D");

  tOutput->Branch("mc_mem_tthsl_weight",&mc_mem_tthsl_weight,"mc_mem_tthsl_weight/D");
  tOutput->Branch("mc_mem_tthsl_weight_log",&mc_mem_tthsl_weight_log,"mc_mem_tthsl_weight_log/D");
  tOutput->Branch("mc_mem_tthsl_weight_err",&mc_mem_tthsl_weight_err,"mc_mem_tthsl_weight_err/D");
  tOutput->Branch("mc_mem_tthsl_weight_chi2",&mc_mem_tthsl_weight_chi2,"mc_mem_tthsl_weight_chi2/F");
  tOutput->Branch("mc_mem_tthsl_weight_time",&mc_mem_tthsl_weight_time,"mc_mem_tthsl_weight_time/F");
  tOutput->Branch("mc_mem_tthsl_weight_max",&mc_mem_tthsl_weight_max,"mc_mem_tthsl_weight_max/D");
  tOutput->Branch("mc_mem_tthsl_weight_avg",&mc_mem_tthsl_weight_avg,"mc_mem_tthsl_weight_avg/D");

  tOutput->Branch("mc_mem_tth_weight",&mc_mem_tth_weight,"mc_mem_tth_weight/D");
  tOutput->Branch("mc_mem_tth_weight_log",&mc_mem_tth_weight_log,"mc_mem_tth_weight_log/D");
  tOutput->Branch("mc_mem_tth_weight_err",&mc_mem_tth_weight_err,"mc_mem_tth_weight_err/D");
  tOutput->Branch("mc_mem_tth_weight_chi2",&mc_mem_tth_weight_chi2,"mc_mem_tth_weight_chi2/F");
  tOutput->Branch("mc_mem_tth_weight_time",&mc_mem_tth_weight_time,"mc_mem_tth_weight_time/F");
  tOutput->Branch("mc_mem_tth_weight_max",&mc_mem_tth_weight_max,"mc_mem_tth_weight_max/D");
  tOutput->Branch("mc_mem_tth_weight_avg",&mc_mem_tth_weight_avg,"mc_mem_tth_weight_avg/D");

  tOutput->Branch("mc_mem_ttz_weight",&mc_mem_ttz_weight,"mc_mem_ttz_weight/D");
  tOutput->Branch("mc_mem_ttz_weight_log",&mc_mem_ttz_weight_log,"mc_mem_ttz_weight_log/D");
  tOutput->Branch("mc_mem_ttz_weight_err",&mc_mem_ttz_weight_err,"mc_mem_ttz_weight_err/D");
  tOutput->Branch("mc_mem_ttz_weight_chi2",&mc_mem_ttz_weight_chi2,"mc_mem_ttz_weight_chi2/F");
  tOutput->Branch("mc_mem_ttz_weight_time",&mc_mem_ttz_weight_time,"mc_mem_ttz_weight_time/F");
  tOutput->Branch("mc_mem_ttz_weight_max",&mc_mem_ttz_weight_max,"mc_mem_ttz_weight_max/D");
  tOutput->Branch("mc_mem_ttz_weight_avg",&mc_mem_ttz_weight_avg,"mc_mem_ttz_weight_avg/D");

  tOutput->Branch("mc_mem_ttw_weight",&mc_mem_ttw_weight,"mc_mem_ttw_weight/D");
  tOutput->Branch("mc_mem_ttw_weight_log",&mc_mem_ttw_weight_log,"mc_mem_ttw_weight_log/D");
  tOutput->Branch("mc_mem_ttw_weight_err",&mc_mem_ttw_weight_err,"mc_mem_ttw_weight_err/D");
  tOutput->Branch("mc_mem_ttw_weight_chi2",&mc_mem_ttw_weight_chi2,"mc_mem_ttw_weight_chi2/F");
  tOutput->Branch("mc_mem_ttw_weight_time",&mc_mem_ttw_weight_time,"mc_mem_ttw_weight_time/F");
  tOutput->Branch("mc_mem_ttw_weight_max",&mc_mem_ttw_weight_max,"mc_mem_ttw_weight_max/D");
  tOutput->Branch("mc_mem_ttw_weight_avg",&mc_mem_ttw_weight_avg,"mc_mem_ttw_weight_avg/D");

  tOutput->Branch("mc_mem_ttwjj_weight",&mc_mem_ttwjj_weight,"mc_mem_ttwjj_weight/D");
  tOutput->Branch("mc_mem_ttwjj_weight_log",&mc_mem_ttwjj_weight_log,"mc_mem_ttwjj_weight_log/D");
  tOutput->Branch("mc_mem_ttwjj_weight_err",&mc_mem_ttwjj_weight_err,"mc_mem_ttwjj_weight_err/D");
  tOutput->Branch("mc_mem_ttwjj_weight_chi2",&mc_mem_ttwjj_weight_chi2,"mc_mem_ttwjj_weight_chi2/F");
  tOutput->Branch("mc_mem_ttwjj_weight_time",&mc_mem_ttwjj_weight_time,"mc_mem_ttwjj_weight_time/F");
  tOutput->Branch("mc_mem_ttwjj_weight_max",&mc_mem_ttwjj_weight_max,"mc_mem_ttwjj_weight_max/D");
  tOutput->Branch("mc_mem_ttwjj_weight_avg",&mc_mem_ttwjj_weight_avg,"mc_mem_ttwjj_weight_avg/D");

  tOutput->Branch("mc_mem_ttz_tth_likelihood",&mc_mem_ttz_tth_likelihood,"mc_mem_ttz_tth_likelihood/D");
  tOutput->Branch("mc_mem_ttz_tth_likelihood_nlog",&mc_mem_ttz_tth_likelihood_nlog,"mc_mem_ttz_tth_likelihood_nlog/D");
  tOutput->Branch("mc_mem_ttz_tth_likelihood_max",&mc_mem_ttz_tth_likelihood_max,"mc_mem_ttz_tth_likelihood_max/D");
  tOutput->Branch("mc_mem_ttz_tth_likelihood_avg",&mc_mem_ttz_tth_likelihood_avg,"mc_mem_ttz_tth_likelihood_avg/D");

  tOutput->Branch("mc_mem_ttw_tth_likelihood",&mc_mem_ttw_tth_likelihood,"mc_mem_ttw_tth_likelihood/D");
  tOutput->Branch("mc_mem_ttw_tth_likelihood_nlog",&mc_mem_ttw_tth_likelihood_nlog,"mc_mem_ttw_tth_likelihood_nlog/D");
  tOutput->Branch("mc_mem_ttw_tth_likelihood_max",&mc_mem_ttw_tth_likelihood_max,"mc_mem_ttw_tth_likelihood_max/D");
  tOutput->Branch("mc_mem_ttw_tth_likelihood_avg",&mc_mem_ttw_tth_likelihood_avg,"mc_mem_ttw_tth_likelihood_avg/D");

  tOutput->Branch("mc_mem_ttwjj_tth_likelihood",&mc_mem_ttwjj_tth_likelihood,"mc_mem_ttwjj_tth_likelihood/D");
  tOutput->Branch("mc_mem_ttwjj_tth_likelihood_nlog",&mc_mem_ttwjj_tth_likelihood_nlog,"mc_mem_ttwjj_tth_likelihood_nlog/D");
  tOutput->Branch("mc_mem_ttwjj_tth_likelihood_max",&mc_mem_ttwjj_tth_likelihood_max,"mc_mem_ttwjj_tth_likelihood_max/D");
  tOutput->Branch("mc_mem_ttwjj_tth_likelihood_avg",&mc_mem_ttwjj_tth_likelihood_avg,"mc_mem_ttwjj_tth_likelihood_avg/D");

  tOutput->Branch("mc_mem_ttv_tth_likelihood",&mc_mem_ttv_tth_likelihood,"mc_mem_ttv_tth_likelihood/D");
  tOutput->Branch("mc_mem_ttv_tth_likelihood_nlog",&mc_mem_ttv_tth_likelihood_nlog,"mc_mem_ttv_tth_likelihood_nlog/D");
  tOutput->Branch("mc_mem_ttv_tth_likelihood_max",&mc_mem_ttv_tth_likelihood_max,"mc_mem_ttv_tth_likelihood_max/D");
  tOutput->Branch("mc_mem_ttv_tth_likelihood_avg",&mc_mem_ttv_tth_likelihood_avg,"mc_mem_ttv_tth_likelihood_avg/D");

  tOutput->Branch("mc_mem_ttvjj_tth_likelihood",&mc_mem_ttvjj_tth_likelihood,"mc_mem_ttvjj_tth_likelihood/D");
  tOutput->Branch("mc_mem_ttvjj_tth_likelihood_nlog",&mc_mem_ttvjj_tth_likelihood_nlog,"mc_mem_ttvjj_tth_likelihood_nlog/D");
  tOutput->Branch("mc_mem_ttvjj_tth_likelihood_max",&mc_mem_ttvjj_tth_likelihood_max,"mc_mem_ttvjj_tth_likelihood_max/D");
  tOutput->Branch("mc_mem_ttvjj_tth_likelihood_avg",&mc_mem_ttvjj_tth_likelihood_avg,"mc_mem_ttvjj_tth_likelihood_avg/D");
}


void ReadGenFlatTree::FillGenMultilepton(Long64_t iEvent, MultiLepton* multiLepton){

  cout << "FillGenMultilepton"<<endl;

  TLorentzVector Ptot, Pll, Pl1, PtotNeut(0,0,0,0), Pjj;

  tInput->LoadTree(iEvent);
  tInput->GetEntry(iEvent);

  mc_ttbar_decay = -1;
  mc_boson_decay = -1;

  if (mc_truth_t1_id==-666 || mc_truth_t2_id==-666) return;

  (*multiLepton).Leptons.clear();
  (*multiLepton).Jets.clear();
  (*multiLepton).Bjets.clear();
  (*multiLepton).AllJets.clear();
  (*multiLepton).JetsHighestPt.clear();
  (*multiLepton).JetsClosestMw.clear();
  (*multiLepton).JetsLowestMjj.clear();

    if (mc_truth_h0_id>-666 && mc_truth_h0Wl1_id>-666 && mc_truth_h0Wl2_id>-666){
      (*multiLepton).FillParticle("lepton", mc_truth_h0Wl1_id, *mc_truth_h0Wl1_p4);
      (*multiLepton).FillParticle("lepton", mc_truth_h0Wl2_id, *mc_truth_h0Wl2_p4);
      Ptot = (*mc_truth_h0_p4) + (*mc_truth_t1_p4) + (*mc_truth_t2_p4);
      Pll = (*mc_truth_h0Wl1_p4) + (*mc_truth_h0Wl2_p4);
      mc_boson_decay = 0;
      mc_boson_pt = mc_truth_h0_p4->Pt();
      mc_boson_l1_pt = mc_truth_h0Wl1_p4->Pt();
      mc_boson_l1_eta = mc_truth_h0Wl1_p4->Eta();
      mc_boson_l2_pt = mc_truth_h0Wl2_p4->Pt();
      mc_boson_l2_eta = mc_truth_h0Wl2_p4->Eta();
      mc_boson_ll_mass = Pll.M();
      mc_boson_ll_pt = Pll.Pt();
      mc_boson_ll_dphi = TMath::Abs(mc_truth_h0Wl1_p4->DeltaPhi(*mc_truth_h0Wl2_p4));
      PtotNeut += (*mc_truth_h0Wnu1_p4);
      PtotNeut += (*mc_truth_h0Wnu2_p4);
    }
    else if (mc_truth_h0_id>-666 && mc_truth_h0Wq11_id>-666 && mc_truth_h0Wl2_id>-666){
      (*multiLepton).FillParticle("lepton", mc_truth_h0Wl2_id, *mc_truth_h0Wl2_p4);
      mc_boson_decay = 1;
      Ptot = (*mc_truth_h0_p4) + (*mc_truth_t1_p4) + (*mc_truth_t2_p4);
      Pjj = (*mc_truth_h0Wq11_p4) + (*mc_truth_h0Wq12_p4);
      mc_boson_pt = mc_truth_h0_p4->Pt();
      mc_boson_l1_pt = mc_truth_h0Wl1_p4->Pt();
      mc_boson_l1_eta = mc_truth_h0Wl1_p4->Eta();
      mc_boson_l2_pt = -666;
      mc_boson_l2_eta = -666;
      mc_boson_ll_mass = -666;
      mc_boson_ll_pt = -666;
      mc_boson_j1_pt = mc_truth_h0Wq11_p4->Pt();
      mc_boson_j1_eta = mc_truth_h0Wq11_p4->Eta();
      mc_boson_j2_pt = mc_truth_h0Wq12_p4->Pt();
      mc_boson_j2_eta = mc_truth_h0Wq12_p4->Eta();
      mc_boson_jj_mass = Pjj.M();
      mc_boson_jj_pt = Pjj.Pt();
      mc_boson_jj_dphi = TMath::Abs(mc_truth_h0Wq11_p4->DeltaPhi(*mc_truth_h0Wq12_p4));
      PtotNeut += (*mc_truth_h0Wnu1_p4);
    }
    else if (mc_truth_h0_id>-666 && mc_truth_h0Wq21_id>-666 && mc_truth_h0Wl1_id>-666){
      (*multiLepton).FillParticle("lepton", mc_truth_h0Wl1_id, *mc_truth_h0Wl1_p4);
      mc_boson_decay = 1;
      Ptot = (*mc_truth_h0_p4) + (*mc_truth_t1_p4) + (*mc_truth_t2_p4);
      Pjj = (*mc_truth_h0Wq21_p4) + (*mc_truth_h0Wq22_p4);
      mc_boson_pt = mc_truth_h0_p4->Pt();
      mc_boson_l1_pt = mc_truth_h0Wl1_p4->Pt();
      mc_boson_l1_eta = mc_truth_h0Wl1_p4->Eta();
      mc_boson_l2_pt = -666;
      mc_boson_l2_eta = -666;
      mc_boson_ll_mass = -666;
      mc_boson_ll_pt = -666;
      mc_boson_j1_pt = mc_truth_h0Wq21_p4->Pt();
      mc_boson_j1_eta = mc_truth_h0Wq21_p4->Eta();
      mc_boson_j2_pt = mc_truth_h0Wq22_p4->Pt();
      mc_boson_j2_eta = mc_truth_h0Wq22_p4->Eta();
      mc_boson_jj_mass = Pjj.M();
      mc_boson_jj_pt = Pjj.Pt();
      mc_boson_jj_dphi = TMath::Abs(mc_truth_h0Wq21_p4->DeltaPhi(*mc_truth_h0Wq22_p4));
      PtotNeut += (*mc_truth_h0Wnu1_p4);
    }
    else if (mc_truth_W_id>-666 && mc_truth_Wl_id>-666){
      //cout << "TTW"<<endl;
      (*multiLepton).FillParticle("lepton", mc_truth_Wl_id, *mc_truth_Wl_p4);
      mc_boson_decay = 3;
      Ptot = (*mc_truth_W_p4) + (*mc_truth_t1_p4) + (*mc_truth_t2_p4);
      mc_boson_pt = mc_truth_W_p4->Pt();
      mc_boson_l1_pt = mc_truth_Wl_p4->Pt();
      mc_boson_l1_eta = mc_truth_Wl_p4->Eta();
    }
    else if (mc_truth_Z_id>-666 && mc_truth_Zl1_id>-666 && mc_truth_Zl2_id>-666){
      (*multiLepton).FillParticle("lepton", mc_truth_Zl1_id, *mc_truth_Zl1_p4);
      (*multiLepton).FillParticle("lepton", mc_truth_Zl2_id, *mc_truth_Zl2_p4);
      mc_boson_decay = 2;
      Ptot = (*mc_truth_Zl1_p4) + (*mc_truth_Zl2_p4) + (*mc_truth_t1_p4) + (*mc_truth_t2_p4);
      Pll = (*mc_truth_Zl1_p4) + (*mc_truth_Zl2_p4);
      mc_boson_pt = Pll.Pt();
      mc_boson_l1_pt = mc_truth_Zl1_p4->Pt();
      mc_boson_l1_eta = mc_truth_Zl1_p4->Eta();
      mc_boson_l2_pt = mc_truth_Zl2_p4->Pt();
      mc_boson_l2_eta = mc_truth_Zl2_p4->Eta();
      mc_boson_ll_mass = Pll.M();
      mc_boson_ll_pt = Pll.Pt();
      mc_boson_ll_dphi = TMath::Abs(mc_truth_Zl1_p4->DeltaPhi(*mc_truth_Zl2_p4));
    }
    else if (mc_truth_gammal1_id>-666 && mc_truth_gammal2_id>-666){
      (*multiLepton).FillParticle("lepton", mc_truth_gammal1_id, *mc_truth_gammal1_p4);
      (*multiLepton).FillParticle("lepton", mc_truth_gammal2_id, *mc_truth_gammal2_p4);
      mc_boson_decay = 2;
      Ptot = (*mc_truth_gammal1_p4) + (*mc_truth_gammal2_p4) + (*mc_truth_t1_p4) + (*mc_truth_t2_p4);
      Pll = (*mc_truth_gammal1_p4) + (*mc_truth_gammal2_p4);
      mc_boson_pt = Pll.Pt();
      mc_boson_l1_pt = mc_truth_gammal1_p4->Pt();
      mc_boson_l1_eta = mc_truth_gammal1_p4->Eta();
      mc_boson_l2_pt = mc_truth_gammal2_p4->Pt();
      mc_boson_l2_eta = mc_truth_gammal2_p4->Eta();
      mc_boson_ll_mass = Pll.M();
      mc_boson_ll_pt = Pll.Pt();
      mc_boson_ll_dphi = TMath::Abs(mc_truth_gammal1_p4->DeltaPhi(*mc_truth_gammal2_p4));
    }
    else return;

    if (mc_truth_tWq11_id>-666 && mc_truth_tWl2_id>-666){
      //cout << "Semileptonic ttbar"<<endl;
      (*multiLepton).FillParticle("lepton", mc_truth_tWl2_id, *mc_truth_tWl2_p4);
      (*multiLepton).FillParticle("bjet", mc_truth_tb1_id, *mc_truth_tb1_p4);
      (*multiLepton).FillParticle("bjet", mc_truth_tb2_id, *mc_truth_tb2_p4);
      mc_ttbar_decay = 1;
      mc_thad_pt = mc_truth_t1_p4->Pt();
      mc_thad_b_pt = mc_truth_tb1_p4->Pt();
      mc_thad_b_eta = mc_truth_tb1_p4->Eta();
      mc_thad_j1_pt = mc_truth_tWq11_p4->Pt();
      mc_thad_j1_eta = mc_truth_tWq11_p4->Eta();
      mc_thad_j2_pt = mc_truth_tWq12_p4->Pt();
      mc_thad_j2_eta = mc_truth_tWq12_p4->Eta();
      mc_tlep_pt = mc_truth_t2_p4->Pt();
      mc_tlep_b_pt = mc_truth_tb2_p4->Pt();
      mc_tlep_b_eta = mc_truth_tb2_p4->Eta();
      mc_tlep_l_pt = mc_truth_tWl2_p4->Pt();
      mc_tlep_l_eta = mc_truth_tWl2_p4->Eta();
      PtotNeut += (*mc_truth_tWnu2_p4);
   }
   else if (mc_truth_tWq21_id>-666 && mc_truth_tWl1_id>-666){
      //cout << "Semileptonic ttbar"<<endl;
      (*multiLepton).FillParticle("lepton", mc_truth_tWl1_id, *mc_truth_tWl1_p4);
      (*multiLepton).FillParticle("bjet", mc_truth_tb1_id, *mc_truth_tb1_p4);
      (*multiLepton).FillParticle("bjet", mc_truth_tb2_id, *mc_truth_tb2_p4);
      mc_ttbar_decay = 1;
      mc_thad_pt = mc_truth_t2_p4->Pt();
      mc_thad_b_pt = mc_truth_tb2_p4->Pt();
      mc_thad_b_eta = mc_truth_tb2_p4->Eta();
      mc_thad_j1_pt = mc_truth_tWq21_p4->Pt();
      mc_thad_j1_eta = mc_truth_tWq21_p4->Eta();
      mc_thad_j2_pt = mc_truth_tWq22_p4->Pt();
      mc_thad_j2_eta = mc_truth_tWq22_p4->Eta();
      mc_tlep_pt = mc_truth_t1_p4->Pt();
      mc_tlep_b_pt = mc_truth_tb1_p4->Pt();
      mc_tlep_b_eta = mc_truth_tb1_p4->Eta();
      mc_tlep_l_pt = mc_truth_tWl1_p4->Pt();
      mc_tlep_l_eta = mc_truth_tWl1_p4->Eta();
      PtotNeut += (*mc_truth_tWnu1_p4);
   }
   else if (mc_truth_tWl1_id>-666 && mc_truth_tWl2_id>-666 && (mc_boson_decay==1 || mc_boson_decay==3)){
      //cout << "Fully leptonic ttbar"<<endl;
      (*multiLepton).FillParticle("bjet", mc_truth_tb1_id, *mc_truth_tb1_p4);
      (*multiLepton).FillParticle("bjet", mc_truth_tb2_id, *mc_truth_tb2_p4);
      (*multiLepton).FillParticle("lepton", mc_truth_tWl1_id, *mc_truth_tWl1_p4);
      (*multiLepton).FillParticle("lepton", mc_truth_tWl2_id, *mc_truth_tWl2_p4);
      mc_ttbar_decay = 2;
      mc_thad_pt = -666;
      mc_thad_b_pt = -666;
      mc_thad_b_eta = -666;
      mc_thad_j1_pt = -666;
      mc_thad_j1_eta = -666;
      mc_thad_j2_pt = -666;
      mc_thad_j2_eta = -666;
      mc_tlep_pt = mc_truth_t1_p4->Pt();
      mc_tlep_b_pt = mc_truth_tb1_p4->Pt();
      mc_tlep_b_eta = mc_truth_tb1_p4->Eta();
      mc_tlep_l_pt = mc_truth_tWl1_p4->Pt();
      mc_tlep_l_eta = mc_truth_tWl1_p4->Eta();
      mc_tlep2_pt = mc_truth_t2_p4->Pt();
      mc_tlep2_b_eta = mc_truth_tb2_p4->Eta();
      mc_tlep2_l_pt = mc_truth_tWl2_p4->Pt();
      mc_tlep2_l_eta = mc_truth_tWl2_p4->Eta();
      PtotNeut =PtotNeut + (*mc_truth_tWnu1_p4)+ (*mc_truth_tWnu2_p4);
   }
   else return;


   TLorentzVector Pjet;
   for (int i=0; i<genJet_n; i++){
     //cout << "genJet "<<i<<" pt="<<genJet_pt->at(i)<<endl;
     if (TMath::Abs(genJet_eta->at(i))>2.5) continue;
     if (genJet_pt->at(i)<25) continue;
     Pjet.SetPtEtaPhiE(genJet_pt->at(i), genJet_eta->at(i), genJet_phi->at(i), genJet_E->at(i));
     if ((*multiLepton).Bjets[0].P4.DeltaR(Pjet)<0.4 || (*multiLepton).Bjets[1].P4.DeltaR(Pjet)<0.4 || (*multiLepton).Leptons[0].P4.DeltaR(Pjet)<0.4 || (*multiLepton).Leptons[1].P4.DeltaR(Pjet)<0.4 || (*multiLepton).Leptons[2].P4.DeltaR(Pjet)<0.4) continue;
     (*multiLepton).FillParticle("alljet", 0, Pjet);
   }

   if ((*multiLepton).AllJets.size()==0) mc_3l_category = 0;
   else if ((*multiLepton).AllJets.size()==1) {
     mc_3l_category = 1;
     //(*multiLepton).FillParticle("jet", 0, Pjet);
   }
   else if ((*multiLepton).AllJets.size()>=2) {
      mc_3l_category = 2;

      TLorentzVector Pjet2;
      float pt_max=0, pt_max2=0; int ij1=-1, ij2=-1;
      float diffmass_min = 10000, mass_min = 10000; int ik1=-1, ik2=-1, il1=-1, il2=-1;
      for (unsigned int ij=0; ij<(*multiLepton).AllJets.size(); ij++){
        if ((*multiLepton).AllJets[ij].P4.Pt() > pt_max ) {
           pt_max2 = pt_max;
           ij2 = ij1;
           pt_max = (*multiLepton).AllJets[ij].P4.Pt();
           ij1 = ij;
         }
         if ((*multiLepton).AllJets[ij].P4.Pt() < pt_max && (*multiLepton).AllJets[ij].P4.Pt() > pt_max2){
           pt_max2 = (*multiLepton).AllJets[ij].P4.Pt();
           ij2 = ij;
         }
         for (unsigned int ik=0; ik<(*multiLepton).AllJets.size(); ik++){
           if (ik==ij) continue;
           if (TMath::Abs(((*multiLepton).AllJets[ij].P4+(*multiLepton).AllJets[ik].P4).M()-80.419)<diffmass_min){
             ik1=ij;
             ik2=ik;
             diffmass_min = TMath::Abs(((*multiLepton).AllJets[ij].P4+(*multiLepton).AllJets[ik].P4).M()-80.419);
           }
	   if (((*multiLepton).AllJets[ij].P4+(*multiLepton).AllJets[ik].P4).M()<mass_min){
	     il1=ij;
             il2=ik;
	     mass_min = ((*multiLepton).AllJets[ij].P4+(*multiLepton).AllJets[ik].P4).M();
	   }
         }
      }
      if (ij1!=-1 && ij2!=-1) {
        (*multiLepton).FillParticle("jetHighestPt", 1, (*multiLepton).AllJets[ij1].P4);
        (*multiLepton).FillParticle("jetHighestPt", 1, (*multiLepton).AllJets[ij2].P4);
      }
      if (ik1!=-1 && ik2!=-1){
        (*multiLepton).FillParticle("jetClosestMw", 2, (*multiLepton).AllJets[ik1].P4);
        (*multiLepton).FillParticle("jetClosestMw", 2, (*multiLepton).AllJets[ik2].P4);
      }
      if (il1!=-1 && il2!=-1){
        (*multiLepton).FillParticle("jetLowestMjj", 3, (*multiLepton).AllJets[il1].P4);
        (*multiLepton).FillParticle("jetLowestMjj", 3, (*multiLepton).AllJets[il2].P4);
      }
   }

   mc_totp4_px = Ptot.Px();
   mc_totp4_py = Ptot.Py();
   mc_totp4_pt = Ptot.Pt();
   mc_met = PtotNeut.Pt();

   mc_njets25 = (*multiLepton).AllJets.size();

   (*multiLepton).Ptot = Ptot;
   (*multiLepton).mET = PtotNeut;

}

int ReadGenFlatTree::ApplyGenSelection(Long64_t iEvent, MultiLepton* multiLepton){

  cout << "iEvent="<<iEvent<<" BosonDecay="<<mc_boson_decay<< " TTBarDecay="<< mc_ttbar_decay<<" mc_weigth="<<mc_weight<<" PtTot="<< (*multiLepton).Ptot.Pt()<<endl;
  //cout << "MultiLepton nLepton="<<(*multiLepton).Leptons.size()<<" nJets="<<(*multiLepton).Jets.size()<<" nBjets="<<(*multiLepton).Bjets.size()<<endl;


  if (!(((mc_boson_decay==0 || mc_boson_decay==2) && mc_ttbar_decay==1) //ttHfl, ttZ
      || (mc_boson_decay==1 && mc_ttbar_decay==2) //ttHsl
      || (mc_boson_decay==3 && mc_ttbar_decay==2) //ttW
      )) return 0;
 
  if (!((*multiLepton).Leptons.size()==3 && (*multiLepton).Bjets.size()==2 )) return 0;

  if (mc_3l_category!=2) return 0;

   mc_ttZhypAllowed = 0;
   mc_hasLLcombZpeak = -1;
   mc_passMllGt12 = -1;

   mc_passMllGt12 = 1;
   if (((*multiLepton).Leptons.at(0).P4+(*multiLepton).Leptons.at(1).P4).M()<12) mc_passMllGt12 = 0;
   if (((*multiLepton).Leptons.at(0).P4+(*multiLepton).Leptons.at(2).P4).M()<12) mc_passMllGt12 = 0; 
   if (((*multiLepton).Leptons.at(1).P4+(*multiLepton).Leptons.at(2).P4).M()<12) mc_passMllGt12 = 0; 

   if ((*multiLepton).Leptons.at(0).Id==-(*multiLepton).Leptons.at(1).Id || (*multiLepton).Leptons.at(0).Id==-(*multiLepton).Leptons.at(2).Id || (*multiLepton).Leptons.at(1).Id==-(*multiLepton).Leptons.at(2).Id){

     mc_ttZhypAllowed = 1;
     mc_hasLLcombZpeak = 0;
     mc_passMllGt12 = 1;
     if ((*multiLepton).Leptons.at(0).Id==-(*multiLepton).Leptons.at(1).Id){
        if (TMath::Abs(((*multiLepton).Leptons.at(0).P4+(*multiLepton).Leptons.at(1).P4).M()-91.188)<10) mc_hasLLcombZpeak = 1;
        //if (((*multiLepton).Leptons.at(0).P4+(*multiLepton).Leptons.at(1).P4).M()<12) mc_passMllGt12 = 0;
     }
     if ((*multiLepton).Leptons.at(0).Id==-(*multiLepton).Leptons.at(2).Id){
        if (TMath::Abs(((*multiLepton).Leptons.at(0).P4+(*multiLepton).Leptons.at(2).P4).M()-91.188)<10) mc_hasLLcombZpeak = 1;
        //if (((*multiLepton).Leptons.at(0).P4+(*multiLepton).Leptons.at(2).P4).M()<12) mc_passMllGt12 = 0;
     }
     if ((*multiLepton).Leptons.at(1).Id==-(*multiLepton).Leptons.at(2).Id){
        if (TMath::Abs(((*multiLepton).Leptons.at(1).P4+(*multiLepton).Leptons.at(2).P4).M()-91.188)<10) mc_hasLLcombZpeak = 1;
        //if (((*multiLepton).Leptons.at(1).P4+(*multiLepton).Leptons.at(2).P4).M()<12) mc_passMllGt12 = 0;
     }
   }

   //cout << "Sort particles"<<endl;
   (*multiLepton).DoSort(&(*multiLepton).Leptons);
   (*multiLepton).DoSort(&(*multiLepton).Jets);
   (*multiLepton).DoSort(&(*multiLepton).Bjets);

   cout << "After sorting, Lepton0Pt="<<(*multiLepton).Leptons.at(0).P4.Pt()<<" Lepton1Pt="<<(*multiLepton).Leptons.at(1).P4.Pt() << " Lepton2Pt="<<(*multiLepton).Leptons.at(2).P4.Pt()<<endl;

   mc_passLepPresel = true;
   if (TMath::Abs((*multiLepton).Leptons.at(0).P4.Eta())>2.5 || TMath::Abs((*multiLepton).Leptons.at(1).P4.Eta())>2.5 || TMath::Abs((*multiLepton).Leptons.at(2).P4.Eta())>2.5) mc_passLepPresel=false;
   if ((*multiLepton).Leptons.at(0).P4.Pt()<10 || (*multiLepton).Leptons.at(1).P4.Pt()<10 || (*multiLepton).Leptons.at(2).P4.Pt()<10) mc_passLepPresel=false;
   if ((*multiLepton).Leptons.at(0).P4.Pt()<20) mc_passLepPresel=false;

   mc_passJetPresel25 = true;
   if (mc_3l_category!=2) mc_passJetPresel25 = false;

   mc_passBjetPresel25 = true;
   //if (TMath::Abs((*multiLepton).Jets.at(0).P4.Eta())>2.5 || TMath::Abs((*multiLepton).Jets.at(1).P4.Eta())>2.5) mc_passJetPresel25 = false;
   //if ((*multiLepton).Jets.at(0).P4.Pt()<25 || (*multiLepton).Jets.at(1).P4.Pt()<25) mc_passJetPresel25 = false;
   if (TMath::Abs((*multiLepton).Bjets.at(0).P4.Eta())>2.5 || TMath::Abs((*multiLepton).Bjets.at(1).P4.Eta())>2.5) mc_passBjetPresel25 = false;
   if ((*multiLepton).Bjets.at(0).P4.Pt()<25 || (*multiLepton).Bjets.at(1).P4.Pt()<25) mc_passBjetPresel25 = false;

   if (mc_passLepPresel && mc_passBjetPresel25 && mc_passMllGt12 && (mc_ttZhypAllowed==0 || (mc_ttZhypAllowed==1 && mc_hasLLcombZpeak!=1)))
     return 1;
   else return 0;

}

void ReadGenFlatTree::WriteMultilepton(MultiLepton* multiLepton){

  cout << "WriteMultilepton"<<endl;

  multilepton_Bjet1_Id = (*multiLepton).Bjets[0].Id;
  multilepton_Bjet1_P4 = (*multiLepton).Bjets[0].P4;
  multilepton_Bjet2_Id = (*multiLepton).Bjets[1].Id;
  multilepton_Bjet2_P4 = (*multiLepton).Bjets[1].P4;
  multilepton_Lepton1_Id = (*multiLepton).Leptons[0].Id;
  multilepton_Lepton1_P4 = (*multiLepton).Leptons[0].P4;
  multilepton_Lepton2_Id = (*multiLepton).Leptons[1].Id;
  multilepton_Lepton2_P4 = (*multiLepton).Leptons[1].P4;
  multilepton_Lepton3_Id = (*multiLepton).Leptons[2].Id;
  multilepton_Lepton3_P4 = (*multiLepton).Leptons[2].P4;
  multilepton_JetHighestPt1_Id = (*multiLepton).JetsHighestPt[0].Id;
  multilepton_JetHighestPt1_P4 = (*multiLepton).JetsHighestPt[0].P4;
  multilepton_JetHighestPt2_Id = (*multiLepton).JetsHighestPt[1].Id;
  multilepton_JetHighestPt2_P4 = (*multiLepton).JetsHighestPt[1].P4;
  multilepton_JetClosestMw1_Id = (*multiLepton).JetsClosestMw[0].Id;
  multilepton_JetClosestMw1_P4 = (*multiLepton).JetsClosestMw[0].P4;
  multilepton_JetClosestMw2_Id = (*multiLepton).JetsClosestMw[1].Id;
  multilepton_JetClosestMw2_P4 = (*multiLepton).JetsClosestMw[1].P4;
  multilepton_JetLowestMjj1_Id = (*multiLepton).JetsLowestMjj[0].Id;
  multilepton_JetLowestMjj1_P4 = (*multiLepton).JetsLowestMjj[0].P4;
  multilepton_JetLowestMjj2_Id = (*multiLepton).JetsLowestMjj[1].Id;
  multilepton_JetLowestMjj2_P4 = (*multiLepton).JetsLowestMjj[1].P4;
  multilepton_mET = (*multiLepton).mET;
  multilepton_Ptot = (*multiLepton).Ptot;

  multilepton_JetHighestPt_Mjj = (multilepton_JetHighestPt1_P4+multilepton_JetHighestPt2_P4).M();
  multilepton_JetClosestMw_Mjj = (multilepton_JetClosestMw1_P4+multilepton_JetClosestMw2_P4).M();
  multilepton_JetLowestMjj_Mjj = (multilepton_JetLowestMjj1_P4+multilepton_JetLowestMjj2_P4).M();


  return;
}

void ReadGenFlatTree::ReadMultilepton(Long64_t iEvent, MultiLepton* multiLepton){

  cout << "ReadMultilepton"<<endl;

  tInput->LoadTree(iEvent);
  tInput->GetEntry(iEvent);

  (*multiLepton).Leptons.clear();
  (*multiLepton).FillParticle("lepton", multilepton_Lepton1_Id, *multilepton_Lepton1_P4_ptr);
  (*multiLepton).FillParticle("lepton", multilepton_Lepton2_Id, *multilepton_Lepton2_P4_ptr);
  (*multiLepton).FillParticle("lepton", multilepton_Lepton3_Id, *multilepton_Lepton3_P4_ptr);

  (*multiLepton).JetsHighestPt.clear();
  (*multiLepton).FillParticle("jetHighestPt", multilepton_JetHighestPt1_Id, *multilepton_JetHighestPt1_P4_ptr);
  (*multiLepton).FillParticle("jetHighestPt", multilepton_JetHighestPt2_Id, *multilepton_JetHighestPt2_P4_ptr);

  (*multiLepton).JetsClosestMw.clear();
  (*multiLepton).FillParticle("jetClosestMw", multilepton_JetClosestMw1_Id, *multilepton_JetClosestMw1_P4_ptr);
  (*multiLepton).FillParticle("jetClosestMw", multilepton_JetClosestMw2_Id, *multilepton_JetClosestMw2_P4_ptr);

  (*multiLepton).JetsLowestMjj.clear();
  (*multiLepton).FillParticle("jetLowestMjj", multilepton_JetLowestMjj1_Id, *multilepton_JetLowestMjj1_P4_ptr);
  (*multiLepton).FillParticle("jetLowestMjj", multilepton_JetLowestMjj2_Id, *multilepton_JetLowestMjj2_P4_ptr);

  (*multiLepton).Bjets.clear();
  (*multiLepton).FillParticle("bjet", multilepton_Bjet1_Id, *multilepton_Bjet1_P4_ptr);
  (*multiLepton).FillParticle("bjet", multilepton_Bjet2_Id, *multilepton_Bjet2_P4_ptr);

  (*multiLepton).Ptot = *multilepton_Ptot_ptr;
  (*multiLepton).mET = *multilepton_mET_ptr;

  //cout << "Lepton0Pt="<<(*multiLepton).Leptons.at(0).P4.Pt()<<" Lepton1Pt="<<(*multiLepton).Leptons.at(1).P4.Pt() << " Lepton2Pt="<<(*multiLepton).Leptons.at(2).P4.Pt()<<endl;  
  //cout << "Bjet0Pt="<<(*multiLepton).Bjets.at(0).P4.Pt()<<" Bjet1Pt="<<(*multiLepton).Bjets.at(1).P4.Pt() << endl;
  //cout << "JetHighestPt0Pt="<<(*multiLepton).JetsHighestPt.at(0).P4.Pt() << " JetHighestPt1Pt="<<(*multiLepton).JetsHighestPt.at(1).P4.Pt() << endl;
  //cout << "JetClosestMw0Pt="<<(*multiLepton).JetsClosestMw.at(0).P4.Pt() << " JetClosestMw1Pt="<<(*multiLepton).JetsClosestMw.at(1).P4.Pt() << endl;
  //cout << "JetLowestMjj0Pt="<<(*multiLepton).JetsLowestMjj.at(0).P4.Pt() << " JetLowestMjj1Pt="<<(*multiLepton).JetsLowestMjj.at(1).P4.Pt() << endl;

  return;
}


#endif
