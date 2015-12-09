
#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

#include "Math/Integrator.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/IOptions.h"
#include "Math/IntegratorOptions.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/AdaptiveIntegratorMultiDim.h"

#include "MEPhaseSpace.h"
#include "HypIntegrator.h"
#include "MultiLepton.h"
 
#include <iostream>
#include <ctime>

using namespace std;

int main(int argc, char *argv[])
{

  cout << "Test MEM weight computation" <<endl;
  cout << "Will run on "<<argv[1]<<" from event "<<argv[2]<<" to "<<argv[3]<<endl;

  HypIntegrator hypIntegrator;
  hypIntegrator.InitializeIntegrator(13000, kMadgraph, kTFGaussian, 10000);

  double xsTTslHl = hypIntegrator.meIntegrator->xsTTH * hypIntegrator.meIntegrator->brTopHad * hypIntegrator.meIntegrator->brTopLep * hypIntegrator.meIntegrator->brHiggsFullLep;
  double xsTTLLsl = hypIntegrator.meIntegrator->xsTTLL * hypIntegrator.meIntegrator->brTopHad * hypIntegrator.meIntegrator->brTopLep;

  int nParam = 10;
  unsigned int intPoints = 1000000;

  hypIntegrator.meIntegrator->SetVerbosity(2);
  //hypIntegrator.meIntegrator->SetVerbosity(1);

  double mTop = 173;
  double mHiggs = 125;
  double mB = 4.7;

  double *xL = new double[nParam];
  double *xU = new double[nParam];

    xL[0] = mB; //TopHad, Bjet_E
    xL[1] = 0; //TopHad, Jet1_E
    xL[2] = 0; //TopLep, Bjet_E
    xL[3] = 0; //TopLep, Neut_Theta
    xL[4] = 0; //TopLep, Neut_Phi
    xL[5] = 0; //HiggsFullLep, Neut1_E
    xL[6] = 0; //HiggsFullLep, Neut1_Theta
    xL[7] = 0; //HiggsFullLep, Neut1_Phi
    xL[8] = 0; //HiggsFullLep, Neut2_Theta
    xL[9] = 0; //HiggsFullLep, Neut2_Phi

    xU[0] = 1000;//meIntegrator->comEnergy; //TopHad, Bjet_E
    xU[1] = 1000;//meIntegrator->comEnergy; //TopHad, Jet1_E
    xU[2] = 1000;//meIntegrator->comEnergy; //TopLep, Bjet_E
    xU[3] = TMath::Pi(); //TopLep, Neut_Theta
    xU[4] = 2*TMath::Pi(); //TopLep, Neut_Phi
    xU[5] = 1000;//meIntegrator->comEnergy; //HiggsFullLep, Neut1_E
    xU[6] = TMath::Pi(); //HiggsFullLep, Neut1_Theta
    xU[7] = 2*TMath::Pi(); //HiggsFullLep, Neut1_Phi
    xU[8] = TMath::Pi(); //HiggsFullLep, Neut2_Theta
    xU[9] = 2*TMath::Pi(); //HiggsFullLep, Neut2_Phi
 
  string InputFileName = string(argv[1]);
  TFile* f = TFile::Open(InputFileName.c_str());
  //TFile* f = new TFile("output_TTHaMCatNLO_Slimmed.root");
  //TFile* f = new TFile("../test/output_TTHaMCatNLO.root");
  //TFile* f = new TFile("../test/output_TTZMadgraph.root");
  TTree* t = (TTree*)f->Get("FlatTree/tree");

  Float_t mc_weight;
  Int_t mc_truth_h0_id;
  TLorentzVector* mc_truth_h0_p4 = 0;
  Int_t mc_truth_h0Wl1_id;
  TLorentzVector* mc_truth_h0Wl1_p4 = 0;
  Int_t mc_truth_h0Wl2_id;
  TLorentzVector* mc_truth_h0Wl2_p4 = 0;
  //Int_t mc_truth_h0Wnu1_id;
  TLorentzVector* mc_truth_h0Wnu1_p4 = 0;
  //Int_t mc_truth_h0Wnu2_id;
  TLorentzVector* mc_truth_h0Wnu2_p4 = 0;
  Int_t mc_truth_h0Wq11_id;
  TLorentzVector* mc_truth_h0Wq11_p4 = 0;
  //Int_t mc_truth_h0Wq12_id;
  TLorentzVector* mc_truth_h0Wq12_p4 = 0;
  Int_t mc_truth_h0Wq21_id;
  TLorentzVector* mc_truth_h0Wq21_p4 = 0;
  //Int_t mc_truth_h0Wq22_id;
  TLorentzVector* mc_truth_h0Wq22_p4 = 0;
  Int_t mc_truth_t1_id;
  TLorentzVector* mc_truth_t1_p4 = 0;
  Int_t mc_truth_t2_id;
  TLorentzVector* mc_truth_t2_p4 = 0;
  Int_t mc_truth_tb1_id;
  TLorentzVector* mc_truth_tb1_p4 = 0;
  Int_t mc_truth_tb2_id;
  TLorentzVector* mc_truth_tb2_p4 = 0;
  Int_t mc_truth_tWl1_id; 
  TLorentzVector* mc_truth_tWl1_p4 = 0;
  //Int_t mc_truth_tWnu1_id;
  TLorentzVector* mc_truth_tWnu1_p4 = 0;
  Int_t mc_truth_tWl2_id;
  TLorentzVector* mc_truth_tWl2_p4 = 0;
  //Int_t mc_truth_tWnu2_id;
  TLorentzVector* mc_truth_tWnu2_p4 = 0;
  Int_t mc_truth_tWq11_id;
  TLorentzVector* mc_truth_tWq11_p4 = 0;
  Int_t mc_truth_tWq21_id;
  TLorentzVector* mc_truth_tWq21_p4 = 0;
  Int_t mc_truth_tWq12_id;
  TLorentzVector* mc_truth_tWq12_p4 = 0;
  Int_t mc_truth_tWq22_id;
  TLorentzVector* mc_truth_tWq22_p4 = 0;
  Int_t mc_truth_Z_id;
  TLorentzVector* mc_truth_Z_p4 = 0;
  Int_t mc_truth_Zl1_id;
  TLorentzVector* mc_truth_Zl1_p4 = 0;
  Int_t mc_truth_Zl2_id;
  TLorentzVector* mc_truth_Zl2_p4 = 0;


  TBranch* b_mc_weight;
  TBranch* b_mc_truth_h0_id;
  TBranch* b_mc_truth_h0_p4;
  TBranch* b_mc_truth_h0Wl1_id;
  TBranch* b_mc_truth_h0Wl1_p4;
  TBranch* b_mc_truth_h0Wl2_id;
  TBranch* b_mc_truth_h0Wl2_p4;
  //TBranch* b_mc_truth_h0Wnu1_id;
  TBranch* b_mc_truth_h0Wnu1_p4;
  //TBranch* b_mc_truth_h0Wnu2_id;
  TBranch* b_mc_truth_h0Wnu2_p4;
  TBranch* b_mc_truth_h0Wq11_id;
  TBranch* b_mc_truth_h0Wq11_p4;
  //TBranch* b_mc_truth_h0Wq12_id;
  TBranch* b_mc_truth_h0Wq12_p4;
  TBranch* b_mc_truth_h0Wq21_id;
  TBranch* b_mc_truth_h0Wq21_p4;
  //TBranch* b_mc_truth_h0Wq22_id;
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
  //TBranch* b_mc_truth_tWnu1_id;
  TBranch* b_mc_truth_tWnu1_p4;
  //TBranch* b_mc_truth_tWnu2_id;
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


  t->SetBranchAddress("mc_weight",&mc_weight,&b_mc_weight);
  t->SetBranchAddress("mc_truth_h0_id",&mc_truth_h0_id,&b_mc_truth_h0_id);
  t->SetBranchAddress("mc_truth_h0_p4",&mc_truth_h0_p4,&b_mc_truth_h0_p4);
  t->SetBranchAddress("mc_truth_h0Wl1_id",&mc_truth_h0Wl1_id,&b_mc_truth_h0Wl1_id);
  t->SetBranchAddress("mc_truth_h0Wl1_p4",&mc_truth_h0Wl1_p4,&b_mc_truth_h0Wl1_p4);
  t->SetBranchAddress("mc_truth_h0Wl2_id",&mc_truth_h0Wl2_id,&b_mc_truth_h0Wl2_id);
  t->SetBranchAddress("mc_truth_h0Wl2_p4",&mc_truth_h0Wl2_p4,&b_mc_truth_h0Wl2_p4);
  t->SetBranchAddress("mc_truth_h0Wnu1_p4",&mc_truth_h0Wnu1_p4,&b_mc_truth_h0Wnu1_p4);
  t->SetBranchAddress("mc_truth_h0Wnu2_p4",&mc_truth_h0Wnu2_p4,&b_mc_truth_h0Wnu2_p4);
  t->SetBranchAddress("mc_truth_h0Wq11_id",&mc_truth_h0Wq11_id,&b_mc_truth_h0Wq11_id);
  t->SetBranchAddress("mc_truth_h0Wq11_p4",&mc_truth_h0Wq11_p4,&b_mc_truth_h0Wq11_p4);
  t->SetBranchAddress("mc_truth_h0Wq12_p4",&mc_truth_h0Wq12_p4,&b_mc_truth_h0Wq12_p4);
  t->SetBranchAddress("mc_truth_h0Wq21_id",&mc_truth_h0Wq21_id,&b_mc_truth_h0Wq21_id);
  t->SetBranchAddress("mc_truth_h0Wq21_p4",&mc_truth_h0Wq21_p4,&b_mc_truth_h0Wq21_p4);
  t->SetBranchAddress("mc_truth_h0Wq22_p4",&mc_truth_h0Wq22_p4,&b_mc_truth_h0Wq22_p4);
  t->SetBranchAddress("mc_truth_t1_id",&mc_truth_t1_id,&b_mc_truth_t1_id);
  t->SetBranchAddress("mc_truth_t1_p4",&mc_truth_t1_p4,&b_mc_truth_t1_p4);
  t->SetBranchAddress("mc_truth_t2_id",&mc_truth_t2_id,&b_mc_truth_t2_id);
  t->SetBranchAddress("mc_truth_t2_p4",&mc_truth_t2_p4,&b_mc_truth_t2_p4);
  t->SetBranchAddress("mc_truth_tb1_id",&mc_truth_tb1_id,&b_mc_truth_tb1_id);
  t->SetBranchAddress("mc_truth_tb1_p4",&mc_truth_tb1_p4,&b_mc_truth_tb1_p4);
  t->SetBranchAddress("mc_truth_tb2_id",&mc_truth_tb2_id,&b_mc_truth_tb2_id);
  t->SetBranchAddress("mc_truth_tb2_p4",&mc_truth_tb2_p4,&b_mc_truth_tb2_p4);
  t->SetBranchAddress("mc_truth_tWl1_id",&mc_truth_tWl1_id,&b_mc_truth_tWl1_id);
  t->SetBranchAddress("mc_truth_tWl1_p4",&mc_truth_tWl1_p4,&b_mc_truth_tWl1_p4);
  t->SetBranchAddress("mc_truth_tWl2_id",&mc_truth_tWl2_id,&b_mc_truth_tWl2_id);
  t->SetBranchAddress("mc_truth_tWl2_p4",&mc_truth_tWl2_p4,&b_mc_truth_tWl2_p4);
  t->SetBranchAddress("mc_truth_tWnu1_p4",&mc_truth_tWnu1_p4,&b_mc_truth_tWnu1_p4);
  t->SetBranchAddress("mc_truth_tWnu2_p4",&mc_truth_tWnu2_p4,&b_mc_truth_tWnu2_p4);
  t->SetBranchAddress("mc_truth_tWq11_id",&mc_truth_tWq11_id,&b_mc_truth_tWq11_id);
  t->SetBranchAddress("mc_truth_tWq11_p4",&mc_truth_tWq11_p4,&b_mc_truth_tWq11_p4);
  t->SetBranchAddress("mc_truth_tWq21_id",&mc_truth_tWq21_id,&b_mc_truth_tWq21_id);
  t->SetBranchAddress("mc_truth_tWq21_p4",&mc_truth_tWq21_p4,&b_mc_truth_tWq21_p4);
  t->SetBranchAddress("mc_truth_tWq12_id",&mc_truth_tWq12_id,&b_mc_truth_tWq12_id);
  t->SetBranchAddress("mc_truth_tWq12_p4",&mc_truth_tWq12_p4,&b_mc_truth_tWq12_p4);
  t->SetBranchAddress("mc_truth_tWq22_id",&mc_truth_tWq22_id,&b_mc_truth_tWq22_id);
  t->SetBranchAddress("mc_truth_tWq22_p4",&mc_truth_tWq22_p4,&b_mc_truth_tWq22_p4);
  t->SetBranchAddress("mc_truth_Z_id",&mc_truth_Z_id,&b_mc_truth_Z_id);
  t->SetBranchAddress("mc_truth_Z_p4",&mc_truth_Z_p4,&b_mc_truth_Z_p4);
  t->SetBranchAddress("mc_truth_Zl1_id",&mc_truth_Zl1_id,&b_mc_truth_Zl1_id);
  t->SetBranchAddress("mc_truth_Zl1_p4",&mc_truth_Zl1_p4,&b_mc_truth_Zl1_p4);
  t->SetBranchAddress("mc_truth_Zl2_id",&mc_truth_Zl2_id,&b_mc_truth_Zl2_id);
  t->SetBranchAddress("mc_truth_Zl2_p4",&mc_truth_Zl2_p4,&b_mc_truth_Zl2_p4);

  TFile* fOutput = new TFile("output.root", "RECREATE");
  TTree* t2 = new TTree("Tree", "Tree");

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
  Float_t mc_mem_tth_weight;
  Float_t mc_mem_tth_weight_err;
  Float_t mc_mem_tth_weight_chi2;
  Float_t mc_mem_tth_weight_time;
  Float_t mc_mem_tth_weight2;
  Float_t mc_mem_tth_weight2_err;
  Float_t mc_mem_tth_weight2_chi2;
  Float_t mc_mem_tth_weight2_time;
  Float_t mc_mem_ttz_weight;
  Float_t mc_mem_ttz_weight_err;
  Float_t mc_mem_ttz_weight_chi2;
  Float_t mc_mem_ttz_weight_time;
  Float_t mc_mem_ttz_weight2;
  Float_t mc_mem_ttz_weight2_err;
  Float_t mc_mem_ttz_weight2_chi2;
  Float_t mc_mem_ttz_weight2_time;
  Float_t mc_mem_likelihood;
  Float_t mc_mem_likelihood2;

  t2->Branch("mc_weight",&mc_weight,"mc_weight/F");
  t2->Branch("mc_totp4_px",&mc_totp4_px,"mc_totp4_px/F"); 
  t2->Branch("mc_totp4_py",&mc_totp4_py,"mc_totp4_py/F");
  t2->Branch("mc_totp4_pt",&mc_totp4_pt,"mc_totp4_pt/F");
  t2->Branch("mc_thad_pt",&mc_thad_pt,"mc_thad_pt/F");
  t2->Branch("mc_thad_b_pt",&mc_thad_b_pt,"mc_thad_b_pt/F");
  t2->Branch("mc_thad_b_eta",&mc_thad_b_eta,"mc_thad_b_eta/F");
  t2->Branch("mc_thad_j1_pt",&mc_thad_j1_pt,"mc_thad_j1_pt/F");
  t2->Branch("mc_thad_j1_eta",&mc_thad_j1_eta,"mc_thad_j1_eta/F");
  t2->Branch("mc_thad_j2_pt",&mc_thad_j2_pt,"mc_thad_j2_pt/F");
  t2->Branch("mc_thad_j2_eta",&mc_thad_j2_eta,"mc_thad_j2_eta/F");
  t2->Branch("mc_tlep_pt",&mc_tlep_pt,"mc_tlep_pt/F");
  t2->Branch("mc_tlep_b_pt",&mc_tlep_b_pt,"mc_tlep_b_pt/F");
  t2->Branch("mc_tlep_b_eta",&mc_tlep_b_eta,"mc_tlep_b_eta/F");
  t2->Branch("mc_tlep_l_pt",&mc_tlep_l_pt,"mc_tlep_l_pt/F");
  t2->Branch("mc_tlep_l_eta",&mc_tlep_l_eta,"mc_tlep_l_eta/F");
  t2->Branch("mc_tlep2_pt",&mc_tlep2_pt,"mc_tlep2_pt/F");
  t2->Branch("mc_tlep2_b_pt",&mc_tlep2_b_pt,"mc_tlep2_b_pt/F");
  t2->Branch("mc_tlep2_b_eta",&mc_tlep2_b_eta,"mc_tlep2_b_eta/F");
  t2->Branch("mc_tlep2_l_pt",&mc_tlep2_l_pt,"mc_tlep2_l_pt/F");
  t2->Branch("mc_tlep2_l_eta",&mc_tlep2_l_eta,"mc_tlep2_l_eta/F");
  t2->Branch("mc_boson_decay",&mc_boson_decay,"mc_boson_decay/I");
  t2->Branch("mc_boson_pt",&mc_boson_pt,"mc_boson_pt/F");
  t2->Branch("mc_boson_l1_pt",&mc_boson_l1_pt,"mc_boson_l1_pt/F");
  t2->Branch("mc_boson_l1_eta",&mc_boson_l1_eta,"mc_boson_l1_eta/F");
  t2->Branch("mc_boson_l2_pt",&mc_boson_l2_pt,"mc_boson_l2_pt/F");
  t2->Branch("mc_boson_l2_eta",&mc_boson_l2_eta,"mc_boson_l2_eta/F");
  t2->Branch("mc_boson_ll_mass",&mc_boson_ll_mass,"mc_boson_ll_mass/F");
  t2->Branch("mc_boson_ll_pt",&mc_boson_ll_pt,"mc_boson_ll_pt/F");
  t2->Branch("mc_boson_ll_dphi",&mc_boson_ll_dphi,"mc_boson_ll_dphi/F");
  t2->Branch("mc_boson_j1_pt",&mc_boson_j1_pt,"mc_boson_j1_pt/F");
  t2->Branch("mc_boson_j1_eta",&mc_boson_j1_eta,"mc_boson_j1_eta/F");
  t2->Branch("mc_boson_j2_pt",&mc_boson_j2_pt,"mc_boson_j2_pt/F");
  t2->Branch("mc_boson_j2_eta",&mc_boson_j2_eta,"mc_boson_j2_eta/F");
  t2->Branch("mc_boson_jj_mass",&mc_boson_jj_mass,"mc_boson_jj_mass/F");
  t2->Branch("mc_boson_jj_pt",&mc_boson_jj_pt,"mc_boson_jj_pt/F");
  t2->Branch("mc_boson_jj_dphi",&mc_boson_jj_dphi,"mc_boson_jj_dphi/F");
  t2->Branch("mc_met",&mc_met,"mc_met/F");
  t2->Branch("mc_mem_tth_weight",&mc_mem_tth_weight,"mc_mem_tth_weight/F");
  t2->Branch("mc_mem_tth_weight_err",&mc_mem_tth_weight_err,"mc_mem_tth_weight_err/F");
  t2->Branch("mc_mem_tth_weight_chi2",&mc_mem_tth_weight_chi2,"mc_mem_tth_weight_chi2/F");
  t2->Branch("mc_mem_tth_weight_time",&mc_mem_tth_weight_time,"mc_mem_tth_weight_time/F");
  t2->Branch("mc_mem_tth_weight2",&mc_mem_tth_weight2,"mc_mem_tth_weight2/F");
  t2->Branch("mc_mem_tth_weight2_err",&mc_mem_tth_weight2_err,"mc_mem_tth_weight2_err/F");
  t2->Branch("mc_mem_tth_weight2_chi2",&mc_mem_tth_weight2_chi2,"mc_mem_tth_weight2_chi2/F");
  t2->Branch("mc_mem_tth_weight2_time",&mc_mem_tth_weight2_time,"mc_mem_tth_weight2_time/F");
  t2->Branch("mc_mem_ttz_weight",&mc_mem_ttz_weight,"mc_mem_ttz_weight/F");
  t2->Branch("mc_mem_ttz_weight_err",&mc_mem_ttz_weight_err,"mc_mem_ttz_weight_err/F");
  t2->Branch("mc_mem_ttz_weight_chi2",&mc_mem_ttz_weight_chi2,"mc_mem_ttz_weight_chi2/F");
  t2->Branch("mc_mem_ttz_weight_time",&mc_mem_ttz_weight_time,"mc_mem_ttz_weight_time/F");
  t2->Branch("mc_mem_ttz_weight2",&mc_mem_ttz_weight2,"mc_mem_ttz_weight2/F");
  t2->Branch("mc_mem_ttz_weight2_err",&mc_mem_ttz_weight2_err,"mc_mem_ttz_weight2_err/F");
  t2->Branch("mc_mem_ttz_weight2_chi2",&mc_mem_ttz_weight2_chi2,"mc_mem_ttz_weight2_chi2/F");
  t2->Branch("mc_mem_ttz_weight2_time",&mc_mem_ttz_weight2_time,"mc_mem_ttz_weight2_time/F");
  t2->Branch("mc_mem_likelihood",&mc_mem_likelihood,"mc_mem_likelihood/F");
  t2->Branch("mc_mem_likelihood2",&mc_mem_likelihood2,"mc_mem_likelihood2/F");

  int nEvents = t->GetEntries();
  //double *wTTH = new double[nEvents];

  cout << "Start running on events"<<endl;
  double click1=0;
  double res=0;
  double time=0;
  double pErr=0;
  double chi2=0;

  for (Long64_t iEvent=atoi(argv[2]); iEvent<atoi(argv[3]); iEvent++){
    t->LoadTree(iEvent);
    t->GetEntry(iEvent);
    //cout << "iEvent="<<iEvent<<" mc_weigth="<<mc_weight<<endl;

    if (mc_truth_t1_id==-666 || mc_truth_t2_id==-666) continue; 

    //Reconstruct Ptot (needed for recoil)
    TLorentzVector Ptot, Pll, Pl1, PtotNeut(0,0,0,0), Pjj; 

    if (mc_truth_h0_id>-666 && mc_truth_h0Wl1_id>-666 && mc_truth_h0Wl2_id>-666){
      //Fully leptonic Higgs
      hypIntegrator.meIntegrator->MEMFix_HiggsFullLep.Lep1_E = mc_truth_h0Wl1_p4->E();
      hypIntegrator.meIntegrator->MEMFix_HiggsFullLep.Lep1_Theta = mc_truth_h0Wl1_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_HiggsFullLep.Lep1_Phi = mc_truth_h0Wl1_p4->Phi();
      hypIntegrator.meIntegrator->MEMFix_HiggsFullLep.Lep2_E = mc_truth_h0Wl2_p4->E();
      hypIntegrator.meIntegrator->MEMFix_HiggsFullLep.Lep2_Theta = mc_truth_h0Wl2_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_HiggsFullLep.Lep2_Phi = mc_truth_h0Wl2_p4->Phi();
      //cout << "FullLep Higgs id="<<mc_truth_h0_id<<" Pt="<<mc_truth_h0_p4->Pt() <<endl;
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
      //Semi-leptonic Higgs W1->jets
      hypIntegrator.meIntegrator->MEMFix_HiggsSemiLep.Lep1_E = mc_truth_h0Wl2_p4->E();
      hypIntegrator.meIntegrator->MEMFix_HiggsSemiLep.Lep1_Theta = mc_truth_h0Wl2_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_HiggsSemiLep.Lep1_Phi = mc_truth_h0Wl2_p4->Phi();
      hypIntegrator.meIntegrator->MEMFix_HiggsSemiLep.Jet1_Theta = mc_truth_h0Wq11_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_HiggsSemiLep.Jet1_Phi = mc_truth_h0Wq11_p4->Phi();
      hypIntegrator.meIntegrator->MEMFix_HiggsSemiLep.Jet2_Theta = mc_truth_h0Wq12_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_HiggsSemiLep.Jet2_Phi = mc_truth_h0Wq12_p4->Phi();
      mc_boson_decay = 1;
      Ptot = (*mc_truth_h0_p4) + (*mc_truth_t1_p4) + (*mc_truth_t2_p4);
      Pjj = (*mc_truth_h0Wq11_p4) + (*mc_truth_h0Wq12_p4);
      mc_boson_pt = mc_truth_h0_p4->Pt();
      mc_boson_l1_pt = mc_truth_h0Wl2_p4->Pt();
      mc_boson_l1_eta = mc_truth_h0Wl2_p4->Eta();
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
      PtotNeut += (*mc_truth_h0Wnu2_p4);
      hypIntegrator.meIntegrator->MeasuredVarForTF.Jet_E = mc_truth_h0Wq11_p4->E();
      hypIntegrator.meIntegrator->MeasuredVarForTF.Jet2_E = mc_truth_h0Wq12_p4->E();
   }
   else if (mc_truth_h0_id>-666 && mc_truth_h0Wq21_id>-666 && mc_truth_h0Wl1_id>-666){
     //Semi-leptonic Higgs W2->jets
      hypIntegrator.meIntegrator->MEMFix_HiggsSemiLep.Lep1_E = mc_truth_h0Wl1_p4->E();
      hypIntegrator.meIntegrator->MEMFix_HiggsSemiLep.Lep1_Theta = mc_truth_h0Wl1_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_HiggsSemiLep.Lep1_Phi = mc_truth_h0Wl1_p4->Phi();
      hypIntegrator.meIntegrator->MEMFix_HiggsSemiLep.Jet1_Theta = mc_truth_h0Wq21_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_HiggsSemiLep.Jet1_Phi = mc_truth_h0Wq21_p4->Phi();
      hypIntegrator.meIntegrator->MEMFix_HiggsSemiLep.Jet2_Theta = mc_truth_h0Wq22_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_HiggsSemiLep.Jet2_Phi = mc_truth_h0Wq22_p4->Phi();
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
      hypIntegrator.meIntegrator->MeasuredVarForTF.Jet_E = mc_truth_h0Wq21_p4->E();
      hypIntegrator.meIntegrator->MeasuredVarForTF.Jet2_E = mc_truth_h0Wq22_p4->E();
    }
    else if (mc_truth_Z_id>-666 && mc_truth_Zl1_id>-666 && mc_truth_Zl2_id>-666){
      //Z or Gamma*->ll
      hypIntegrator.meIntegrator->MEMFix_HiggsFullLep.Lep1_E = mc_truth_Zl1_p4->E();
      hypIntegrator.meIntegrator->MEMFix_HiggsFullLep.Lep1_Theta = mc_truth_Zl1_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_HiggsFullLep.Lep1_Phi = mc_truth_Zl1_p4->Phi();
      hypIntegrator.meIntegrator->MEMFix_HiggsFullLep.Lep2_E = mc_truth_Zl2_p4->E();
      hypIntegrator.meIntegrator->MEMFix_HiggsFullLep.Lep2_Theta = mc_truth_Zl2_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_HiggsFullLep.Lep2_Phi = mc_truth_Zl2_p4->Phi();
      mc_boson_decay = 2;
      Ptot = (*mc_truth_Zl1_p4) + (*mc_truth_Zl2_p4) + (*mc_truth_t1_p4) + (*mc_truth_t2_p4);
      mc_boson_pt = mc_truth_Z_p4->Pt();
      mc_boson_l1_pt = mc_truth_Zl1_p4->Pt();
      mc_boson_l1_eta = mc_truth_Zl1_p4->Eta();
      mc_boson_l2_pt = mc_truth_Zl2_p4->Pt();
      mc_boson_l2_eta = mc_truth_Zl2_p4->Eta();
      mc_boson_ll_mass = mc_truth_Z_p4->M();
      mc_boson_ll_pt = mc_truth_Z_p4->Pt();
      mc_boson_ll_dphi = TMath::Abs(mc_truth_Zl1_p4->DeltaPhi(*mc_truth_Zl2_p4));
    }
    else continue;
 
    //TTbar
    if (mc_truth_tWq11_id>-666 && mc_truth_tWl2_id>-666){
      //Top1 had, top2 lep
      hypIntegrator.meIntegrator->MEMFix_TopHad.Bjet_Theta = mc_truth_tb1_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_TopHad.Bjet_Phi = mc_truth_tb1_p4->Phi(); 
      hypIntegrator.meIntegrator->MEMFix_TopHad.Jet1_Theta = mc_truth_tWq11_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_TopHad.Jet1_Phi = mc_truth_tWq11_p4->Phi();
      hypIntegrator.meIntegrator->MEMFix_TopHad.Jet2_Theta = mc_truth_tWq12_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_TopHad.Jet2_Phi = mc_truth_tWq12_p4->Phi();

      hypIntegrator.meIntegrator->MeasuredVarForTF.Jet_E = mc_truth_tWq11_p4->E();

      hypIntegrator.meIntegrator->MEMFix_TopLep.Bjet_Theta = mc_truth_tb2_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_TopLep.Bjet_Phi = mc_truth_tb2_p4->Phi();
      hypIntegrator.meIntegrator->MEMFix_TopLep.Lep_E = mc_truth_tWl2_p4->E();
      hypIntegrator.meIntegrator->MEMFix_TopLep.Lep_Theta = mc_truth_tWl2_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_TopLep.Lep_Phi = mc_truth_tWl2_p4->Phi();

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
      //Top2 had, top1 lep
      hypIntegrator.meIntegrator->MEMFix_TopHad.Bjet_Theta = mc_truth_tb2_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_TopHad.Bjet_Phi = mc_truth_tb2_p4->Phi();
      hypIntegrator.meIntegrator->MEMFix_TopHad.Jet1_Theta = mc_truth_tWq21_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_TopHad.Jet1_Phi = mc_truth_tWq21_p4->Phi();
      hypIntegrator.meIntegrator->MEMFix_TopHad.Jet2_Theta = mc_truth_tWq22_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_TopHad.Jet2_Phi = mc_truth_tWq22_p4->Phi();

      hypIntegrator.meIntegrator->MeasuredVarForTF.Jet_E = mc_truth_tWq21_p4->E();

      hypIntegrator.meIntegrator->MEMFix_TopLep.Bjet_Theta = mc_truth_tb1_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_TopLep.Bjet_Phi = mc_truth_tb1_p4->Phi();
      hypIntegrator.meIntegrator->MEMFix_TopLep.Lep_E = mc_truth_tWl1_p4->E();
      hypIntegrator.meIntegrator->MEMFix_TopLep.Lep_Theta = mc_truth_tWl1_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_TopLep.Lep_Phi = mc_truth_tWl1_p4->Phi();

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
   else if (mc_truth_tWl1_id>-666 && mc_truth_tWl2_id>-666 && mc_boson_decay==2){
      //Top1 lep, top2 lep (cas higgs semi-lep)

      hypIntegrator.meIntegrator->MEMFix_TopLep.Bjet_Theta = mc_truth_tb1_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_TopLep.Bjet_Phi = mc_truth_tb1_p4->Phi();
      hypIntegrator.meIntegrator->MEMFix_TopLep.Lep_E = mc_truth_tWl1_p4->E();
      hypIntegrator.meIntegrator->MEMFix_TopLep.Lep_Theta = mc_truth_tWl1_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_TopLep.Lep_Phi = mc_truth_tWl1_p4->Phi();

      hypIntegrator.meIntegrator->MEMFix_TopLep2.Bjet_Theta = mc_truth_tb2_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_TopLep2.Bjet_Phi = mc_truth_tb2_p4->Phi();
      hypIntegrator.meIntegrator->MEMFix_TopLep2.Lep_E = mc_truth_tWl2_p4->E();
      hypIntegrator.meIntegrator->MEMFix_TopLep2.Lep_Theta = mc_truth_tWl2_p4->Theta();
      hypIntegrator.meIntegrator->MEMFix_TopLep2.Lep_Phi = mc_truth_tWl2_p4->Phi();

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
      mc_tlep2_b_pt = mc_truth_tb2_p4->Pt();
      mc_tlep2_b_eta = mc_truth_tb2_p4->Eta();
      mc_tlep2_l_pt = mc_truth_tWl2_p4->Pt();
      mc_tlep2_l_eta = mc_truth_tWl2_p4->Eta();

      PtotNeut =PtotNeut + (*mc_truth_tWnu1_p4)+ (*mc_truth_tWnu2_p4);

   }
   else continue;

   mc_totp4_px = Ptot.Px();
   mc_totp4_py = Ptot.Py();
   mc_totp4_pt = Ptot.Pt();
   mc_met = PtotNeut.Pt();

   cout << "iEvent="<<iEvent<<" mc_weigth="<<mc_weight<<" PtTot="<< Ptot.Pt()<<endl;

   //Measured values for Transfer functions
   hypIntegrator.meIntegrator->MeasuredVarForTF.Bjet1_E =  mc_truth_tb1_p4->E();
   hypIntegrator.meIntegrator->MeasuredVarForTF.Bjet2_E =  mc_truth_tb2_p4->E();
   hypIntegrator.meIntegrator->MeasuredVarForTF.Recoil_Px = -Ptot.Px();
   hypIntegrator.meIntegrator->MeasuredVarForTF.Recoil_Py = -Ptot.Py();
 
  //TODO: check which one is a top/antitop

  xU[0] = 1.6*hypIntegrator.meIntegrator->MeasuredVarForTF.Bjet1_E; //TopHad, Bjet_E
  xU[1] = 1.6*hypIntegrator.meIntegrator->MeasuredVarForTF.Jet_E; //TopHad, Jet1_E
  xU[2] = 1.6*hypIntegrator.meIntegrator->MeasuredVarForTF.Bjet2_E; //TopLep, Bjet_E

  //Perform integration TTZ hypothesis
  hypIntegrator.SetupIntegrationHypothesis(kMEM_TTLL_TopAntitopDecay, 0, 50000);
  IntegrationResult res = hypIntegrator.DoIntegration(xL, xU);
  mc_mem_ttz_weight = res.weight / xsTTLLsl;
  mc_mem_ttz_weight_time = res.time;
  mc_mem_ttz_weight_err = res.err / xsTTLLsl;
  mc_mem_ttz_weight_chi2 = res.chi2;
  cout << "TTLL Vegas Warm-up: Ncall="<<intPoints <<" Cross section (pb) : " << mc_mem_ttz_weight<< " +/- "<< mc_mem_ttz_weight_err<<" chi2/ndof="<< mc_mem_ttz_weight_chi2<<" Time(s)="<<mc_mem_ttz_weight_time<<endl;

  //Perform integration TTH hypothesis
  hypIntegrator.SetupIntegrationHypothesis(kMEM_TTH_TopAntitopHiggsDecay, 0, 1000000);
  res = hypIntegrator.DoIntegration(xL, xU);
  mc_mem_tth_weight = res.weight / xsTTslHl;
  mc_mem_tth_weight_time = res.time;
  mc_mem_tth_weight_err = res.err / xsTTslHl;
  mc_mem_tth_weight_chi2 = res.chi2;
  cout << "TTH Vegas Warm-up: Ncall="<<intPoints <<" Cross section (pb) : " << mc_mem_tth_weight << " +/- "<< mc_mem_tth_weight_err<<" chi2/ndof="<< mc_mem_tth_weight_chi2<<" Time(s)="<<mc_mem_tth_weight_time<<endl;

  mc_mem_likelihood = mc_mem_ttz_weight/(mc_mem_tth_weight+mc_mem_ttz_weight);

  t2->Fill();

 }

  t2->Write();
  fOutput->Close();

  return 0;
}
