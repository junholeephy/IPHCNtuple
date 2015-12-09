#ifndef MEINTEGRATOR_H
#define MEINTEGRATOR_H

#include <iostream>
#include <iomanip>

#include "TLorentzVector.h"
#include "TMath.h"

#include "LHAPDF/LHAPDF.h"
//#include "/grid_mnt/home/nchanon/LHAPDF-6.1.5-install/include/LHAPDF/LHAPDF.h"
#include "rambo.h"
#include "../Madgraph/PROC_SA_CPP_sm_4/SubProcesses/P0_Sigma_sm_gg_ttxh/CPPProcess.h"
#include "../Madgraph/PROC_SA_CPP_sm_DECAY_tbwjj/SubProcesses/P0_Sigma_sm_t_budx/CPPProcess_tbwjj.h"
#include "../Madgraph/PROC_SA_CPP_sm_DECAY_tbwlnu/SubProcesses/P0_Sigma_sm_t_bepve/CPPProcess_tbwlnu.h"
#include "../Madgraph/PROC_SA_CPP_sm_DECAY_hw2l2nu/SubProcesses/P0_Sigma_sm_h_epveemvex/CPPProcess_hw2l2nu.h"
#include "../Madgraph/PROC_SA_CPP_sm_DECAY_antitbwjj/SubProcesses/P0_Sigma_sm_tx_bxdux/CPPProcess_antitbwjj.h"
#include "../Madgraph/PROC_SA_CPP_sm_DECAY_ggttll/SubProcesses/P0_Sigma_sm_gg_ttxepem/CPPProcess_ggttll.h"
//#include "../Madgraph/PROC_SA_CPP_sm_DECAY_qqttlpvl/SubProcesses/P0_Sigma_sm_udx_ttxepve/CPPProcess_qqttlpvl.h"
//#include "../Madgraph/PROC_SA_CPP_sm_DECAY_qqttlmvl/SubProcesses/P0_Sigma_sm_dux_ttxemvex/CPPProcess_qqttlmvl.h"
#include "../Madgraph/PROC_SA_CPP_sm_ckm_DECAY_qqttlpvl/SubProcesses/P0_Sigma_sm_ckm_cdx_ttxepve/CPPProcess_P0_Sigma_sm_ckm_cdx_ttxepve.h"
#include "../Madgraph/PROC_SA_CPP_sm_ckm_DECAY_qqttlpvl/SubProcesses/P0_Sigma_sm_ckm_udx_ttxepve/CPPProcess_P0_Sigma_sm_ckm_udx_ttxepve.h"
#include "../Madgraph/PROC_SA_CPP_sm_ckm_DECAY_qqttlpvl/SubProcesses/P0_Sigma_sm_ckm_usx_ttxepve/CPPProcess_P0_Sigma_sm_ckm_usx_ttxepve.h"
#include "../Madgraph/PROC_SA_CPP_sm_ckm_DECAY_qqttlmvl/SubProcesses/P0_Sigma_sm_ckm_dcx_ttxemvex/CPPProcess_P0_Sigma_sm_ckm_dcx_ttxemvex.h"
#include "../Madgraph/PROC_SA_CPP_sm_ckm_DECAY_qqttlmvl/SubProcesses/P0_Sigma_sm_ckm_dux_ttxemvex/CPPProcess_P0_Sigma_sm_ckm_dux_ttxemvex.h"
#include "../Madgraph/PROC_SA_CPP_sm_ckm_DECAY_qqttlmvl/SubProcesses/P0_Sigma_sm_ckm_sux_ttxemvex/CPPProcess_P0_Sigma_sm_ckm_sux_ttxemvex.h"
#include "../Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlpvljj/SubProcesses/P0_Sigma_sm_ckm_gc_ttxepvegd/CPPProcess_P0_Sigma_sm_ckm_gc_ttxepvegd.h"
#include "../Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlpvljj/SubProcesses/P0_Sigma_sm_ckm_gdx_ttxepvegcx/CPPProcess_P0_Sigma_sm_ckm_gdx_ttxepvegcx.h"
#include "../Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlpvljj/SubProcesses/P0_Sigma_sm_ckm_gdx_ttxepvegux/CPPProcess_P0_Sigma_sm_ckm_gdx_ttxepvegux.h"
#include "../Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlpvljj/SubProcesses/P0_Sigma_sm_ckm_gsx_ttxepvegux/CPPProcess_P0_Sigma_sm_ckm_gsx_ttxepvegux.h"
#include "../Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlpvljj/SubProcesses/P0_Sigma_sm_ckm_gu_ttxepvegd/CPPProcess_P0_Sigma_sm_ckm_gu_ttxepvegd.h"
#include "../Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlpvljj/SubProcesses/P0_Sigma_sm_ckm_gu_ttxepvegs/CPPProcess_P0_Sigma_sm_ckm_gu_ttxepvegs.h"
#include "../Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlmvljj/SubProcesses/P0_Sigma_sm_ckm_gcx_ttxemvexgdx/CPPProcess_P0_Sigma_sm_ckm_gcx_ttxemvexgdx.h"
#include "../Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlmvljj/SubProcesses/P0_Sigma_sm_ckm_gd_ttxemvexgc/CPPProcess_P0_Sigma_sm_ckm_gd_ttxemvexgc.h"
#include "../Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlmvljj/SubProcesses/P0_Sigma_sm_ckm_gd_ttxemvexgu/CPPProcess_P0_Sigma_sm_ckm_gd_ttxemvexgu.h"
#include "../Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlmvljj/SubProcesses/P0_Sigma_sm_ckm_gs_ttxemvexgu/CPPProcess_P0_Sigma_sm_ckm_gs_ttxemvexgu.h"
#include "../Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlmvljj/SubProcesses/P0_Sigma_sm_ckm_gux_ttxemvexgdx/CPPProcess_P0_Sigma_sm_ckm_gux_ttxemvexgdx.h"
#include "../Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlmvljj/SubProcesses/P0_Sigma_sm_ckm_gux_ttxemvexgsx/CPPProcess_P0_Sigma_sm_ckm_gux_ttxemvexgsx.h"


#define kNone -1 //No Integration, just evaluation
#define kInitialPartons 0 //Integration over bjorken x, given the final states
#define kAllPartonsTTH 1 //Integration over everything (total cross-section)
//#define kAllPartonsTTH_FSonly 2 //Integration over final state phase space only
#define kAllPartonsTopHad 3
#define kAllPartonsTopHad_FiniteWidth 4
#define kAllPartonsTopLep 5
#define kAllPartonsTopLep_FixedTopM_FiniteWwidth 6
#define kAllPartonsTopHad_FixedTopM_FiniteWwidth 7
#define kAllPartonsHiggsWWLep_FixedHiggsM_FiniteWwidth 8
#define kAllPartonsAntiTopHad_FixedTopM_FiniteWwidth 9
#define kAllPartonsTTH_TopDecay 10
#define kAllPartonsTTLL 11

#define kMEM_TTH_TopAntitopHiggsDecay 12
#define kMEM_TTLL_TopAntitopDecay 13
#define kMEM_TTH_TopAntitopHiggsSemiLepDecay 14
#define kMEM_TTW_TopAntitopDecay 15
#define kMEM_TTWJJ_TopAntitopDecay 16

#define kFixMw 1
#define kFixBenergy 2
#define kTwoBjorken 3
#define kOneBjorken 4
#define kTwoBjorkenSubMasses 5
#define kTwoBjorken1to3 6
#define kNoBjorken 7
#define kNoBjorkenMomentum 8

#define kMadgraph 0
#define kGosam 1

#define kTTH 0
#define kTTLL 1
#define kTTW 2
#define kTTWJJ 3
#define kTop 4
#define kAntitop 5

#define kTopLepDecay 0
#define kTopHadDecay 1

#define kTFNone 0
#define kTFGaussian 1
#define kTFHistoBnonB_GausRecoil 2
#define kTFHistoBnonBmET 3

#define kTFRecoilPtot 0
#define kTFRecoilmET 1

//#define kErr_PS_Top 0
//#define kErr_PS_Boson 1
#define kErr_PS_Product 0
//#define kErr_TF_OutOfRange 3
#define kErr_TF_Product 1
#define kErr_ME 2
#define kErr_PDF 3
#define kErr_Weight_Product 4

using namespace std;
using namespace LHAPDF;

//GoSam
//extern "C" void OLP_Start(const char * filename, int* success);
//extern "C" void OLP_EvalSubProcess(int,double*,double,double*,double*);

const string EtaRangeLabel[] = {"00eta08", "08eta16", "16eta24"};
const string EnergyRangeLabel[] = {"25E50", "50E80", "80E120", "120E200", "200E300", "EGT300"};
const float EtaRange[] = {0.0, 0.8, 1.6, 2.4};
const float EnergyRange[] = {25, 50, 80, 120, 200, 300, 7000};
const string MetRangeLabel[] = {"0E100","100E150","EGT150"};
const float MetRange[] = {0, 100, 150, 7000};

class MEPhaseSpace 
{

  public:
    MEPhaseSpace();
    ~MEPhaseSpace();

    double Eval(const double* ) const;

    CPPProcess* process;
    CPPProcess_tbwjj* process_tbwjj;
    CPPProcess_tbwlnu* process_tbwlnu;
    CPPProcess_hw2l2nu* process_hw2l2nu;
    CPPProcess_antitbwjj* process_antitbwjj;
    CPPProcess_ggttll* process_ggttll;
    //CPPProcess_qqttlpvl* process_qqttlpvl;
    //CPPProcess_qqttlmvl* process_qqttlmvl;
    CPPProcess_P0_Sigma_sm_ckm_cdx_ttxepve* process_qqttlpvl_cdx;
    CPPProcess_P0_Sigma_sm_ckm_udx_ttxepve* process_qqttlpvl_udx;
    CPPProcess_P0_Sigma_sm_ckm_usx_ttxepve* process_qqttlpvl_usx;
    CPPProcess_P0_Sigma_sm_ckm_dcx_ttxemvex* process_qqttlmvl_dcx;
    CPPProcess_P0_Sigma_sm_ckm_dux_ttxemvex* process_qqttlmvl_dux;
    CPPProcess_P0_Sigma_sm_ckm_sux_ttxemvex* process_qqttlmvl_sux;
    CPPProcess_P0_Sigma_sm_ckm_gc_ttxepvegd* process_P0_Sigma_sm_ckm_gc_ttxepvegd;
    CPPProcess_P0_Sigma_sm_ckm_gdx_ttxepvegcx* process_P0_Sigma_sm_ckm_gdx_ttxepvegcx;
    CPPProcess_P0_Sigma_sm_ckm_gdx_ttxepvegux* process_P0_Sigma_sm_ckm_gdx_ttxepvegux;
    CPPProcess_P0_Sigma_sm_ckm_gsx_ttxepvegux* process_P0_Sigma_sm_ckm_gsx_ttxepvegux;
    CPPProcess_P0_Sigma_sm_ckm_gu_ttxepvegd* process_P0_Sigma_sm_ckm_gu_ttxepvegd;
    CPPProcess_P0_Sigma_sm_ckm_gu_ttxepvegs* process_P0_Sigma_sm_ckm_gu_ttxepvegs;
    CPPProcess_P0_Sigma_sm_ckm_gcx_ttxemvexgdx* process_P0_Sigma_sm_ckm_gcx_ttxemvexgdx;
    CPPProcess_P0_Sigma_sm_ckm_gd_ttxemvexgc* process_P0_Sigma_sm_ckm_gd_ttxemvexgc;
    CPPProcess_P0_Sigma_sm_ckm_gd_ttxemvexgu* process_P0_Sigma_sm_ckm_gd_ttxemvexgu;
    CPPProcess_P0_Sigma_sm_ckm_gs_ttxemvexgu* process_P0_Sigma_sm_ckm_gs_ttxemvexgu;
    CPPProcess_P0_Sigma_sm_ckm_gux_ttxemvexgdx* process_P0_Sigma_sm_ckm_gux_ttxemvexgdx;
    CPPProcess_P0_Sigma_sm_ckm_gux_ttxemvexgsx* process_P0_Sigma_sm_ckm_gux_ttxemvexgsx;


    const LHAPDF::PDF* pdf;

    double GeV2barn;              // Conversion factor 1/GeV^2 ->pb

    int iMode;
    int iOption;
    int iGen;
    int iCore;
    int nExternals;
    mutable int iCall;
    mutable int iIteration;
    int iTF;
    int iTFOption;

    mutable int isTopAntitop;

    double comEnergy;
    mutable double muF;
    double mTop;
    double mHiggs;
    double mW;
    double mZ;
    double mB;
    double gammaTop;
    double gammaHiggs;
    double gammaW;

    double xsTTH;
    double xsTTLL;
    double xsTTW;
    double brTopLep;
    double brTopHad;
    double brHiggsFullLep;

    mutable int errorCounter[20];

    int CheckMomentum(TLorentzVector&, double) const;

    mutable vector<double*> * pCore; //phase-space point
    mutable vector<double*> * pTop; //phase-space point
    mutable vector<double*> * pAntitop; //phase-space point
    mutable vector<double*> * pHiggs; //phase-space point

    void InitializePhaseSpacePoint(vector<double*> **, int);
    void FillTTHPhaseSpacePoint(TLorentzVector&, TLorentzVector&, TLorentzVector&) const;
    void FillTTLLPhaseSpacePoint(TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector&) const;
    void FillTTLNuJJPhaseSpacePoint(TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector&) const;
    void FillTopDecayPhaseSpacePoint(TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector&, int) const;
    void FillHiggsDecayPhaseSpacePoint(TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector&) const;


    double Wenergy;
    double Benergy;

    void SetComEnergy(double);
    void SetInitialPartonMomenta(const double, const double) const;
    void SetFinalStatePartonMomenta(double**, int) const;
    void SetAllPartonMomenta(vector<double*> *, vector<double*>);
    void ReadPartonMomenta(vector<double*> *, int ) const;
    vector<double*> GetPhaseSpacePoint();

    void SetIntegrationMode(int );
    void SetOption(int );
    void SetGenerator(int);
    void SetTFChoice(int);
    void SetTFOption(int);

    double SetupKinematicsTTH(const double *) const;
    double SetupKinematicsTopHad(const double *) const;
    double SetupKinematicsTopHad_FiniteWidths(const double *) const;
    double SetupKinematicsTopLep(const double *) const; 
    double SetupKinematicsTopLep_FixedTopM_FiniteWwidth(const double *) const; 
    double SetupKinematicsTopHad_FixedTopM_FiniteWwidth(const double *) const;
    double SetupKinematicsHiggsWWLep_FixedHiggsM_FiniteWwidth(const double *) const;
    double SetupKinematics1to3_LabFrame(const double *) const;
    double SetupKinematics2to3_LabFrame_OneBjorken(const double *) const;
    double SetupKinematics2to3_LabFrame_NoBjorken(const double *) const;
    double SetupKinematicsTTH_NoBjorken_TopHadDecay(const double*) const;
    double SetupKinematics_TopHadDecay_WithTopPhaseSpace(const double*, TLorentzVector*, int, int) const;
    double SetupKinematics_Higgs2l2nuDecay_WithHiggsPhaseSpace(const double*, TLorentzVector*) const;
    double SetupKinematics2to4_LabFrame_NoBjorken(const double *) const;
    double SetupKinematicsTTH_NoBjorken_TopAntitopHiggsDecay(const double*) const;
    double SetupKinematicsTTLL_NoBjorken_TopAntitopDecay(const double*) const;
    double SetupKinematicsTTW_NoBjorken_TopAntitopDecay(const double*) const;
    double SetupKinematicsTTWJJ_NoBjorken_TopAntitopDecay(const double*) const;

    void ApplyTotalTransverseBoost() const;

    void CheckMatrixElement();
    double ComputeMatrixElement() const;
    double ComputeSubMatrixElement(int , int , int ) const;	
    double BreitWigner(double, double, double) const;
    double KallenFunction(double , double , double ) const;
    double ComputeDecayMomenta(TLorentzVector& , double, double, double , double , TLorentzVector* , TLorentzVector* ) const;
    int SetMomentumFromEThetaPhi(double, double, double, double, TLorentzVector*) const;

    double ComputePDF(double, double, double) const;
    double ConvolvePdfCrossSection(double, double, double) const;

    double ComputeTFProduct(double, double, double, double, double, double) const;
    double ComputeSingleTF_Gaus(double, double, double) const; 
    double ComputeSingleTF_Histo(string, double, double, double) const;
    void LoadTFfromHisto(string) const;

    TH1F*** TF_nonB;
    TH1F*** TF_B;
    TH1F** TF_MetPx;
    TH1F** TF_MetPy;

    struct MeasuredVarForTF {
      double Bjet1_E;
      double Bjet1_Eta;
      double Bjet2_E;
      double Bjet2_Eta;
      double Jet1_E;
      double Jet1_Eta;
      double Jet2_E;
      double Jet2_Eta;
      double Recoil_Px;
      double Recoil_Py;
      double mET_Px;
      double mET_Py;
    } MeasuredVarForTF;
 
    mutable struct ComputedVarForTF {
      double Bjet1_E;
      double Bjet2_E;
      double Jet1_E;
      double Jet2_E;
      double Recoil_Px;
      double Recoil_Py;
      double mET_Px;
      double mET_Py;
    } ComputedVarForTF;

    mutable TLorentzVector Computed_mETvect;
 
    double* y;
    double* xMEM;

    struct MEMFix_TopHad {
      double Bjet_Theta;
      double Bjet_Phi;
      double Jet1_Theta;
      double Jet1_Phi;
      double Jet2_Theta;
      double Jet2_Phi;
      int TopSign;
    } MEMFix_TopHad;

    struct MEMFix_TopLep {
      double Bjet_Theta;
      double Bjet_Phi;
      double Lep_E;
      double Lep_Theta;
      double Lep_Phi;
      int TopSign;
    } MEMFix_TopLep, MEMFix_TopLep2;

    struct MEMFix_HiggsFullLep {
      double Lep1_E;
      double Lep1_Theta;
      double Lep1_Phi;
      double Lep2_E;
      double Lep2_Theta;
      double Lep2_Phi;
    } MEMFix_HiggsFullLep;

    struct MEMFix_HiggsSemiLep {
      double Lep1_E;
      double Lep1_Theta;
      double Lep1_Phi;
      double Jet1_Theta;
      double Jet1_Phi;
      double Jet2_Theta;
      double Jet2_Phi;
      int LepSign;
    } MEMFix_HiggsSemiLep;


/*
    double* integrationBounds_Down;
    double* integrationBounds_Up; 
    ROOT::Math::Functor* toIntegrate;
    ROOT::Math::GSLMCIntegrator* ig;
    double* y;
    mutable int iSuperCall;
    ROOT::Math::VegasParameters* param;
    mutable ROOT::Math::IntegratorMultiDimOptions opts;
    int intPoints;
*/
    int verbosity;
    void SetVerbosity(int);
};

MEPhaseSpace::MEPhaseSpace(){

  if (verbosity>=1) cout << "MEPhaseSpace constructor" << endl;

  //Initialise all the partons
  pCore = new vector<double*>();
  pTop = new vector<double*>();
  pAntitop = new vector<double*>();
  pHiggs = new vector<double*>();
    
    process = new CPPProcess();
    process->initProc("/afs/cern.ch/work/c/chanon/MEM/Madgraph/PROC_SA_CPP_sm_4/Cards/param_card.dat");
    if (verbosity>=1) cout << "TTH Process nexternal="<<process->nexternal<<endl;
/*
    int check=0;
    std::string fname("../gosam/ttH_ME/OLE_order.olc");
    OLP_Start(fname.c_str(),&check);
    if (verbosity>=1) cout << "GoSam initialization: "<<check<<endl;
*/
  process_ggttll = new CPPProcess_ggttll();
  process_ggttll->initProc("/afs/cern.ch/work/c/chanon/MEM/Madgraph/PROC_SA_CPP_sm_DECAY_ggttll/Cards/param_card.dat");
  if (verbosity>=1) cout << "TTLL Process nexternal="<<process_ggttll->nexternal<<endl;


  process_qqttlpvl_cdx = new CPPProcess_P0_Sigma_sm_ckm_cdx_ttxepve();
  process_qqttlpvl_cdx->initProc("/afs/cern.ch/work/c/chanon/MEM/Madgraph/PROC_SA_CPP_sm_ckm_DECAY_qqttlpvl/Cards/param_card.dat");
  if (verbosity>=1) cout << "TTWP cdbar Process nexternal="<<process_qqttlpvl_cdx->nexternal<<endl;

  process_qqttlpvl_udx = new CPPProcess_P0_Sigma_sm_ckm_udx_ttxepve();
  process_qqttlpvl_udx->initProc("/afs/cern.ch/work/c/chanon/MEM/Madgraph/PROC_SA_CPP_sm_ckm_DECAY_qqttlpvl/Cards/param_card.dat");
  if (verbosity>=1) cout << "TTWP udbar Process nexternal="<<process_qqttlpvl_udx->nexternal<<endl;

  process_qqttlpvl_usx = new CPPProcess_P0_Sigma_sm_ckm_usx_ttxepve();
  process_qqttlpvl_usx->initProc("/afs/cern.ch/work/c/chanon/MEM/Madgraph/PROC_SA_CPP_sm_ckm_DECAY_qqttlpvl/Cards/param_card.dat");
  if (verbosity>=1) cout << "TTWP usbar Process nexternal="<<process_qqttlpvl_usx->nexternal<<endl;

  process_qqttlmvl_dcx = new CPPProcess_P0_Sigma_sm_ckm_dcx_ttxemvex();
  process_qqttlmvl_dcx->initProc("/afs/cern.ch/work/c/chanon/MEM/Madgraph/PROC_SA_CPP_sm_ckm_DECAY_qqttlmvl/Cards/param_card.dat");
  if (verbosity>=1) cout << "TTWM dcbar Process nexternal="<<process_qqttlmvl_dcx->nexternal<<endl;

  process_qqttlmvl_dux = new CPPProcess_P0_Sigma_sm_ckm_dux_ttxemvex();
  process_qqttlmvl_dux->initProc("/afs/cern.ch/work/c/chanon/MEM/Madgraph/PROC_SA_CPP_sm_ckm_DECAY_qqttlmvl/Cards/param_card.dat");
  if (verbosity>=1) cout << "TTWM dubar Process nexternal="<<process_qqttlmvl_dux->nexternal<<endl;

  process_qqttlmvl_sux = new CPPProcess_P0_Sigma_sm_ckm_sux_ttxemvex();
  process_qqttlmvl_sux->initProc("/afs/cern.ch/work/c/chanon/MEM/Madgraph/PROC_SA_CPP_sm_ckm_DECAY_qqttlmvl/Cards/param_card.dat");
  if (verbosity>=1) cout << "TTWM subar Process nexternal="<<process_qqttlmvl_sux->nexternal<<endl;

process_P0_Sigma_sm_ckm_gc_ttxepvegd = new CPPProcess_P0_Sigma_sm_ckm_gc_ttxepvegd();
process_P0_Sigma_sm_ckm_gc_ttxepvegd->initProc("/afs/cern.ch/work/c/chanon/MEM/Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlpvljj/Cards/param_card.dat");
if (verbosity>=1) cout << "TTWPJJ Process nexternal=" << process_P0_Sigma_sm_ckm_gc_ttxepvegd->nexternal << endl;

process_P0_Sigma_sm_ckm_gdx_ttxepvegcx = new CPPProcess_P0_Sigma_sm_ckm_gdx_ttxepvegcx();
process_P0_Sigma_sm_ckm_gdx_ttxepvegcx->initProc("/afs/cern.ch/work/c/chanon/MEM/Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlpvljj/Cards/param_card.dat");
if (verbosity>=1) cout << "TTWPJJ Process nexternal=" << process_P0_Sigma_sm_ckm_gdx_ttxepvegcx->nexternal << endl;

process_P0_Sigma_sm_ckm_gdx_ttxepvegux = new CPPProcess_P0_Sigma_sm_ckm_gdx_ttxepvegux();
process_P0_Sigma_sm_ckm_gdx_ttxepvegux->initProc("/afs/cern.ch/work/c/chanon/MEM/Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlpvljj/Cards/param_card.dat");
if (verbosity>=1) cout << "TTWPJJ Process nexternal=" << process_P0_Sigma_sm_ckm_gdx_ttxepvegux->nexternal << endl;

process_P0_Sigma_sm_ckm_gsx_ttxepvegux = new CPPProcess_P0_Sigma_sm_ckm_gsx_ttxepvegux();
process_P0_Sigma_sm_ckm_gsx_ttxepvegux->initProc("/afs/cern.ch/work/c/chanon/MEM/Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlpvljj/Cards/param_card.dat");
if (verbosity>=1) cout << "TTWPJJ Process nexternal=" << process_P0_Sigma_sm_ckm_gsx_ttxepvegux->nexternal << endl;

process_P0_Sigma_sm_ckm_gu_ttxepvegd = new CPPProcess_P0_Sigma_sm_ckm_gu_ttxepvegd();
process_P0_Sigma_sm_ckm_gu_ttxepvegd->initProc("/afs/cern.ch/work/c/chanon/MEM/Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlpvljj/Cards/param_card.dat");
if (verbosity>=1) cout << "TTWPJJ Process nexternal=" << process_P0_Sigma_sm_ckm_gu_ttxepvegd->nexternal << endl;

process_P0_Sigma_sm_ckm_gu_ttxepvegs = new CPPProcess_P0_Sigma_sm_ckm_gu_ttxepvegs();
process_P0_Sigma_sm_ckm_gu_ttxepvegs->initProc("/afs/cern.ch/work/c/chanon/MEM/Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlpvljj/Cards/param_card.dat");
if (verbosity>=1) cout << "TTWPJJ Process nexternal=" << process_P0_Sigma_sm_ckm_gu_ttxepvegs->nexternal << endl;

process_P0_Sigma_sm_ckm_gcx_ttxemvexgdx = new CPPProcess_P0_Sigma_sm_ckm_gcx_ttxemvexgdx();
process_P0_Sigma_sm_ckm_gcx_ttxemvexgdx->initProc("/afs/cern.ch/work/c/chanon/MEM/Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlmvljj/Cards/param_card.dat");
if (verbosity>=1) cout << "TTWMJJ Process nexternal=" << process_P0_Sigma_sm_ckm_gcx_ttxemvexgdx->nexternal << endl;

process_P0_Sigma_sm_ckm_gd_ttxemvexgc = new CPPProcess_P0_Sigma_sm_ckm_gd_ttxemvexgc();
process_P0_Sigma_sm_ckm_gd_ttxemvexgc->initProc("/afs/cern.ch/work/c/chanon/MEM/Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlmvljj/Cards/param_card.dat");
if (verbosity>=1) cout << "TTWMJJ Process nexternal=" << process_P0_Sigma_sm_ckm_gd_ttxemvexgc->nexternal << endl;

process_P0_Sigma_sm_ckm_gd_ttxemvexgu = new CPPProcess_P0_Sigma_sm_ckm_gd_ttxemvexgu();
process_P0_Sigma_sm_ckm_gd_ttxemvexgu->initProc("/afs/cern.ch/work/c/chanon/MEM/Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlmvljj/Cards/param_card.dat");
if (verbosity>=1) cout << "TTWMJJ Process nexternal=" << process_P0_Sigma_sm_ckm_gd_ttxemvexgu->nexternal << endl;

process_P0_Sigma_sm_ckm_gs_ttxemvexgu = new CPPProcess_P0_Sigma_sm_ckm_gs_ttxemvexgu();
process_P0_Sigma_sm_ckm_gs_ttxemvexgu->initProc("/afs/cern.ch/work/c/chanon/MEM/Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlmvljj/Cards/param_card.dat");
if (verbosity>=1) cout << "TTWMJJ Process nexternal=" << process_P0_Sigma_sm_ckm_gs_ttxemvexgu->nexternal << endl;

process_P0_Sigma_sm_ckm_gux_ttxemvexgdx = new CPPProcess_P0_Sigma_sm_ckm_gux_ttxemvexgdx();
process_P0_Sigma_sm_ckm_gux_ttxemvexgdx->initProc("/afs/cern.ch/work/c/chanon/MEM/Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlmvljj/Cards/param_card.dat");
if (verbosity>=1) cout << "TTWMJJ Process nexternal=" << process_P0_Sigma_sm_ckm_gux_ttxemvexgdx->nexternal << endl;

process_P0_Sigma_sm_ckm_gux_ttxemvexgsx = new CPPProcess_P0_Sigma_sm_ckm_gux_ttxemvexgsx();
process_P0_Sigma_sm_ckm_gux_ttxemvexgsx->initProc("/afs/cern.ch/work/c/chanon/MEM/Madgraph/PROC_SA_CPP_sm_ckm_DECAY_ppttlmvljj/Cards/param_card.dat");
if (verbosity>=1) cout << "TTWMJJ Process nexternal=" << process_P0_Sigma_sm_ckm_gux_ttxemvexgsx->nexternal << endl;


  InitializePhaseSpacePoint(&pCore, 8);

  process_tbwjj = new CPPProcess_tbwjj();
  process_tbwjj->initProc("/afs/cern.ch/work/c/chanon/MEM/Madgraph/PROC_SA_CPP_sm_DECAY_tbwjj/Cards/param_card.dat");
  if (verbosity>=1) cout << "T->bWjj Process nexternal="<<process_tbwjj->nexternal<<endl;
  InitializePhaseSpacePoint(&pTop, 4);

  process_antitbwjj = new CPPProcess_antitbwjj();
  process_antitbwjj->initProc("/afs/cern.ch/work/c/chanon/MEM/Madgraph/PROC_SA_CPP_sm_DECAY_antitbwjj/Cards/param_card.dat");
  if (verbosity>=1) cout << "Anti T->bWjj Process nexternal="<<process_antitbwjj->nexternal<<endl;
  InitializePhaseSpacePoint(&pAntitop, 4);

  process_hw2l2nu = new CPPProcess_hw2l2nu();
  process_hw2l2nu->initProc("/afs/cern.ch/work/c/chanon/MEM/Madgraph/PROC_SA_CPP_sm_DECAY_hw2l2nu/Cards/param_card.dat");
  if (verbosity>=1) cout << "H->WW->2l2nu Process nexternal="<<process_hw2l2nu->nexternal<<endl;
  InitializePhaseSpacePoint(&pHiggs, 5);

  //Initialise LHAPDF
  pdf = LHAPDF::mkPDF("NNPDF23_lo_as_0119_qed",0);

  //Initialize transfert function
  TF_nonB = new TH1F**[3];
  TF_B = new TH1F**[3];
  for (int iEta=0; iEta<3; iEta++) TF_nonB[iEta] = new TH1F*[6];
  for (int iEta=0; iEta<3; iEta++) TF_B[iEta] = new TH1F*[6];
  TF_MetPx = new TH1F*[3];
  TF_MetPy = new TH1F*[3];
  LoadTFfromHisto("/afs/cern.ch/work/c/chanon/MEM/data/TF_Daniel_05112015.root");

  //set constants from madgraph
  mTop = 173.;
  mHiggs = 125.;
  mW = 80.419;
  mZ = 91.188;
  mB = 4.7;
  gammaTop = 1.491500;
  gammaHiggs = 0.006382339; 
  gammaW = 2.0476;
 
  //total cross sections
  xsTTH = 0.284034;
  xsTTLL = 0.02722653;
  xsTTW = 0.04378 + 0.0217;
  brTopHad = 0.97706 / gammaTop;
  brTopLep = brTopHad/3.;
  brHiggsFullLep = 1.570177e-05 / gammaHiggs;

  //set factorization scale
  muF = (2*mTop+mHiggs)/2.;

  iCall = 0;
  iIteration = 0;
  for (int i=0; i<20; i++) errorCounter[i]=0;

  //GeV-2 to pb translation factor
  GeV2barn = 0.38937966 * pow (10, 9.);

  y = new double[15];
  xMEM = new double[28];

/*
  //Setting double integration
  integrationBounds_Down = new double[5];
  integrationBounds_Up = new double[5];
  y = new double[7];
  intPoints = 20000;
  toIntegrate = new ROOT::Math::Functor(this, &MEPhaseSpace::Eval, 5);
  ig = new ROOT::Math::GSLMCIntegrator( ROOT::Math::IntegrationMultiDim::kVEGAS, 1.e-14, 1.e-7, intPoints);
  ig->SetFunction(*toIntegrate);
  iSuperCall = 0;
  param = new ROOT::Math::VegasParameters( *(ig->ExtraOptions()) );
*/
}

MEPhaseSpace::~MEPhaseSpace(){

  //delete ig;
  //delete toIntegrate;

  return;
}

void MEPhaseSpace::SetOption(int ioption){

  iOption = ioption;
  return;
}

void MEPhaseSpace::SetGenerator(int igen){

  iGen = igen;
  return;
}

void MEPhaseSpace::SetTFChoice(int iChoice){

  iTF = iChoice;
  return;
}

void MEPhaseSpace::SetTFOption(int iOption){

  iTFOption = iOption;
  return;
}


void MEPhaseSpace::SetIntegrationMode(int imode){

    if (verbosity>=1) cout << "SetIntegrationMode "<<imode<<endl;
  iMode = imode;

  if (iMode==kMEM_TTH_TopAntitopHiggsDecay || iMode==kMEM_TTH_TopAntitopHiggsSemiLepDecay) iCore = kTTH;
  if (iMode==kMEM_TTLL_TopAntitopDecay) iCore = kTTLL;
  if (iMode==kMEM_TTW_TopAntitopDecay) iCore = kTTW;
  if (iMode==kMEM_TTWJJ_TopAntitopDecay) iCore = kTTWJJ;

  return;
}

void MEPhaseSpace::SetComEnergy(double energy){

  comEnergy = energy;
  if (verbosity>=1) cout << "comEnergy="<<comEnergy<<endl;

  return;
}

void MEPhaseSpace::InitializePhaseSpacePoint(vector<double*> ** p0, int nPart){

  double** part = new double*[nPart];
  for (int i=0; i<nPart; i++){
    part[i] = new double[4];
    part[i][0] = 0;
    part[i][1] = 0;
    part[i][2] = 0;
    part[i][3] = 0;
  }
  for (int i=0; i<nPart; i++){
    (*p0)->push_back(part[i]);
  }

  return;
}

void MEPhaseSpace::SetInitialPartonMomenta(const double x1, const double x2) const{

    if (verbosity>=2) cout << "Setting initial parton momenta" << endl;

  //x[0], x[1]:  bjorken x of initial partons
  pCore->at(0)[0] = x1 * comEnergy / 2.;
  pCore->at(0)[1] = 0;
  pCore->at(0)[2] = 0;
  pCore->at(0)[3] = x1 * comEnergy / 2.;

  //cout << "Initial parton z "<<p->at(0)[3]<<endl;

  pCore->at(1)[0] = x2 * comEnergy / 2.;
  pCore->at(1)[1] = 0;
  pCore->at(1)[2] = 0;
  pCore->at(1)[3] = - x2 * comEnergy / 2.;//(-1) * (*p)[1][0];

    //if (verbosity>=2) cout << "Initial parton z "<<pCore->at(0)[3]<< " and "<<pCore->at(1)[3]<<endl;

 return;
}

void MEPhaseSpace::SetFinalStatePartonMomenta(double** parton, int nparton) const {

  for (int ipart=0; ipart<nparton; ipart++){
      (*pCore)[ipart+2] = parton[ipart];
  }
  return;
}

void MEPhaseSpace::SetAllPartonMomenta(vector<double*> * p0, vector<double*> part){

  (*p0) = part;  
  return;
}

void MEPhaseSpace::ReadPartonMomenta(vector<double*> * p0, int npart) const{

  for (int i=0; i<npart; i++){

        if (verbosity>=2) cout << "Parton "<<i<<" energy="<<p0->at(i)[0]<<" px="<<p0->at(i)[1]<<" py="<<p0->at(i)[2]<<" pz="<<p0->at(i)[3]<<endl;

  }
  return;
}

int MEPhaseSpace::CheckMomentum(TLorentzVector& P, double mass) const {

   double P2 = P.E()*P.E()-(P.Px()*P.Px()+P.Py()*P.Py()+P.Pz()*P.Pz());

   if (P2<-1.e-3)  {
     if (verbosity>=2) cout << "P2="<<P2<<" negative"<<endl;
     return 0;
   }

   else if (abs(sqrt(TMath::Abs(P2))-mass)<1.e-3){
     if (verbosity>=2)  cout << "|p4|="<<sqrt(TMath::Abs(P2))<<" matches expected mass"<<endl;
   }

   //else if (P2<0) {
   //  if (verbosity>=2) cout << "P2="<<P2<<" negative"<<endl;
   //  return 0;
   //}

   else if (abs(sqrt(P2)-mass)>1.e-3){
     if (verbosity>=2) cout << "|p4|="<< sqrt(P2)<<" but different from expected mass"<<endl;
     return 0;
   }

   else {
     if (verbosity>=2) cout << "|p4|="<< sqrt(P2)<<endl;

   }

  return 1;
}

void MEPhaseSpace::FillTTHPhaseSpacePoint(TLorentzVector& Top, TLorentzVector& Antitop, TLorentzVector& Higgs) const{

  pCore->at(2)[0] = Top.E();
  pCore->at(2)[1] = Top.Px();
  pCore->at(2)[2] = Top.Py();
  pCore->at(2)[3] = Top.Pz();

  pCore->at(3)[0] = Antitop.E();
  pCore->at(3)[1] = Antitop.Px();
  pCore->at(3)[2] = Antitop.Py();
  pCore->at(3)[3] = Antitop.Pz();

  pCore->at(4)[0] = Higgs.E();
  pCore->at(4)[1] = Higgs.Px();
  pCore->at(4)[2] = Higgs.Py();
  pCore->at(4)[3] = Higgs.Pz();

  return;
}

void MEPhaseSpace::FillTTLLPhaseSpacePoint(TLorentzVector& Top, TLorentzVector& Antitop, TLorentzVector& Lep1, TLorentzVector& Lep2) const{

  pCore->at(2)[0] = Top.E();
  pCore->at(2)[1] = Top.Px();
  pCore->at(2)[2] = Top.Py();
  pCore->at(2)[3] = Top.Pz();

  pCore->at(3)[0] = Antitop.E();
  pCore->at(3)[1] = Antitop.Px();
  pCore->at(3)[2] = Antitop.Py();
  pCore->at(3)[3] = Antitop.Pz();

  pCore->at(4)[0] = Lep1.E();
  pCore->at(4)[1] = Lep1.Px();
  pCore->at(4)[2] = Lep1.Py();
  pCore->at(4)[3] = Lep1.Pz();

  pCore->at(5)[0] = Lep2.E();
  pCore->at(5)[1] = Lep2.Px();
  pCore->at(5)[2] = Lep2.Py();
  pCore->at(5)[3] = Lep2.Pz();

  return;
}

void MEPhaseSpace::FillTTLNuJJPhaseSpacePoint(TLorentzVector& Top, TLorentzVector& Antitop, TLorentzVector& Lep, TLorentzVector& Nu, TLorentzVector& Jet1, TLorentzVector& Jet2) const {

  pCore->at(2)[0] = Top.E();
  pCore->at(2)[1] = Top.Px();
  pCore->at(2)[2] = Top.Py();
  pCore->at(2)[3] = Top.Pz();

  pCore->at(3)[0] = Antitop.E();
  pCore->at(3)[1] = Antitop.Px();
  pCore->at(3)[2] = Antitop.Py();
  pCore->at(3)[3] = Antitop.Pz();

  pCore->at(4)[0] = Lep.E();
  pCore->at(4)[1] = Lep.Px();
  pCore->at(4)[2] = Lep.Py();
  pCore->at(4)[3] = Lep.Pz();

  pCore->at(5)[0] = Nu.E();
  pCore->at(5)[1] = Nu.Px();
  pCore->at(5)[2] = Nu.Py();
  pCore->at(5)[3] = Nu.Pz();

  pCore->at(6)[0] = Jet1.E();
  pCore->at(6)[1] = Jet1.Px();
  pCore->at(6)[2] = Jet1.Py();
  pCore->at(6)[3] = Jet1.Pz();

  pCore->at(7)[0] = Jet2.E();
  pCore->at(7)[1] = Jet2.Px();
  pCore->at(7)[2] = Jet2.Py();
  pCore->at(7)[3] = Jet2.Pz();

  return;
}

void MEPhaseSpace::FillTopDecayPhaseSpacePoint(TLorentzVector& Top, TLorentzVector& Bjet, TLorentzVector& Decay1, TLorentzVector& Decay2, int TopType) const{

  if (TopType==kAntitop){
    pAntitop->at(0)[0] = Top.E();
    pAntitop->at(0)[1] = Top.Px();
    pAntitop->at(0)[2] = Top.Py();
    pAntitop->at(0)[3] = Top.Pz();

    pAntitop->at(1)[0] = Bjet.E();
    pAntitop->at(1)[1] = Bjet.Px();
    pAntitop->at(1)[2] = Bjet.Py();
    pAntitop->at(1)[3] = Bjet.Pz();

    pAntitop->at(2)[0] = Decay1.E();
    pAntitop->at(2)[1] = Decay1.Px();
    pAntitop->at(2)[2] = Decay1.Py();
    pAntitop->at(2)[3] = Decay1.Pz();

    pAntitop->at(3)[0] = Decay2.E();
    pAntitop->at(3)[1] = Decay2.Px();
    pAntitop->at(3)[2] = Decay2.Py();
    pAntitop->at(3)[3] = Decay2.Pz();
  }
  else if (TopType==kTop){
    pTop->at(0)[0] = Top.E();
    pTop->at(0)[1] = Top.Px();
    pTop->at(0)[2] = Top.Py();
    pTop->at(0)[3] = Top.Pz();

    pTop->at(1)[0] = Bjet.E();
    pTop->at(1)[1] = Bjet.Px();
    pTop->at(1)[2] = Bjet.Py();
    pTop->at(1)[3] = Bjet.Pz();

    pTop->at(2)[0] = Decay1.E();
    pTop->at(2)[1] = Decay1.Px();
    pTop->at(2)[2] = Decay1.Py();
    pTop->at(2)[3] = Decay1.Pz();

    pTop->at(3)[0] = Decay2.E();
    pTop->at(3)[1] = Decay2.Px();
    pTop->at(3)[2] = Decay2.Py();
    pTop->at(3)[3] = Decay2.Pz();
  }

  return;
}

void MEPhaseSpace::FillHiggsDecayPhaseSpacePoint(TLorentzVector& Higgs, TLorentzVector& P1, TLorentzVector& P2, TLorentzVector& P3, TLorentzVector& P4) const{

  pHiggs->at(0)[0] = Higgs.E();
  pHiggs->at(0)[1] = Higgs.Px();
  pHiggs->at(0)[2] = Higgs.Py();
  pHiggs->at(0)[3] = Higgs.Pz();

  pHiggs->at(1)[0] = P1.E();
  pHiggs->at(1)[1] = P1.Px();
  pHiggs->at(1)[2] = P1.Py();
  pHiggs->at(1)[3] = P1.Pz();

  pHiggs->at(2)[0] = P2.E();
  pHiggs->at(2)[1] = P2.Px();
  pHiggs->at(2)[2] = P2.Py();
  pHiggs->at(2)[3] = P2.Pz();

  pHiggs->at(3)[0] = P3.E();
  pHiggs->at(3)[1] = P3.Px();
  pHiggs->at(3)[2] = P3.Py();
  pHiggs->at(3)[3] = P3.Pz();

  pHiggs->at(4)[0] = P4.E();
  pHiggs->at(4)[1] = P4.Px();
  pHiggs->at(4)[2] = P4.Py();
  pHiggs->at(4)[3] = P4.Pz();

  if (iMode==kMEM_TTH_TopAntitopHiggsSemiLepDecay && MEMFix_HiggsSemiLep.LepSign == 1){
    pHiggs->at(1)[0] = P3.E();
    pHiggs->at(1)[1] = P3.Px();
    pHiggs->at(1)[2] = P3.Py();
    pHiggs->at(1)[3] = P3.Pz();

    pHiggs->at(2)[0] = P4.E();
    pHiggs->at(2)[1] = P4.Px();
    pHiggs->at(2)[2] = P4.Py();
    pHiggs->at(2)[3] = P4.Pz();

    pHiggs->at(3)[0] = P1.E();
    pHiggs->at(3)[1] = P1.Px();
    pHiggs->at(3)[2] = P1.Py();
    pHiggs->at(3)[3] = P1.Pz();

    pHiggs->at(4)[0] = P2.E();
    pHiggs->at(4)[1] = P2.Px();
    pHiggs->at(4)[2] = P2.Py();
    pHiggs->at(4)[3] = P2.Pz();
  }



  return;
}


void MEPhaseSpace::ApplyTotalTransverseBoost() const {

  if (verbosity>=2) cout << "Applying boost"<<endl;

  TLorentzVector Top(pCore->at(2)[1], pCore->at(2)[2], pCore->at(2)[3], pCore->at(2)[0]);
  TLorentzVector Antitop(pCore->at(3)[1], pCore->at(3)[2], pCore->at(3)[3], pCore->at(3)[0]);
 
  if (iCore==kTTH){
    TLorentzVector Higgs(pCore->at(4)[1], pCore->at(4)[2], pCore->at(4)[3], pCore->at(4)[0]);
    TLorentzVector Ptot = Top + Antitop + Higgs;
    TVector3 BoostPt = Ptot.BoostVector();
    BoostPt.SetZ(0);

    Top.Boost( -BoostPt);
    Antitop.Boost( -BoostPt);
    Higgs.Boost( -BoostPt);

    Ptot = Top + Antitop + Higgs;
    double x1 = (Ptot.Pz()+Ptot.E())/(comEnergy);
    double x2 = (-Ptot.Pz()+Ptot.E())/(comEnergy);
    SetInitialPartonMomenta(x1, x2);

    FillTTHPhaseSpacePoint(Top, Antitop, Higgs); 
    ReadPartonMomenta(pCore, 5);
  }
  if (iCore==kTTLL || iCore==kTTW){
    TLorentzVector Lep1(pCore->at(4)[1], pCore->at(4)[2], pCore->at(4)[3], pCore->at(4)[0]);
    TLorentzVector Lep2(pCore->at(5)[1], pCore->at(5)[2], pCore->at(5)[3], pCore->at(5)[0]);
    TLorentzVector Ptot = Top + Antitop + Lep1 + Lep2;
    TVector3 BoostPt = Ptot.BoostVector();
    BoostPt.SetZ(0);
    if (verbosity>=2) cout << "Boost Px="<<BoostPt.X()<<" Py="<<BoostPt.Py()<<" Pz="<<BoostPt.Z()<<endl;

    Top.Boost( -BoostPt);
    Antitop.Boost( -BoostPt);
    Lep1.Boost( -BoostPt);

    Lep2.Boost( -BoostPt);
    if (verbosity>=2) cout << "Boost Lep2 Px="<< Lep2.Px()<<" Py="<<Lep2.Py()<<" Pz="<<Lep2.Pz()<<" E="<< Lep2.E()<<endl;

    double Lep2_Px = -Top.Px()-Antitop.Px()-Lep1.Px();
    double Lep2_Py = -Top.Py()-Antitop.Py()-Lep1.Py();
    double Lep2_Pz = Lep2.Pz();
    Lep2.SetPxPyPzE(Lep2_Px, Lep2_Py, Lep2_Pz, sqrt(Lep2_Px*Lep2_Px+Lep2_Py*Lep2_Py+Lep2_Pz*Lep2_Pz));
    if (verbosity>=2) cout << "Conservation Lep2 Px="<<Lep2_Px<<" Py="<<Lep2_Py<<" Pz="<<Lep2_Pz<<" E="<<Lep2.E()<<endl;

    Ptot = Top + Antitop + Lep1 + Lep2;
    double x1 = (Ptot.Pz()+Ptot.E())/(comEnergy);
    double x2 = (-Ptot.Pz()+Ptot.E())/(comEnergy);
    SetInitialPartonMomenta(x1, x2);
    FillTTLLPhaseSpacePoint(Top, Antitop, Lep1, Lep2);
    ReadPartonMomenta(pCore, 6);
  }
  if (iCore==kTTWJJ){
    TLorentzVector Lep1(pCore->at(4)[1], pCore->at(4)[2], pCore->at(4)[3], pCore->at(4)[0]);
    TLorentzVector Lep2(pCore->at(5)[1], pCore->at(5)[2], pCore->at(5)[3], pCore->at(5)[0]);
    TLorentzVector Jet1(pCore->at(6)[1], pCore->at(6)[2], pCore->at(6)[3], pCore->at(6)[0]);
    TLorentzVector Jet2(pCore->at(7)[1], pCore->at(7)[2], pCore->at(7)[3], pCore->at(7)[0]);
 
   TLorentzVector Ptot = Top + Antitop + Lep1 + Lep2 + Jet1 + Jet2;
    TVector3 BoostPt = Ptot.BoostVector();
    BoostPt.SetZ(0);
    if (verbosity>=2) cout << "Boost Px="<<BoostPt.X()<<" Py="<<BoostPt.Py()<<" Pz="<<BoostPt.Z()<<endl;

    Top.Boost( -BoostPt);
    Antitop.Boost( -BoostPt);
    Lep1.Boost( -BoostPt);
    Lep2.Boost( -BoostPt);
    Jet1.Boost( -BoostPt);
    Jet2.Boost( -BoostPt);
    if (verbosity>=2) cout << "Boost Jet2 Px="<< Jet2.Px()<<" Py="<<Jet2.Py()<<" Pz="<<Jet2.Pz()<<" E="<< Jet2.E()<<endl;

    double Jet2_Px = -Top.Px()-Antitop.Px()-Lep1.Px()-Lep2.Px()-Jet1.Px();
    double Jet2_Py = -Top.Py()-Antitop.Py()-Lep1.Py()-Lep2.Py()-Jet1.Py();
    double Jet2_Pz = Jet2.Pz();
    Jet2.SetPxPyPzE(Jet2_Px, Jet2_Py, Jet2_Pz, sqrt(Jet2_Px*Jet2_Px+Jet2_Py*Jet2_Py+Jet2_Pz*Jet2_Pz));
    if (verbosity>=2) cout << "Conservation Jet2 Px="<<Jet2_Px<<" Py="<<Jet2_Py<<" Pz="<<Jet2_Pz<<" E="<<Jet2.E()<<endl;

    Ptot = Top + Antitop + Lep1 + Lep2 + Jet1 + Jet2;
    double x1 = (Ptot.Pz()+Ptot.E())/(comEnergy);
    double x2 = (-Ptot.Pz()+Ptot.E())/(comEnergy);
    SetInitialPartonMomenta(x1, x2);
    FillTTLNuJJPhaseSpacePoint(Top, Antitop, Lep1, Lep2, Jet1, Jet2);
    ReadPartonMomenta(pCore, 8);
  }

}

double MEPhaseSpace::Eval(const double* x) const {

  if (verbosity>=1 && iIteration % 1000==0) cout << "Iteration "<<iIteration<<", non-zero "<< iCall<<endl;
  iIteration++;

  //cout << "Bjorken x: x1="<<x[0]<<" x2="<<x[1]<<endl;
  if (verbosity>=2) 
  cout << "***************** Weight evaluation iCall="<<iCall<<" *****************"<<endl; 


  double weight = 0;

  if (iMode == kNone){
  //Evaluate weight only
    double weightME = ComputeMatrixElement();
    double weightPDF = ComputePDF(pCore->at(0)[0]*2/comEnergy, pCore->at(1)[0]*2/comEnergy, muF);
    weight = weightME * weightPDF;
  }

  if (iMode == kInitialPartons){
    //MEM mode: integrate over x1,x2
    SetInitialPartonMomenta(x[0], x[1]);
    double weightME = ComputeMatrixElement();
    double weightPDF = ComputePDF(x[0], x[1], muF);
    weight = weightME * weightPDF;
  }
  if (iMode==kMEM_TTH_TopAntitopHiggsDecay){

    muF = (2*mTop+mHiggs)/2.;
 
    xMEM[0] = x[0]; //TopHad, Bjet_E
    xMEM[1] = MEMFix_TopHad.Bjet_Theta;
    xMEM[2] = MEMFix_TopHad.Bjet_Phi; 
    xMEM[3] = x[1]; //TopHad, Jet1_E
    xMEM[4] = MEMFix_TopHad.Jet1_Theta; 
    xMEM[5] = MEMFix_TopHad.Jet1_Phi;
    xMEM[6] = MEMFix_TopHad.Jet2_Theta;
    xMEM[7] = MEMFix_TopHad.Jet2_Phi;

    xMEM[8] = x[2]; //TopLep, Bjet_E
    xMEM[9] = MEMFix_TopLep.Bjet_Theta; 
    xMEM[10] = MEMFix_TopLep.Bjet_Phi;  
    xMEM[11] = MEMFix_TopLep.Lep_E;
    xMEM[12] = MEMFix_TopLep.Lep_Theta;
    xMEM[13] = MEMFix_TopLep.Lep_Phi;
    xMEM[14] = x[3]; //TopLep, Neut_Theta
    xMEM[15] = x[4]; //TopLep, Neut_Phi

    xMEM[16] = MEMFix_HiggsFullLep.Lep1_E;
    xMEM[17] = MEMFix_HiggsFullLep.Lep1_Theta;
    xMEM[18] = MEMFix_HiggsFullLep.Lep1_Phi;
    xMEM[19] = MEMFix_HiggsFullLep.Lep2_E;
    xMEM[20] = MEMFix_HiggsFullLep.Lep2_Theta;
    xMEM[21] = MEMFix_HiggsFullLep.Lep2_Phi;
    xMEM[22] = x[5]; //HiggsFullLep, Neut1_E
    xMEM[23] = x[6]; //HiggsFullLep, Neut1_Theta
    xMEM[24] = x[7]; //HiggsFullLep, Neut1_Phi
    xMEM[25] = x[8]; //HiggsFullLep, Neut2_Theta
    xMEM[26] = x[9]; //HiggsFullLep, Neut2_Phi

    Computed_mETvect.SetPxPyPzE(0,0,0,0);

    double weightPS = SetupKinematicsTTH_NoBjorken_TopAntitopHiggsDecay(xMEM);
    if (verbosity>=2) cout << "weightPS="<<weightPS<<endl;
    if (weightPS==0) {
      errorCounter[kErr_PS_Product]++;
      return 0;
    }

 if (isTopAntitop==0){
    ComputedVarForTF.Bjet1_E = pTop->at(1)[0];
    ComputedVarForTF.Bjet2_E = pAntitop->at(1)[0];
  }
  else if (isTopAntitop==1){
    ComputedVarForTF.Bjet1_E = pAntitop->at(1)[0];
    ComputedVarForTF.Bjet2_E = pTop->at(1)[0];
  }
    if (isTopAntitop==0) {
      ComputedVarForTF.Jet1_E = pTop->at(2)[0];
      ComputedVarForTF.Jet2_E = pTop->at(3)[0];
    }
    else if (isTopAntitop==1){
      ComputedVarForTF.Jet1_E = pAntitop->at(2)[0];
      ComputedVarForTF.Jet2_E = pAntitop->at(3)[0];
    }

    ComputedVarForTF.mET_Px = Computed_mETvect.Px();
    ComputedVarForTF.mET_Py = Computed_mETvect.Py();
    double Gen_Bjet1_E = ComputedVarForTF.Bjet1_E;
    double Gen_Bjet2_E = ComputedVarForTF.Bjet2_E;
    double Gen_Jet1_E = ComputedVarForTF.Jet1_E;
    double Gen_Jet2_E = ComputedVarForTF.Jet2_E; 
    double Gen_Recoil_Px = -(pCore->at(2)[1]+pCore->at(3)[1]+pCore->at(4)[1]);
    double Gen_Recoil_Py = -(pCore->at(2)[2]+pCore->at(3)[2]+pCore->at(4)[2]);
    double weightTF = ComputeTFProduct(Gen_Bjet1_E, Gen_Bjet2_E, Gen_Jet1_E, Gen_Jet2_E, Gen_Recoil_Px, Gen_Recoil_Py);
    if (weightTF==0) { errorCounter[kErr_TF_Product]++; return 0; }

    ApplyTotalTransverseBoost();
    double weightME = ComputeMatrixElement() * GeV2barn;
    if (weightME==0) { errorCounter[kErr_ME]++; return 0; }

    muF = (2*mTop+mHiggs)/2.;
    double weightPDF = ComputePDF(pCore->at(0)[0]*2/comEnergy, pCore->at(1)[0]*2/comEnergy, muF);
    if (weightPDF==0) { errorCounter[kErr_PDF]++; return 0; }

    weight = weightME * weightPDF * weightPS * weightTF;
    if (weight==0) { errorCounter[kErr_Weight_Product]++; return 0;}
  }
  if (iMode==kMEM_TTH_TopAntitopHiggsSemiLepDecay){

    xMEM[0] = x[0]; //TopLep, Bjet_E
    xMEM[1] = MEMFix_TopLep.Bjet_Theta;
    xMEM[2] = MEMFix_TopLep.Bjet_Phi;
    xMEM[3] = MEMFix_TopLep.Lep_E;
    xMEM[4] = MEMFix_TopLep.Lep_Theta;
    xMEM[5] = MEMFix_TopLep.Lep_Phi;
    xMEM[6] = x[1]; //TopLep, Neut_Theta
    xMEM[7] = x[2]; //TopLep, Neut_Phi

    xMEM[8] = x[3]; //TopLep2, Bjet_E
    xMEM[9] = MEMFix_TopLep2.Bjet_Theta;
    xMEM[10] = MEMFix_TopLep2.Bjet_Phi;
    xMEM[11] = MEMFix_TopLep2.Lep_E;
    xMEM[12] = MEMFix_TopLep2.Lep_Theta;
    xMEM[13] = MEMFix_TopLep2.Lep_Phi;
    xMEM[14] = x[4]; //TopLep2, Neut_Theta
    xMEM[15] = x[5]; //TopLep2, Neut_Phi

    xMEM[16] = x[8]; //HiggsSemiLep, Jet1_E
    xMEM[17] = MEMFix_HiggsSemiLep.Jet1_Theta;
    xMEM[18] = MEMFix_HiggsSemiLep.Jet1_Phi; 
    xMEM[19] = x[9]; //HiggsSemiLep, Jet2_E
    xMEM[20] = MEMFix_HiggsSemiLep.Jet2_Theta;
    xMEM[21] = MEMFix_HiggsSemiLep.Jet2_Phi;
    xMEM[22] = MEMFix_HiggsSemiLep.Lep1_E;
    xMEM[23] = MEMFix_HiggsSemiLep.Lep1_Theta;
    xMEM[24] = MEMFix_HiggsSemiLep.Lep1_Phi;
    xMEM[25] = x[6]; //HiggsSemiLep, Neut_Theta
    xMEM[26] = x[7]; //HiggsSemiLep, Neut_Phi

    Computed_mETvect.SetPxPyPzE(0,0,0,0);

    double weightPS = SetupKinematicsTTH_NoBjorken_TopAntitopHiggsDecay(xMEM);
    if (verbosity>=2) cout << "weightPS="<<weightPS<<endl;
    if (weightPS==0) { errorCounter[kErr_PS_Product]++; return 0; }

  if (isTopAntitop==0){
    ComputedVarForTF.Bjet1_E = pTop->at(1)[0];
    ComputedVarForTF.Bjet2_E = pAntitop->at(1)[0];
  }
  else if (isTopAntitop==1){
    ComputedVarForTF.Bjet1_E = pAntitop->at(1)[0];
    ComputedVarForTF.Bjet2_E = pTop->at(1)[0];
  }
    if (MEMFix_HiggsSemiLep.LepSign == 1){
      ComputedVarForTF.Jet1_E = pHiggs->at(3)[0];
      ComputedVarForTF.Jet2_E = pHiggs->at(4)[0];
    }
    else if (MEMFix_HiggsSemiLep.LepSign == -1){
      ComputedVarForTF.Jet1_E = pHiggs->at(1)[0];
      ComputedVarForTF.Jet2_E = pHiggs->at(2)[0];
    }

    ComputedVarForTF.mET_Px = Computed_mETvect.Px();
    ComputedVarForTF.mET_Py = Computed_mETvect.Py();
    double Gen_Bjet1_E = ComputedVarForTF.Bjet1_E;
    double Gen_Bjet2_E = ComputedVarForTF.Bjet2_E;
    double Gen_Jet1_E = ComputedVarForTF.Jet1_E; 
    double Gen_Jet2_E = ComputedVarForTF.Jet2_E; 
    double Gen_Recoil_Px = -(pCore->at(2)[1]+pCore->at(3)[1]+pCore->at(4)[1]);
    double Gen_Recoil_Py = -(pCore->at(2)[2]+pCore->at(3)[2]+pCore->at(4)[2]);
    double weightTF = ComputeTFProduct(Gen_Bjet1_E, Gen_Bjet2_E, Gen_Jet1_E, Gen_Jet2_E, Gen_Recoil_Px, Gen_Recoil_Py);
    if (weightTF==0) { errorCounter[kErr_TF_Product]++; return 0; }

    ApplyTotalTransverseBoost();
    double weightME = ComputeMatrixElement() * GeV2barn;
    if (weightME==0) { errorCounter[kErr_ME]++; return 0; }

    muF = (2*mTop+mHiggs)/2.;
    double weightPDF = ComputePDF(pCore->at(0)[0]*2/comEnergy, pCore->at(1)[0]*2/comEnergy, muF);
    if (weightPDF==0) { errorCounter[kErr_PDF]++; return 0; }

    weight = weightME * weightPDF * weightPS * weightTF;
    if (weight==0) { errorCounter[kErr_Weight_Product]++; return 0;}
  }
  if (iMode==kMEM_TTW_TopAntitopDecay){

    xMEM[0] = x[0]; //TopLep, Bjet_E
    xMEM[1] = MEMFix_TopLep.Bjet_Theta;
    xMEM[2] = MEMFix_TopLep.Bjet_Phi;
    xMEM[3] = MEMFix_TopLep.Lep_E;
    xMEM[4] = MEMFix_TopLep.Lep_Theta;
    xMEM[5] = MEMFix_TopLep.Lep_Phi;
    xMEM[6] = x[1]; //TopLep, Neut_Theta
    xMEM[7] = x[2]; //TopLep, Neut_Phi

    xMEM[8] = x[3]; //TopLep2, Bjet_E
    xMEM[9] = MEMFix_TopLep2.Bjet_Theta;
    xMEM[10] = MEMFix_TopLep2.Bjet_Phi;
    xMEM[11] = MEMFix_TopLep2.Lep_E;
    xMEM[12] = MEMFix_TopLep2.Lep_Theta;
    xMEM[13] = MEMFix_TopLep2.Lep_Phi;
    xMEM[14] = x[4]; //TopLep2, Neut_Theta
    xMEM[15] = x[5]; //TopLep2, Neut_Phi

    xMEM[16] = MEMFix_HiggsSemiLep.Lep1_E;
    xMEM[17] = MEMFix_HiggsSemiLep.Lep1_Theta;
    xMEM[18] = MEMFix_HiggsSemiLep.Lep1_Phi;
    xMEM[19] = x[6]; //Neut_E (//tW, ie mW transform later)
    xMEM[20] = x[7]; //Neut_Theta
    xMEM[21] = x[8]; //Neut_Phi

    Computed_mETvect.SetPxPyPzE(0,0,0,0);

    double weightPS = SetupKinematicsTTW_NoBjorken_TopAntitopDecay(xMEM);
    if (verbosity>=2) cout << "weightPS="<<weightPS<<endl;
    if (weightPS==0) { errorCounter[kErr_PS_Product]++; return 0; }

  if (isTopAntitop==0){
    ComputedVarForTF.Bjet1_E = pTop->at(1)[0];
    ComputedVarForTF.Bjet2_E = pAntitop->at(1)[0];
  }
  else if (isTopAntitop==1){
    ComputedVarForTF.Bjet1_E = pAntitop->at(1)[0];
    ComputedVarForTF.Bjet2_E = pTop->at(1)[0];
  }
  if (isTopAntitop==0) {
      ComputedVarForTF.Jet1_E = pTop->at(2)[0];
      ComputedVarForTF.Jet2_E = pTop->at(3)[0];
  }
  else if (isTopAntitop==1){
      ComputedVarForTF.Jet1_E = pAntitop->at(2)[0];
      ComputedVarForTF.Jet2_E = pAntitop->at(3)[0];
  }

    ComputedVarForTF.mET_Px = Computed_mETvect.Px();
    ComputedVarForTF.mET_Py = Computed_mETvect.Py();
    double Gen_Bjet1_E = ComputedVarForTF.Bjet1_E;
    double Gen_Bjet2_E = ComputedVarForTF.Bjet2_E;
    double Gen_Jet1_E = ComputedVarForTF.Jet1_E;
    double Gen_Jet2_E = ComputedVarForTF.Jet2_E;
    double Gen_Recoil_Px = -(pCore->at(2)[1]+pCore->at(3)[1]+pCore->at(4)[1]+pCore->at(5)[1]);
    double Gen_Recoil_Py = -(pCore->at(2)[2]+pCore->at(3)[2]+pCore->at(4)[2]+pCore->at(5)[2]);
    double weightTF = ComputeTFProduct(Gen_Bjet1_E, Gen_Bjet2_E, Gen_Jet1_E, Gen_Jet2_E, Gen_Recoil_Px, Gen_Recoil_Py);
    if (weightTF==0) { errorCounter[kErr_TF_Product]++; return 0; }

    ApplyTotalTransverseBoost();
    muF = (2*mTop+mW)/2.;
    double weightMEPDF = ConvolvePdfCrossSection(pCore->at(0)[0]*2/comEnergy, pCore->at(1)[0]*2/comEnergy, muF) * GeV2barn;
    if (weightMEPDF==0) { errorCounter[kErr_ME]++; errorCounter[kErr_PDF]++; return 0; }

    weight = weightMEPDF * weightPS * weightTF;
    if (weight==0) { errorCounter[kErr_Weight_Product]++; return 0;}
  }
  if (iMode==kMEM_TTWJJ_TopAntitopDecay){

    xMEM[0] = x[0]; //TopLep, Bjet_E
    xMEM[1] = MEMFix_TopLep.Bjet_Theta;
    xMEM[2] = MEMFix_TopLep.Bjet_Phi;
    xMEM[3] = MEMFix_TopLep.Lep_E;
    xMEM[4] = MEMFix_TopLep.Lep_Theta;
    xMEM[5] = MEMFix_TopLep.Lep_Phi;
    xMEM[6] = x[1]; //TopLep, Neut_Theta
    xMEM[7] = x[2]; //TopLep, Neut_Phi

    xMEM[8] = x[3]; //TopLep2, Bjet_E
    xMEM[9] = MEMFix_TopLep2.Bjet_Theta;
    xMEM[10] = MEMFix_TopLep2.Bjet_Phi;
    xMEM[11] = MEMFix_TopLep2.Lep_E;
    xMEM[12] = MEMFix_TopLep2.Lep_Theta;
    xMEM[13] = MEMFix_TopLep2.Lep_Phi;
    xMEM[14] = x[4]; //TopLep2, Neut_Theta
    xMEM[15] = x[5]; //TopLep2, Neut_Phi

    xMEM[16] = MEMFix_HiggsSemiLep.Lep1_E;
    xMEM[17] = MEMFix_HiggsSemiLep.Lep1_Theta;
    xMEM[18] = MEMFix_HiggsSemiLep.Lep1_Phi;
    xMEM[19] = x[6]; //Neut_E (//tW, ie mW transform later)
    xMEM[20] = x[7]; //Neut_Theta
    xMEM[21] = x[8]; //Neut_Phi

    xMEM[22] = x[9]; //Jet1_E
    xMEM[23] = MEMFix_HiggsSemiLep.Jet1_Theta; //Jet1_Theta
    xMEM[24] = MEMFix_HiggsSemiLep.Jet1_Phi; //Jet1_Phi    
    xMEM[25] = x[10]; //Jet2_E
    xMEM[26] = MEMFix_HiggsSemiLep.Jet2_Theta; //Jet2_Theta
    xMEM[27] = MEMFix_HiggsSemiLep.Jet2_Phi; //Jet2_Phi

    Computed_mETvect.SetPxPyPzE(0,0,0,0);

    double weightPS = SetupKinematicsTTWJJ_NoBjorken_TopAntitopDecay(xMEM);
    if (verbosity>=2) cout << "weightPS="<<weightPS<<endl;
    if (weightPS==0) { errorCounter[kErr_PS_Product]++; return 0; }

  if (isTopAntitop==0){
    ComputedVarForTF.Bjet1_E = pTop->at(1)[0];
    ComputedVarForTF.Bjet2_E = pAntitop->at(1)[0];
  }
  else if (isTopAntitop==1){
    ComputedVarForTF.Bjet1_E = pAntitop->at(1)[0];
    ComputedVarForTF.Bjet2_E = pTop->at(1)[0];
  }
  ComputedVarForTF.Jet1_E = pCore->at(6)[0];
  ComputedVarForTF.Jet2_E = pCore->at(7)[0];

  if (verbosity>=2) {
        cout << "Bjet1_E computed "<<ComputedVarForTF.Bjet1_E<<" input "<<x[0]<<endl;
        cout << "Bjet2_E computed "<<ComputedVarForTF.Bjet2_E<<" input "<<x[3]<<endl;
	cout << "Jet1_E computed "<<ComputedVarForTF.Jet1_E<<" input "<<x[9]<<endl;
	cout << "Jet2_E computed "<<ComputedVarForTF.Jet2_E<<" input "<<x[10]<<endl;
  }

    ComputedVarForTF.mET_Px = Computed_mETvect.Px();
    ComputedVarForTF.mET_Py = Computed_mETvect.Py();
    double Gen_Bjet1_E = ComputedVarForTF.Bjet1_E;
    double Gen_Bjet2_E = ComputedVarForTF.Bjet2_E;
    double Gen_Jet1_E = ComputedVarForTF.Jet1_E;
    double Gen_Jet2_E = ComputedVarForTF.Jet2_E;
    double Gen_Recoil_Px = -(pCore->at(2)[1]+pCore->at(3)[1]+pCore->at(4)[1]+pCore->at(5)[1]+pCore->at(6)[1]+pCore->at(7)[1]);
    double Gen_Recoil_Py = -(pCore->at(2)[2]+pCore->at(3)[2]+pCore->at(4)[2]+pCore->at(5)[2]+pCore->at(6)[2]+pCore->at(7)[2]);
    double weightTF = ComputeTFProduct(Gen_Bjet1_E, Gen_Bjet2_E, Gen_Jet1_E, Gen_Jet2_E, Gen_Recoil_Px, Gen_Recoil_Py);
    if (weightTF==0) { errorCounter[kErr_TF_Product]++; return 0; }

    ApplyTotalTransverseBoost();
    muF = (2*mTop+mW)/2.;
    double weightMEPDF = ConvolvePdfCrossSection(pCore->at(0)[0]*2/comEnergy, pCore->at(1)[0]*2/comEnergy, muF) * GeV2barn;
    if (weightMEPDF==0) { errorCounter[kErr_ME]++; errorCounter[kErr_PDF]++; return 0; }

    weight = weightMEPDF * weightPS * weightTF;
    if (weight==0) { errorCounter[kErr_Weight_Product]++; return 0;}
  }
  if (iMode==kMEM_TTLL_TopAntitopDecay){

    xMEM[0] = x[0]; //TopHad, Bjet_E
    xMEM[1] = MEMFix_TopHad.Bjet_Theta;
    xMEM[2] = MEMFix_TopHad.Bjet_Phi;
    xMEM[3] = x[1]; //TopHad, Jet1_E
    xMEM[4] = MEMFix_TopHad.Jet1_Theta;
    xMEM[5] = MEMFix_TopHad.Jet1_Phi;
    xMEM[6] = MEMFix_TopHad.Jet2_Theta;
    xMEM[7] = MEMFix_TopHad.Jet2_Phi;

    xMEM[8] = x[2]; //TopLep, Bjet_E
    xMEM[9] = MEMFix_TopLep.Bjet_Theta;
    xMEM[10] = MEMFix_TopLep.Bjet_Phi;
    xMEM[11] = MEMFix_TopLep.Lep_E;
    xMEM[12] = MEMFix_TopLep.Lep_Theta;
    xMEM[13] = MEMFix_TopLep.Lep_Phi;
    xMEM[14] = x[3]; //TopLep, Neut_Theta
    xMEM[15] = x[4]; //TopLep, Neut_Phi

    xMEM[16] = MEMFix_HiggsFullLep.Lep1_E;
    xMEM[17] = MEMFix_HiggsFullLep.Lep1_Theta;
    xMEM[18] = MEMFix_HiggsFullLep.Lep1_Phi;
    xMEM[19] = MEMFix_HiggsFullLep.Lep2_E;
    xMEM[20] = MEMFix_HiggsFullLep.Lep2_Theta;
    xMEM[21] = MEMFix_HiggsFullLep.Lep2_Phi;

    Computed_mETvect.SetPxPyPzE(0,0,0,0);

    double weightPS = SetupKinematicsTTLL_NoBjorken_TopAntitopDecay(xMEM);
    if (verbosity>=2) cout << "weightPS="<<weightPS<<endl;
    if (weightPS==0) { errorCounter[kErr_PS_Product]++; return 0; }

  if (isTopAntitop==0){
    ComputedVarForTF.Bjet1_E = pTop->at(1)[0];
    ComputedVarForTF.Bjet2_E = pAntitop->at(1)[0];
  }
  else if (isTopAntitop==1){
    ComputedVarForTF.Bjet1_E = pAntitop->at(1)[0];
    ComputedVarForTF.Bjet2_E = pTop->at(1)[0];
  }
  if (isTopAntitop==0) {
      ComputedVarForTF.Jet1_E = pTop->at(2)[0];
      ComputedVarForTF.Jet2_E = pTop->at(3)[0];
  }
  else if (isTopAntitop==1){
      ComputedVarForTF.Jet1_E = pAntitop->at(2)[0];
      ComputedVarForTF.Jet2_E = pAntitop->at(3)[0];
  }

    ComputedVarForTF.mET_Px = Computed_mETvect.Px();
    ComputedVarForTF.mET_Py = Computed_mETvect.Py();
    double Gen_Bjet1_E = ComputedVarForTF.Bjet1_E;
    double Gen_Bjet2_E = ComputedVarForTF.Bjet2_E;
    double Gen_Jet1_E = ComputedVarForTF.Jet1_E; 
    double Gen_Jet2_E = ComputedVarForTF.Jet2_E;
    double Gen_Recoil_Px = -(pCore->at(2)[1]+pCore->at(3)[1]+pCore->at(4)[1]+pCore->at(5)[1]);
    double Gen_Recoil_Py = -(pCore->at(2)[2]+pCore->at(3)[2]+pCore->at(4)[2]+pCore->at(5)[2]);
    double weightTF = ComputeTFProduct(Gen_Bjet1_E, Gen_Bjet2_E, Gen_Jet1_E, Gen_Jet2_E, Gen_Recoil_Px, Gen_Recoil_Py);
    if (weightTF==0) { errorCounter[kErr_TF_Product]++; return 0; }

    ApplyTotalTransverseBoost();
    double weightME = ComputeMatrixElement() * GeV2barn;
    if (weightME==0) { errorCounter[kErr_ME]++; return 0; }

    TLorentzVector Pl1, Pl2;
    SetMomentumFromEThetaPhi(0, MEMFix_HiggsFullLep.Lep1_E, MEMFix_HiggsFullLep.Lep1_Theta, MEMFix_HiggsFullLep.Lep1_Phi, &Pl1);
    SetMomentumFromEThetaPhi(0, MEMFix_HiggsFullLep.Lep2_E, MEMFix_HiggsFullLep.Lep2_Theta, MEMFix_HiggsFullLep.Lep2_Phi, &Pl2);
    double mll = (Pl1+Pl2).M();
    muF = (2*mTop+mll)/2.;
    double weightPDF = ComputePDF(pCore->at(0)[0]*2/comEnergy, pCore->at(1)[0]*2/comEnergy, muF);
    if (weightPDF==0) { errorCounter[kErr_PDF]++; return 0; }

    weight = weightME * weightPDF * weightPS * weightTF;
    if (weight==0) { errorCounter[kErr_Weight_Product]++; return 0;}
  }
  if (iMode==kAllPartonsTTH_TopDecay){
    double weightPS = SetupKinematicsTTH_NoBjorken_TopHadDecay(x); 
    if (weightPS==0) return 0;
    double weightME = ComputeMatrixElement() * GeV2barn;
    double weightPDF = ComputePDF(pCore->at(0)[0]*2/comEnergy, pCore->at(1)[0]*2/comEnergy, muF);
    weight = weightME * weightPDF * weightPS;
  }
  if (iMode == kAllPartonsTTLL){
     double weightPS = SetupKinematics2to4_LabFrame_NoBjorken(x);
      if (weightPS==0) return 0;
      double weightME = ComputeMatrixElement() * GeV2barn;
      double weightPDF = ComputePDF(pCore->at(0)[0]*2/comEnergy, pCore->at(1)[0]*2/comEnergy, muF);
      weight = weightME * weightPDF * weightPS;
  }
  if (iMode == kAllPartonsTTH){
    if (iOption==kTwoBjorken){
      SetInitialPartonMomenta(x[0], x[1]);
      double weightPS = SetupKinematicsTTH(x);
      if (weightPS==0) return 0;
      double weightME = ComputeMatrixElement() * GeV2barn;
      double weightPDF = ComputePDF(x[0], x[1], muF);
      weight = weightME * weightPDF * weightPS;
    }
    if (iOption==kOneBjorken){
      double weightPS = SetupKinematics2to3_LabFrame_OneBjorken(x);
      if (weightPS==0) return 0;
      double weightME = ComputeMatrixElement() * GeV2barn;
      double weightPDF = ComputePDF(pCore->at(0)[0]*2/comEnergy, pCore->at(1)[0]*2/comEnergy, muF);
      weight = weightME * weightPDF * weightPS;
    }
    if (iOption==kTwoBjorken1to3){
      SetInitialPartonMomenta(x[5], x[6]);
      double weightPS = SetupKinematics1to3_LabFrame(x); 
      if (weightPS==0) return 0;
      double weightME = ComputeMatrixElement() * GeV2barn;
      //double weightPDF = ComputePDF(x[5], x[6], muF);
      double weightPDF = ComputePDF(pCore->at(0)[0]*2/comEnergy, pCore->at(1)[0]*2/comEnergy, muF);
      weight = weightME * weightPDF * weightPS;
    }
    if (iOption==kNoBjorken || iOption==kNoBjorkenMomentum){
      double weightPS = SetupKinematics2to3_LabFrame_NoBjorken(x);
      if (weightPS==0) return 0;
      double weightME = ComputeMatrixElement() * GeV2barn;
      double weightPDF = ComputePDF(pCore->at(0)[0]*2/comEnergy, pCore->at(1)[0]*2/comEnergy, muF);
      weight = weightME * weightPDF * weightPS;
    }
    /*
    if (iOption==kTwoBjorkenSubMasses){
      SetInitialPartonMomenta(x[5], x[6]);
      double weightPS = SetupKinematicsTopLep_FixedTopM_FiniteWwidth(x);
      if (weightPS==0) return 0;
      double weightME = ComputeMatrixElement() * GeV2barn;
      double weightPDF = ComputePDF(p->at(0)[0]*2/comEnergy, p->at(1)[0]*2/comEnergy, muF);
      weight = weightME * weightPDF * weightPS;
    }
    */
    //else weight = 0;
  }
  /*
  if (iMode == kAllPartonsTTH_FSonly){
     y[0] = p->at(0)[0]*2/comEnergy;
     y[1] = p->at(1)[0]*2/comEnergy;
     y[2] = x[0];
     y[3] = x[1];
     y[4] = x[2];
     y[5] = x[3];
     y[6] = x[4];
     double weightPS = SetupKinematicsTTH(y);
     if (weightPS==0) return 0;
     double weightME = ComputeMatrixElement() * GeV2barn;
     weight = weightME * weightPS;
  }
  */
  if (iMode==kAllPartonsTopHad){
    double weightPS = SetupKinematicsTopHad(x);
    if (weightPS==0) return 0;
    double weightME = ComputeMatrixElement();
    weight = weightME * weightPS;
  }
  if (iMode==kAllPartonsTopHad_FiniteWidth){
    double weightPS = SetupKinematicsTopHad_FiniteWidths(x);
    if (weightPS==0) return 0;
    double weightME = ComputeMatrixElement();
    weight = weightME * weightPS;
  }
  if (iMode==kAllPartonsTopLep){
    double weightPS = SetupKinematicsTopLep(x);
    if (weightPS==0) return 0;
    double weightME = ComputeMatrixElement();
    weight = weightME * weightPS;
  }
  if (iMode==kAllPartonsTopLep_FixedTopM_FiniteWwidth){
    double weightPS = 0;
    if (iOption==kFixMw) {
      y[0] = x[0];
      y[1] = x[1];
      y[2] = x[2];
      y[3] = x[3];
      y[4] = mW;
      weightPS = SetupKinematicsTopLep_FixedTopM_FiniteWwidth(y);
    }
    else if (iOption==kFixBenergy){
      y[0] = Benergy;
      y[1] = x[0];
      y[2] = x[1];
      y[3] = x[2];
      y[4] = x[3];
      weightPS = SetupKinematics1to3_LabFrame(y);
    }
    else weightPS = SetupKinematicsTopLep_FixedTopM_FiniteWwidth(x);
    //else weightPS = SetupKinematics1to3_LabFrame(x);
    if (weightPS==0) return 0;
    double weightME = ComputeMatrixElement() ;
    weight = weightME * weightPS;
  }
  if (iMode==kAllPartonsTopHad_FixedTopM_FiniteWwidth || iMode==kAllPartonsAntiTopHad_FixedTopM_FiniteWwidth){
   double weightPS = 0;
    if (iOption==kFixMw) {
      y[0] = x[0];
      y[1] = x[1];
      y[2] = x[2];
      y[3] = x[3];
      y[4] = mW;
      weightPS = SetupKinematicsTopLep_FixedTopM_FiniteWwidth(y);
    }
    else if (iOption==kFixBenergy){
      y[0] = Benergy;
      y[1] = x[0];
      y[2] = x[1];
      y[3] = x[2];
      y[4] = x[3];
      weightPS = SetupKinematics1to3_LabFrame(y);
    }
    //else  weightPS = SetupKinematicsTopHad_FixedTopM_FiniteWwidth(x);
    else weightPS = SetupKinematics1to3_LabFrame(x);
    if (weightPS==0) return 0;
    double weightME = ComputeMatrixElement();
    weight = weightME * weightPS;
  }
  if (iMode==kAllPartonsHiggsWWLep_FixedHiggsM_FiniteWwidth){
    double weightPS = SetupKinematicsHiggsWWLep_FixedHiggsM_FiniteWwidth(x);
    if (weightPS==0) return 0;
    double weightME = ComputeMatrixElement();
    weight = weightME * weightPS;
  }


  if (verbosity>=2) cout << "Evaluation call "<<iCall<<" weight="<< weight<<endl;
  iCall++;

  return weight; 

}

double MEPhaseSpace::SetupKinematicsTTH(const double* x) const{

  //arxiv:0210426 parametrization (Feynarts)
  //x[0]: Bjorken x1
  //x[1]: Bjorken x2
  //x[2]: Top energy E2 
  //x[3]: Top azimuthal angle phi2
  //x[4]: Higgs energy E4
  //x[5]: Higgs cos(theta) angle theta4
  //x[6]: Higgs phi angle phi4

  if (verbosity>=2) {
  cout << "Input parameters:"<<endl;
  cout << "x[0]="<<x[0]<<" bjorken x1"<<endl;
  cout << "x[1]="<<x[1]<<" bjorken x2"<<endl;
  cout << "x[2]="<<x[2]<<" top energy E2"<<endl;
  cout << "x[3]="<<x[3]<<" top azimuthal angle phi2"<<endl;
  cout << "x[4]="<<x[4]<<" Higgs energy E4"<<endl;
  cout << "x[5]="<<x[5]<<" Higgs angle theta4"<<endl;
  cout << "x[6]="<<x[6]<<" Higgs phi angle phi4"<<endl;
  }

  //Check Energy conservation
  if (x[4]+x[2]>(x[0]+x[1])/2.*comEnergy) {
    if (verbosity>=2) cout << "Top + Higgs energy > Tot energy, weight=0"<<endl;
    return 0;
  }
  if (x[4]>(x[0]+x[1])/2.*comEnergy-x[2]-mTop){
    if (verbosity>=2) cout << "Higgs energy too high, weight=0"<<endl;
    return 0;
  }
  if (x[2]>(x[0]+x[1])/2.*comEnergy-x[4]-mTop){
    if (verbosity>=2) cout << "Top energy too high, weight=0"<<endl;
    return 0;
  }

  //Reconstruct total momentum
  double x1 = x[0];
  double x2 = x[1];
  TLorentzVector Ptot(0, 0, (x1-x2)*comEnergy/2., (x1+x2)*comEnergy/2.);

  if (verbosity>=2)   cout << "Ptot energy="<< Ptot.E()<<" norm="<<TMath::Abs((x1-x2)*comEnergy/2.)<<endl;

  //Reconstruct Higgs
  double Higgs_M = mHiggs;
  double Higgs_E = x[4];
  //double Higgs_theta = TMath::ACos(x[5]);
  double Higgs_theta = x[5];
  double Higgs_phi = x[6];
  double Higgs_Eta = -TMath::Log(tan(Higgs_theta/2.));
  //double Higgs_Eta = -TMath::Log(sqrt(1-x[5]*x[5])/x[5]);
  double Higgs_Pt = sqrt(Higgs_E*Higgs_E-Higgs_M*Higgs_M)/cosh(Higgs_Eta);
  TLorentzVector Higgs;
  Higgs.SetPtEtaPhiE(Higgs_Pt, Higgs_Eta, Higgs_phi, Higgs_E);
  if (verbosity>=2) cout << "Higgs E="<<Higgs_E<<" norm="<<Higgs.Vect().Mag()<<endl;

  //Reconstruct Top
  TLorentzVector Top;
  double Top_M = mTop;
  double Top_E = x[2];
  double Top_phi = x[3];
  double Top_p3norm2 = Top_E*Top_E-Top_M*Top_M;
  if (verbosity>=2) cout << "Top E="<<Top_E<<" norm="<<sqrt(Top_p3norm2)<<endl;

 //Reconstruct Antitop
  double Antitop_M = mTop;
  double Antitop_E = Ptot.E()-Higgs.E()-Top_E;
  double Antitop_p3norm2 = Antitop_E*Antitop_E-Antitop_M*Antitop_M;
  if (verbosity>=2) cout << "Antitop E="<<Antitop_E<<" norm="<<sqrt(Antitop_p3norm2)<<endl;

  if (Top_p3norm2<0 || Antitop_E<0){
    if (verbosity>=2) cout << "Top or antitop has too low energy, mass constraint not respected"<<endl;
    return 0;
  } 

  if (Ptot.Vect().Mag()>Higgs.Vect().Mag()+sqrt(Top_p3norm2)+sqrt(Antitop_p3norm2)){
    if (verbosity>=2) cout << "Ptot norm=" << Ptot.Vect().Mag()<<" > Higgs + Top + Antitop in norm"<<endl;
    return 0;
  }
  if (Higgs.Vect().Mag()>Ptot.Vect().Mag()+sqrt(Top_p3norm2)+sqrt(Antitop_p3norm2)){
    if (verbosity>=2) cout << "Higgs norm=" << Higgs.Vect().Mag()<<" > Tot + Top + Antitop"<<endl;
    return 0;
  }
  if (sqrt(Top_p3norm2)>Higgs.Vect().Mag()+Ptot.Vect().Mag()+sqrt(Antitop_p3norm2)){
    if (verbosity>=2) cout << "Top norm="<<sqrt(Top_p3norm2)<<" > Higgs + Tot + Antitop"<<endl;
    return 0;
  }
  if (sqrt(Antitop_p3norm2)>Higgs.Vect().Mag()+Ptot.Vect().Mag()+sqrt(Top_p3norm2)){
    if (verbosity>=2) cout << "Antitop norm="<<sqrt(Antitop_p3norm2)<<" > Higgs + Top + Top"<<endl;
    return 0;
  }
  if ((Ptot-Higgs).Vect().Mag()>sqrt(Top_p3norm2)+sqrt(Antitop_p3norm2)){
    if (verbosity>=2) cout << "Ptot-Higgs norm=" << (Ptot-Higgs).Vect().Mag()<<" > Top + Antitop in norm"<<endl;
    return 0;
  }
  //if ((Ptot+Higgs).Vect().Mag()>4*Higgs.Vect().Mag()+sqrt(Top_p3norm2)+sqrt(Antitop_p3norm2)){
  //  cout << "Ptot+Higgs norm=" << (Ptot+Higgs).Vect().Mag()<<" > Top + Antitop in norm"<<endl;
  //  return 0;
  //}


  double Ptot_norm2 = Ptot.Vect().Mag2();
  double Higgs_norm2 = Higgs.Vect().Mag2();
  //double Higgs_norm_bis = sqrt(Higgs.Px()*Higgs.Px()+Higgs.Py()*Higgs.Py()+Higgs.Pz()*Higgs.Pz());
    if (verbosity>=2) cout << "Higgs_norm="<<sqrt(Higgs_norm2)<<" Ptot_norm="<<sqrt(Ptot_norm2)<<" Top_norm="<<sqrt(Top_p3norm2)<<" Antitop_norm="<<sqrt(Antitop_p3norm2)<<endl;
 
  //Solving 2nd order equation for Top Pt 
  double A = Ptot_norm2-Higgs_norm2+Antitop_p3norm2-Top_p3norm2-2*Ptot.Pz()*(Ptot.Pz()-Higgs.Pz());
  //double A = -(Ptot.Mag2()-Higgs_M*Higgs_M-2*Ptot.E()*Antitop_E-2*Higgs.E()*Top_E+2*Ptot.Pz()*(Ptot.Pz()-Higgs.Pz()));
  //double A = -Top_p3norm2+Antitop_p3norm2-(Ptot.Pz()-Higgs.Pz())*(Ptot.Pz()-Higgs.Pz())-Higgs.Pt()*Higgs.Pt();
  double B = 2*(cos(Top_phi)*Higgs.Px()+sin(Top_phi)*Higgs.Py());
  double C = 2*(Higgs.Pz()-Ptot.Pz());
  //cout << "A="<<A<<" B="<<B<<" C="<<C<<endl;
  double a = 1;
  double b = (-2*A*B)/(B*B+C*C);
  double c = (A*A-C*C*Top_p3norm2)/(B*B+C*C);
  //cout <<"a="<<a<<" b="<<b<<" c="<<c<<endl;
  double delta = b*b-4*a*c;

  if (delta<0) {
    if (verbosity>=2) cout << "Equation 2nd order for Top pT has delta<0"<<endl;
    return 0;
  }

  double X1 = (-b-sqrt(delta))/(2.*a);
  double X2 = (-b+sqrt(delta))/(2.*a);
  if (verbosity>=2) cout << "delta="<<delta<<" X1="<<X1<<" X2="<<X2<<endl;
  double Top_Pt=0;
  double Top_Eta=0;
  TLorentzVector Antitop;

  if (X1<0 && X2<0){
    if (verbosity>=2) cout << "Equation 2nd order for Top pT has two negative solutions"<<endl;
    return 0;
  }
  else if (X1>0 && X2<0) Top_Pt=X1;
  else if (X1<0 && X2>0) Top_Pt=X2;
  else if (X1>0 && X2>0) {
    if (verbosity>=2) cout << "Equation 2nd order for Top pT has two positive solutions, trying both"<<endl;
    Top_Pt = X1;
    double Top_Eta_tmp = TMath::ACosH(sqrt(Top_p3norm2)/Top_Pt);
    double Top_Pz_tmp = Top_Pt*TMath::SinH(Top_Eta_tmp);
    double Antitop_Pt_tmp = sqrt(Higgs.Pt()*Higgs.Pt()+B*Top_Pt+Top_Pt*Top_Pt);
    double Antitop_Pz_tmp = Ptot.Pz()-Higgs.Pz()-Top_Pz_tmp;
    double Antitop_Eta_tmp = TMath::ACosH(sqrt(Antitop_p3norm2)/Antitop_Pt_tmp);
    double Antitop_Pz_bis_tmp = Antitop_Pt_tmp*TMath::SinH(Antitop_Eta_tmp);
    if (TMath::Abs(Antitop_Pz_tmp-Antitop_Pz_bis_tmp)>1.e-3 && TMath::Abs(Antitop_Pz_tmp+Antitop_Pz_bis_tmp)>1.e-3) {
      if (verbosity>=2) cout <<" Antitop Pz="<<Antitop_Pz_tmp<<" or "<<Antitop_Pz_bis_tmp<<", flipping Top eta sign"<<endl;
      Top_Eta = -Top_Eta;
    }
    Top.SetPtEtaPhiE(Top_Pt, Top_Eta, Top_phi, Top_E);
    Antitop = Ptot - Higgs - Top;
    if (CheckMomentum(Antitop, Antitop_M)!=1) {
      if (verbosity>=2) cout << "First solution leads to wrong antitop mass, trying the other"<<endl;
      Top_Pt = X2;
    }
  }

  Top_Eta = TMath::ACosH(sqrt(Top_p3norm2)/Top_Pt);
  double Top_Pz = Top_Pt*TMath::SinH(Top_Eta);
  double Antitop_Pz = Ptot.Pz()-Higgs.Pz()-Top_Pz;
  double Antitop_Pt = sqrt(Higgs.Pt()*Higgs.Pt()+B*Top_Pt+Top_Pt*Top_Pt);
  double Antitop_Eta = TMath::ACosH(sqrt(Antitop_p3norm2)/Antitop_Pt);
  double Antitop_Pz_bis = Antitop_Pt*TMath::SinH(Antitop_Eta);
  if (TMath::Abs(Antitop_Pz-Antitop_Pz_bis)>1.e-3 && TMath::Abs(Antitop_Pz+Antitop_Pz_bis)>1.e-3 ) {
    if (verbosity>=2) cout <<" Antitop Pz="<<Antitop_Pz<<" or "<<Antitop_Pz_bis<<", flipping Top eta sign"<<endl;
    Top_Eta = -Top_Eta;
  }
  
  Top.SetPtEtaPhiE(Top_Pt, Top_Eta, Top_phi, Top_E);
  Antitop = Ptot - Higgs - Top;

  ReadPartonMomenta(pCore, 5);

  int check = 1;
  check =  CheckMomentum(Top, Top_M);
  check *= CheckMomentum(Antitop, Antitop_M);
  check *= CheckMomentum(Higgs, Higgs_M);
  if (check==0) {
    if (verbosity>=2) cout << "Mass constraint is not respected"<<endl;
    return 0;
  }

  FillTTHPhaseSpacePoint(Top, Antitop, Higgs);

  double weightPS_1 = Higgs.Beta()*Higgs.Vect().Mag()*sin(Higgs.Theta());
  double weightPS_3 = 1./(Ptot-Higgs).Vect().Mag();
  //double weightPS_2 = 1./Antitop.E();
  //double weightPS_3 = 2*TMath::Abs(Top.Pz())/Ptot.E();
  //double weight_PS = weightPS_1*weightPS_2*weightPS_3/8./TMath::Power(2*TMath::Pi(), 5);
  double weight_PS = weightPS_1*weightPS_3*1./8./TMath::Power(2*TMath::Pi(), 5);

  return weight_PS;
}

double MEPhaseSpace::SetupKinematicsTopHad(const double *x) const {

  //x[0] theta b-jet
  //x[1] phi b-jet
  //x[2] theta jet1
  //x[3] phi jet1
  //x[4] theta jet2
  //x[5] phi jet2
  //x[6] energy jet2

  double Jet2_Theta = x[4];
  double Jet2_Phi = x[5];
  double Jet2_E = x[6];
  double Jet2_Eta = -log(tan(Jet2_Theta/2.));
  double Jet2_Pt = Jet2_E / TMath::CosH(Jet2_Eta);
  TLorentzVector Jet2;
  Jet2.SetPtEtaPhiE(Jet2_Pt, Jet2_Eta, Jet2_Phi, Jet2_E);

  double Jet1_Phi = x[3];
  double Jet1_Theta = x[2];
  TVector3 Jet1unit;
  Jet1unit.SetMagThetaPhi(1, Jet1_Theta, Jet1_Phi);

  double Theta_Jet12 = Jet1unit.Angle(Jet2.Vect().Unit());
  double Jet1_E = mW * mW / ( 4* Jet2_E * sin(Theta_Jet12/2.) * sin(Theta_Jet12/2.) );
  double Jet1_Eta = -log(tan(Jet1_Theta/2));
  double Jet1_Pt = Jet1_E / TMath::CosH(Jet1_Eta);
  TLorentzVector Jet1;
  Jet1.SetPtEtaPhiE(Jet1_Pt, Jet1_Eta, Jet1_Phi, Jet1_E);

  TLorentzVector W = Jet1 + Jet2;
  if (verbosity>=2) cout << "Checking W momentum"<<endl;
  CheckMomentum(W, mW);

  double Bjet_Theta = x[0];
  double Bjet_Phi = x[1];
  TVector3 BjetUnit;
  BjetUnit.SetMagThetaPhi(1, Bjet_Theta, Bjet_Phi);

  double a = Jet1_E + Jet2_E;
  double b = Jet2_E * (BjetUnit.Dot(Jet2.Vect().Unit())) + Jet1_E * (Jet1unit.Dot(BjetUnit));
  double deltaMtop = (mTop*mTop - mW*mW - mB*mB) / 2.;
  double delta = deltaMtop*deltaMtop - (a*a-b*b)*mB*mB;
  if (delta<0){
    if (verbosity>=2) cout << "delta="<< delta<<" < 0, quit" << endl;
    return 0;
  }
  double X1 = (a*deltaMtop+TMath::Abs(b)*sqrt(delta))/(a*a-b*b);
  double X2 = (a*deltaMtop-TMath::Abs(b)*sqrt(delta))/(a*a-b*b);

  if (verbosity>=2) cout << "X1="<<X1<<" and X2="<< X2<<endl;

  double Bjet_E=0;
  double Bjet_Eta;
  double Bjet_Pt;
  TLorentzVector Bjet;
  TLorentzVector Top;
  if (X1<0 && X2<0) {
    if (verbosity>=2) cout << "X1<0 and X2<0, quit"<<endl;
	return 0;
  }
  if (X1<0 && X2>0) Bjet_E = X2;
  if (X1>0 && X2<0) Bjet_E = X1;
  if (X1>0 && X2>0){
    if (verbosity>=2) cout << "Both X1 and X2 positive, trying first value"<<endl;
    //if (X1<X2)
    //  Bjet_E = X2;
    //else 
    //  Bjet_E = X1;
    Bjet_E = X1;
    Bjet_Eta = -log(tan(Bjet_Theta/2.));
    Bjet_Pt = sqrt(Bjet_E*Bjet_E-mB*mB) / TMath::CosH(Bjet_Eta);
    Bjet.SetPtEtaPhiE(Bjet_Pt, Bjet_Eta, Bjet_Phi, Bjet_E);
    Top = W + Bjet;
    int res = CheckMomentum(Top, mTop);
    if (res==0) {
	if (verbosity>=2) cout << "First solution does not respect mass constraint, trying the other"<<endl;
	Bjet_E = X2;
    }
  }

  Bjet_Eta = -log(tan(Bjet_Theta/2.));
  Bjet_Pt = sqrt(Bjet_E*Bjet_E-mB*mB) / TMath::CosH(Bjet_Eta);
  Bjet.SetPtEtaPhiE(Bjet_Pt, Bjet_Eta, Bjet_Phi, Bjet_E);

  Top = W + Bjet;

  FillTopDecayPhaseSpacePoint(Top, Bjet, Jet1, Jet2, kTop);
  if (iMode==kAllPartonsTopLep) FillTopDecayPhaseSpacePoint(Top, Bjet, Jet2, Jet1, kTop);

  ReadPartonMomenta(pTop, 4);

  int check = 1;
  check =  CheckMomentum(Top, mTop);
  check *= CheckMomentum(Bjet, mB);
  check *= CheckMomentum(Jet1, 0);
  check *= CheckMomentum(Jet2, 0);

  if (check==0) {
    if (verbosity>=2) cout << "Mass constraint is not respected"<<endl;
    return 0;
  }

  double weight_PS1 = 1./(2.*Top.E()*2.*W.E());
  double weight_PS2 = Bjet.Vect().Mag() * sin(Bjet.Theta()) / 2.;
  double weight_PS3 = Jet1.Vect().Mag() * sin(Jet1.Theta()) / 2.;
  double weight_PS4 = Jet2.Vect().Mag() * sin(Jet2.Theta()) / 2.;
  double weight_PS = 1./TMath::Power(2*TMath::Pi(), 7) * weight_PS1 * weight_PS2 * weight_PS3 * weight_PS4;

  if (iMode==kAllPartonsTopHad || iMode==kAllPartonsTopLep) weight_PS *= 1./(2.*mTop); //decay width

  return weight_PS;
}

double MEPhaseSpace::SetupKinematicsTopHad_FiniteWidths(const double *x) const {

  //x[0] theta b-jet
  //x[1] phi b-jet
  //x[2] theta jet1
  //x[3] phi jet1
  //x[4] theta jet2
  //x[5] phi jet2
  //x[6] energy jet2
  //x[7] W energy
  //x[8] top energy
  
  double Jet2_Theta = x[4];
  double Jet2_Phi = x[5];
  double Jet2_E = x[6];
  double Jet2_Eta = -log(tan(Jet2_Theta/2.));
  double Jet2_Pt = Jet2_E / TMath::CosH(Jet2_Eta);
  TLorentzVector Jet2;
  Jet2.SetPtEtaPhiE(Jet2_Pt, Jet2_Eta, Jet2_Phi, Jet2_E);
  
  double Jet1_Phi = x[3];
  double Jet1_Theta = x[2];
  double W_E = x[7];
  double Jet1_E = W_E - Jet2_E;
  double Jet1_Eta = -log(tan(Jet1_Theta/2));
  double Jet1_Pt = Jet1_E / TMath::CosH(Jet1_Eta);
  if (Jet1_E<0){
	  if (verbosity>=2) cout << "Jet1 energy="<<Jet1_E<<" <0, quit"<<endl;
	  return 0;
  }
  
  TLorentzVector Jet1;
  Jet1.SetPtEtaPhiE(Jet1_Pt, Jet1_Eta, Jet1_Phi, Jet1_E);
  
  TLorentzVector W = Jet1 + Jet2;
  if (verbosity>=2) cout << "W mass2 ="<< W.Mag2()<<endl;
  if (W.Mag2()<0){
	  if (verbosity>=2) cout << "W mass2 ="<< W.Mag2()<<" <0, quit"<<endl;
	  return 0;
  }
  
  double Bjet_Theta = x[0];
  double Bjet_Phi = x[1];
  double Top_E = x[8];
  double Bjet_E = Top_E - W_E;
  double Bjet_Eta = -log(tan(Bjet_Theta/2));
  if (Bjet_E<mB){
	  if (verbosity>=2) cout << "Bjet energy="<<Bjet_E<<" <mB, quit"<<endl;
	  return 0;
  }
  double Bjet_Pt = sqrt(Bjet_E*Bjet_E-mB*mB) / TMath::CosH(Bjet_Eta);
  TLorentzVector Bjet;
  Bjet.SetPtEtaPhiE(Bjet_Pt, Bjet_Eta, Bjet_Phi, Bjet_E);

  TLorentzVector Top;
  Top = W + Bjet;
  if (verbosity>=2) cout << "Top mass2 ="<< Top.Mag2()<<endl;
if (Top.Mag2()<0){
	  if (verbosity>=2) cout << "Top mass2 ="<< Top.Mag2()<<" <0, quit"<<endl;
	  return 0;
  }
  

/*
  double Jet2_Theta = x[4];
  double Jet2_Phi = x[5];
  double Jet2_E = x[6];
  double Jet2_Eta = -log(tan(Jet2_Theta/2.));
  double Jet2_Pt = Jet2_E / TMath::CosH(Jet2_Eta);
  TLorentzVector Jet2;
  Jet2.SetPtEtaPhiE(Jet2_Pt, Jet2_Eta, Jet2_Phi, Jet2_E);

  double Jet1_Phi = x[3];
  double Jet1_Theta = x[2];
  TVector3 Jet1unit;
  Jet1unit.SetMagThetaPhi(1, Jet1_Theta, Jet1_Phi);

  double W_M = x[7]; //different from mW
  double Theta_Jet12 = Jet1unit.Angle(Jet2.Vect().Unit());
  double Jet1_E = W_M * W_M / ( 4* Jet2_E * sin(Theta_Jet12/2.) * sin(Theta_Jet12/2.) );
  double Jet1_Eta = -log(tan(Jet1_Theta/2));
  double Jet1_Pt = Jet1_E / TMath::CosH(Jet1_Eta);
  TLorentzVector Jet1;
  Jet1.SetPtEtaPhiE(Jet1_Pt, Jet1_Eta, Jet1_Phi, Jet1_E);

  TLorentzVector W = Jet1 + Jet2;
  if (verbosity>=2) cout << "Checking W momentum"<<endl;
  CheckMomentum(W, W_M);

  double Bjet_Theta = x[0];
  double Bjet_Phi = x[1];
  TVector3 BjetUnit;
  BjetUnit.SetMagThetaPhi(1, Bjet_Theta, Bjet_Phi);

  double Top_M = x[8];
  double a = Jet1_E + Jet2_E;
  double b = Jet2_E * (BjetUnit.Dot(Jet2.Vect().Unit())) + Jet1_E * (Jet1unit.Dot(BjetUnit));
  double deltaMtop = (Top_M*Top_M - W_M*W_M - mB*mB) / 2.;
  double delta = deltaMtop*deltaMtop - (a*a-b*b)*mB*mB;
  if (delta<0){
    if (verbosity>=2) cout << "delta="<< delta<<" < 0, quit" << endl;
    return 0;
  }
  double X1 = (a*deltaMtop+TMath::Abs(b)*sqrt(delta))/(a*a-b*b);
  double X2 = (a*deltaMtop-TMath::Abs(b)*sqrt(delta))/(a*a-b*b);

  if (verbosity>=2) cout << "X1="<<X1<<" and X2="<< X2<<endl;

  double Bjet_E=0;
  double Bjet_Eta;
  double Bjet_Pt;
  double sgn = 0;
  TLorentzVector Bjet;
  TLorentzVector Top;
  if (X1<0 && X2<0) {
    if (verbosity>=2) cout << "X1<0 and X2<0, quit"<<endl;
	return 0;
  }
  if (X1<0 && X2>0) {
	  sgn = -1;
	  Bjet_E = X2;
  }
  if (X1>0 && X2<0) {
	  sgn = 1;
	  Bjet_E = X1;  
  }
  if (X1>0 && X2>0){
    if (verbosity>=2) cout << "Both X1 and X2 positive, trying first value"<<endl;
    //if (X1<X2)
    //  Bjet_E = X2;
    //else 
    //  Bjet_E = X1;
    Bjet_E = X1;
	sgn = 1;
    Bjet_Eta = -log(tan(Bjet_Theta/2.));
    Bjet_Pt = sqrt(Bjet_E*Bjet_E-mB*mB) / TMath::CosH(Bjet_Eta);
    Bjet.SetPtEtaPhiE(Bjet_Pt, Bjet_Eta, Bjet_Phi, Bjet_E);
    Top = W + Bjet;
    int res = CheckMomentum(Top, Top_M);
    if (res==0) {
	if (verbosity>=2) cout << "First solution does not respect mass constraint, trying the other"<<endl;
	Bjet_E = X2;
	sgn = -1;
    }
  }

  Bjet_Eta = -log(tan(Bjet_Theta/2.));
  Bjet_Pt = sqrt(Bjet_E*Bjet_E-mB*mB) / TMath::CosH(Bjet_Eta);
  Bjet.SetPtEtaPhiE(Bjet_Pt, Bjet_Eta, Bjet_Phi, Bjet_E);

  Top = W + Bjet;
*/
  FillTopDecayPhaseSpacePoint(Top, Bjet, Jet1, Jet2, kTop);

  ReadPartonMomenta(pTop, 4);

  int check = 1;
  check =  CheckMomentum(Top, Top.M());
  check *= CheckMomentum(Bjet, mB);
  check *= CheckMomentum(Jet1, 0);
  check *= CheckMomentum(Jet2, 0);

  if (check==0) {
    if (verbosity>=2) cout << "Mass constraint is not respected"<<endl;
    return 0;
  }
  
/*
  //Jacobian Et, Ew --> mT, mW
  double J11=0, J12=0, J21=0, J22=0;
  J11 = Top_M / Top.E();
  J22 = W_M / W.E();
  
  if (verbosity>=2) cout << "J11="<<J11<<" J22="<<J22<<endl;
  
  double ap = 2*Jet1_E / W_M;
  double bp = ap * Jet1unit.Dot(BjetUnit);
  double a2mb2p = 2*ap*a-2*bp*b;
  double deltaMtopp = -W_M;
  double Part1 = (ap*(a*a-b*b)*deltaMtop - a*a2mb2p*deltaMtop + a*(a*a-b*b)*deltaMtopp) / ((a*a-b*b)*(a*a-b*b));
  double deltap = 2*deltaMtopp*deltaMtop - a2mb2p*mB*mB;
  double u=sgn*TMath::Abs(b); 
  double up = bp;
  double v=a*a-b*b;
  double vp = a2mb2p; 
  double w=sqrt(delta);
  double wp = deltap / (2*w);
  double Part2 = (up*v*w - u*vp*w + u*v*wp) / (v*v);
  J21 = W_M / W.E() + Part1 + Part2;
  
 if (verbosity>=2)  cout << "J21="<<J21<<endl;
  
  J12 = Top_M / Top.E() - ( a*Top_M/(a*a-b*b)  + sgn*TMath::Abs(b)/(a*a-b*b)*Top_M*deltaMtop/(2*sqrt(delta)));
  
  if (verbosity>=2) cout << "J12="<<J12<<endl;
  
  double detJ = TMath::Abs(J11*J22 - J12*J21);
  
  if (verbosity>=2) cout << "detJ="<< detJ<<endl;
  
  double weight_PS1 = detJ;
*/

  double gammaTop = 1.491500;
  double gammaW = 2.047600;
  
  double weight_PS1 = BreitWigner(Top_E, mTop, gammaTop) * BreitWigner(W_E, mW, gammaW);
  //double weight_PS1 = 1/(mTop*mTop) * 1/(mW*mW);  

  double weight_PS2 = Bjet.Vect().Mag() * sin(Bjet.Theta()) / 2.;
  double weight_PS3 = Jet1.Vect().Mag() * sin(Jet1.Theta()) / 2.;
  double weight_PS4 = Jet2.Vect().Mag() * sin(Jet2.Theta()) / 2.;
  double weight_PS = 1./TMath::Power(2*TMath::Pi(), 9) * weight_PS1 * weight_PS2 * weight_PS3 * weight_PS4;

  if (verbosity>=2) cout << "weight_PS="<< weight_PS << endl;

  if (iMode==kAllPartonsTopHad_FiniteWidth) weight_PS *= 1./(2.*mTop); //decay width

  return weight_PS;
}

double MEPhaseSpace::BreitWigner(double E, double M, double Gamma) const{

	double g = sqrt(M*M*(M*M+Gamma*Gamma));
	double k = 2*sqrt(2)*M*Gamma*g/(TMath::Pi()*sqrt(M*M+g));
	
	double BW = k/((E*E-M*M)*(E*E-M*M)+M*M*Gamma*Gamma);
	
	return BW;
}

double MEPhaseSpace::KallenFunction(double Ecm, double m1, double m2) const{

  double K = (Ecm*Ecm-m1*m1-m2*m2)*(Ecm*Ecm-m1*m1-m2*m2)-4*m1*m1*m2*m2;

  return K;
}


double MEPhaseSpace::SetupKinematicsTopLep(const double *x) const {

  //x[0] theta b-jet
  //x[1] phi b-jet
  //x[2] theta neutrino
  //x[3] phi neutrino
  //x[4] theta lepton
  //x[5] phi lepton
  //x[6] energy lepton

  double weightPS = SetupKinematicsTopHad(x);

  return weightPS;
}

double MEPhaseSpace::SetupKinematicsTopLep_FixedTopM_FiniteWwidth(const double *x) const {

  //x[0] theta neutrino
  //x[1] phi neutrino
  //x[2] theta b-jet
  //x[3] phi b-jet
  //x[4] W mass

  double Bjet_M = mB;

  TLorentzVector Top;
  Top.SetPxPyPzE(0,0,0,mTop);

  double W_M = x[4];
  double W_P = sqrt(KallenFunction(Top.E(), W_M, Bjet_M)) / (2*Top.E());
  double W_E = sqrt(W_M*W_M+W_P*W_P);

  //double W_E = Wenergy;//x[4];

  TLorentzVector Bjet;
  double Bjet_E = Top.E() - W_E;

  double Bjet_Phi = x[3];
  double Bjet_Theta = x[2];
  double Bjet_Eta = -log(tan(Bjet_Theta/2));
  if (Bjet_E<Bjet_M) {
    if (verbosity>=2) cout << "Bjet E<M, quit"<<endl;
    return 0;
  }
  double Bjet_Pt = sqrt(Bjet_E*Bjet_E-mB*mB) / TMath::CosH(Bjet_Eta);
  Bjet.SetPtEtaPhiE(Bjet_Pt, Bjet_Eta, Bjet_Phi, Bjet_E);

  //Bjet.SetPxPyPzE(0,0,W_P,Bjet_E);
  
  TLorentzVector W = Top - Bjet;
  if (verbosity>=2) {
    cout << "W mass="<<W.M()<<endl;
    cout << "DeltaPhi="<<(W.Phi()-Bjet.Phi())<<" DeltaTheta="<<W.Theta()-Bjet.Theta()<<endl;

  }
  if (W.E()<0 || W.Mag2()<0) {
    if (verbosity>=2) cout << "W energy="<<W.E()<<" P2="<<W.Mag2()<<", quit"<<endl;
    return 0;
  }

  //W decay in its cm frame
  double Neut_P = sqrt(KallenFunction(W.M(), 0, 0)) / (2*W.M()) ;
  double Neut_Theta = x[0];
  double Neut_Eta = -log(tan(Neut_Theta/2));
  double Neut_Phi = x[1];
  double Neut_Pt = Neut_P / TMath::CosH(Neut_Eta);
//  double Neut_Pt = W.M()*W.M()/2./(W.E()*TMath::CosH(Neut_Eta) - W.Px()*cos(Neut_Phi) - W.Py()*sin(Neut_Phi) - W.Pz()*TMath::SinH(Neut_Eta));
  if (Neut_Pt<0) {
    if (verbosity>=2) cout << "Neutrino Pt<0, quit"<<endl;
    return 0;
  }
  //double Neut_E = Neut_Pt * TMath::CosH(Neut_Eta);
  TLorentzVector Neut;
  Neut.SetPtEtaPhiE(Neut_Pt, Neut_Eta, Neut_Phi, Neut_P);

  TLorentzVector Wcm; Wcm.SetPxPyPzE(0,0,0,W.M());
  TLorentzVector Lepton= Wcm - Neut;

  //Boost back to lab frame
  Neut.Boost(W.BoostVector());
  Lepton.Boost(W.BoostVector());

  FillTopDecayPhaseSpacePoint(Top, Bjet, Lepton, Neut, kTop);

  ReadPartonMomenta(pTop, 4);

  int check = 1;
  check =  CheckMomentum(Top, mTop);
  check *= CheckMomentum(Bjet, mB);
  check *= CheckMomentum(Lepton, 0);
  check *= CheckMomentum(Neut, 0);

  if (check==0) {
    if (verbosity>=2) cout << "Mass constraint is not respected"<<endl;
    return 0;
  }

/*
  if (verbosity>=2) {
    cout << "Checking Kallen functions, |Pb|="<<Bjet.Vect().Mag()<<" via Kallen formula "<<
W_P << endl;
    double Pnu = Neut.Vect().Mag();
    Neut.Boost(-W.BoostVector());
    cout << "Checking Kallen functions, |Pnu|="<<Pnu<<" via Kallen formula "<< sqrt(KallenFunction(W.M(), 0, 0)) / (2*W.M()) << " by inverse boost: "<< Neut.Vect().Mag()<<endl;
  }
*/

  double weight_PS1 = 2*W.M()/(16*TMath::Power(2*TMath::Pi(), 5));
  double weight_PS2 = W_P*sin(Bjet_Theta)/ Top.M();
  double weight_PS3 = Neut_P*sin(Neut_Theta)/ W.M();
  //double weight_PS4 = 1./(2.*Lepton.E());
  double weight_PS = weight_PS1*weight_PS2*weight_PS3;//*weight_PS4;

  weight_PS /= (mTop);

  return weight_PS; //ME factor
}

double MEPhaseSpace::SetupKinematicsTopHad_FixedTopM_FiniteWwidth(const double *x) const {

  double weight_PS = SetupKinematicsTopLep_FixedTopM_FiniteWwidth(x);
  return weight_PS;
}

double MEPhaseSpace::SetupKinematics1to3_LabFrame(const double *x) const {

  //x[0] E1
  //x[1] theta1
  //x[2] phi1
  //x[3] theta3
  //x[4] phi3

  if (verbosity>=2){
    cout << "E1="<<x[0]<<endl;
    cout << "Theta1="<<x[1]<<endl;
    cout << "Phi1="<<x[2]<<endl;
    cout << "Theta3="<<x[3]<<endl;
    cout << "Phi3="<<x[4]<<endl;
    }

  //Ptot top
  //P1 b-jet
  //P2 jet1
  //P3 jet2
  
  //TLorentzVector Ptot(1000,1000,1000,sqrt(mTop*mTop+3*1000*1000)); //should also work in all frame
  TLorentzVector Ptot;
  double P1_M=0;
  double P2_M=0;
  double P3_M=0;
  if (iOption==kTwoBjorken1to3){
    Ptot.SetPxPyPzE(0,0,(x[5]-x[6])*comEnergy/2., (x[5]+x[6])*comEnergy/2.);
    P1_M = mTop;
    P2_M = mTop;
    P3_M = mHiggs;
  }
  else {
    Ptot.SetPxPyPzE(0,0,0,mTop);
    P1_M = mB;
    P2_M = 0;
    P3_M = 0;
  }

  //Could also work for TTH instead of top decay

  //Reconstruct P1 (angles not necessarily defined relative to Ptot)
  double P1_Theta = x[1];
  double P1_Phi = x[2];
  double P1_E = x[0];
  if (P1_E<P1_M  || P1_E > Ptot.E()){
    if (verbosity>=2) cout << "P1 energy out of range, quit"<<endl;
    return 0;
  } 
  double P1_Eta = -log(tan(P1_Theta/2.));
  double P1_Pt = sqrt(-P1_M*P1_M+P1_E*P1_E) / TMath::CosH(P1_Eta);
  TLorentzVector P1;
  P1.SetPtEtaPhiE(P1_Pt, P1_Eta, P1_Phi, P1_E);
 
  //P3 angles.
  double P3_Theta = x[3];
  double P3_Phi = x[4];
  TVector3 P3unit;
  P3unit.SetPtThetaPhi(1, P3_Theta, P3_Phi);
  P3unit.SetMag(1);
  double cosAlpha = (((Ptot-P1).Vect()).Unit()).Dot(P3unit);
  if (verbosity>=2) {
    cout << "cosAlpha="<<cosAlpha << endl;
    cout << "Avec P3unit x="<<P3unit.X()<< " y="<<P3unit.Y()<<" z="<<P3unit.Z()<<endl;
    cout << "From P3unit, theta3="<<P3unit.Theta()<<" phi3="<<P3unit.Phi()<<endl;
  }
  double P3_Mag = 0;
  double P3_Eta = -log(tan(P3_Theta/2.));
  double P3_Pt = 0;
  TLorentzVector P3;
  TLorentzVector P2 ;
  P2.Clear(); P3.Clear();

  //Find P3 momentum
  double a = -P2_M*P2_M + Ptot.M2() + P1.M2() + P3_M*P3_M -2*Ptot.Dot(P1);
  double b = 2*(Ptot-P1).Vect().Mag()*cosAlpha ; 
  double c = -2*(Ptot.E()-P1.E());
  double A = b*b-c*c;
  double B = 2*a*b;
  double C = a*a-P3_M*P3_M*c*c;
  double delta = B*B-4*A*C;
  if (delta<0){
    if (verbosity>=2) cout << "delta<0, quit"<<endl;
    return 0;
  }
  double X1 = (-B+sqrt(delta))/(2*A);
  double X2 = (-B-sqrt(delta))/(2*A);
  if (verbosity>=2) cout << "X1="<<X1<<" X2="<<X2<<endl;
  if (X1<0 && X2<0){
    if (verbosity>=2) cout << "Both solutions are negative, quit" << endl;
    return 0;
  }
  if (X1>0 && X2<0){
    P3_Mag = X1;
  }
  else if (X1<0 && X2>0){
    P3_Mag = X2;
  }
  else if (X1>0 && X2>0){
    if (verbosity>=2) cout << "Both solutions are positive, try the first one" << endl;
    P3_Mag = X1;
    //P3unit.SetMag(P3_Mag);
    P3_Pt = P3_Mag / TMath::CosH(P3_Eta);
    P3.SetPtEtaPhiE(P3_Pt, P3_Eta, P3_Phi, sqrt(P3_M*P3_M+P3_Mag*P3_Mag));
    P2 = Ptot - P1 - P3;
    if (verbosity>=2) cout << "Test f(P3_Mag)"<< a + b*P3_Mag + c*sqrt(P3_Mag*P3_Mag+P3_M*P3_M) << endl;
    int res = CheckMomentum(P2, P2_M); 
    if (P2.E()>Ptot.E()-P1.E() || P2.E()<0 || P3.E() > Ptot.E()-P1.E() || P3.E()<0){
      if (verbosity>=2) cout << "Energy out of range"<<endl;
      res=0;
    }
    if (res!=1) { 
      P3_Mag = X2;
      P2.Clear(); P3.Clear();
      if (verbosity>=2) cout << "Try the other one"<<endl;
    }
  }
  //P3unit.SetMag(P3_Mag); P3 = P3unit;
  P3_Pt = P3_Mag / TMath::CosH(P3_Eta); 
  P3.SetPtEtaPhiE(P3_Pt, P3_Eta, P3_Phi, sqrt(P3_M*P3_M+P3_Mag*P3_Mag));

  P2 = Ptot - P1 - P3;

  if (P2.E()>Ptot.E()-P1.E() || P2.E()<0 || P3.E() > Ptot.E()-P1.E() || P3.E()<0){
    if (verbosity>=2) cout << "Energy out of range, quit"<<endl;
    return 0;
  }

  //Test
  if (verbosity>=2) cout << "Test f(P3_Mag)="<< a + b*P3_Mag + c*sqrt(P3_Mag*P3_Mag+P3_M*P3_M) << endl;
 
/* 
  double WidthW = 2.047600e+00;
  if (!((P2+P3).M()<mW+15*WidthW && (P2+P3).M()>mW-15*WidthW)) {
    if (verbosity>=2) cout << "W mass out of range"<<endl;
    //return 0;
  }
*/

  if (iOption!=kTwoBjorken1to3) FillTopDecayPhaseSpacePoint(Ptot, P1, P2, P3, kTop);
  else FillTTHPhaseSpacePoint(P1,P2,P3); 
  if (iMode==kAllPartonsTTH) ReadPartonMomenta(pCore, 5);
  else ReadPartonMomenta(pTop, 4);

  int check = 1;
  check =  CheckMomentum(Ptot, Ptot.M());
  check *= CheckMomentum(P1, P1_M);
  check *= CheckMomentum(P2, P2_M);
  check *= CheckMomentum(P3, P3_M);

  if (check==0) {
    if (verbosity>=2) cout << "Mass constraint is not respected"<<endl;
    return 0;
  }

  double dfdP3 = TMath::Abs(b + c*P3.Vect().Mag()/P3.E());
  double weightPS1 = P1.Beta()*P1.Vect().Mag()*sin(P1_Theta)/2.;
  double weightPS2 = 1/(TMath::Power(2*TMath::Pi(), 5));
  double weightPS3 = P3.Beta()*P3.Vect().Mag2()*sin(P3_Theta)/(2*P3.E())/dfdP3;
  double weightPS = weightPS1*weightPS2*weightPS3;

  if (!(iMode==kAllPartonsTTH)) weightPS /= (2*Ptot.M());
  
  return weightPS * 2;
}

double MEPhaseSpace::SetupKinematicsHiggsWWLep_FixedHiggsM_FiniteWwidth(const double *x) const {

  //x[0] theta W rest frame
  //x[1] phi W rest frame
  //x[2] theta lep1 rest frame
  //x[3] phi lep1 rest frame
  //x[4] theta lep2 rest frame
  //x[5] phi lep2 rest frame
  //x[6] W1 mass
  //x[7] W2 mass

  //Sans integrer sur phi:
  //x[0] theta W rest frame
  //x[1] theta lep1 rest frame
  //x[2] theta lep2 rest frame 
  //x[3] W1 mass
  //x[4] W2 mass 

  TLorentzVector Higgs(0,0,0,mHiggs);
  TLorentzVector W1;
  TLorentzVector W2;
  double W1_M = x[3];//x[6];
  double W2_M = x[4];//x[7];
  double W_theta = x[0];
  double W_phi = 0;//x[1];
  double weight_PS_HWW = ComputeDecayMomenta(Higgs, W1_M, W2_M, W_theta, W_phi, &W1, &W2);
  if (weight_PS_HWW==0) return 0;

  TLorentzVector Lepton1; 
  TLorentzVector Neutrino1;
  double Lepton1_Theta = x[1];//x[2];
  double Lepton1_Phi = 0;//x[3];
  double weight_PS_Wlnu1 = ComputeDecayMomenta(W1, 0, 0, Lepton1_Theta, Lepton1_Phi, &Lepton1, &Neutrino1);
  if (weight_PS_Wlnu1==0) return 0;

  TLorentzVector Lepton2;
  TLorentzVector Neutrino2;
  double Lepton2_Theta = x[2];//x[4];
  double Lepton2_Phi = 0;//x[5];
  double weight_PS_Wlnu2 = ComputeDecayMomenta(W2, 0, 0, Lepton2_Theta, Lepton2_Phi, &Lepton2, &Neutrino2);
  if (weight_PS_Wlnu2==0) return 0;

  FillHiggsDecayPhaseSpacePoint(Higgs, Lepton1, Neutrino1, Lepton2, Neutrino2);
  ReadPartonMomenta(pHiggs, 5);
  int check = 1;
  check =  CheckMomentum(Higgs, mHiggs);
  check *= CheckMomentum(Lepton1, 0);
  check *= CheckMomentum(Neutrino1, 0);
  check *= CheckMomentum(Lepton2, 0);
  check *= CheckMomentum(Neutrino2, 0);
  if (check==0) {
    if (verbosity>=2) cout << "Mass constraint is not respected"<<endl;
    return 0;
  }

  double weight_PS = 4 * W1_M * W2_M * weight_PS_HWW * weight_PS_Wlnu1 * weight_PS_Wlnu2 / TMath::Power(2*TMath::Pi(), 2);

  weight_PS *= TMath::Power(2*TMath::Pi(), 3);  //no integration over the three phi
 
  weight_PS /= (2*mHiggs);

  return weight_PS ; //ME factor
}

double MEPhaseSpace::SetupKinematics2to3_LabFrame_OneBjorken(const double* x) const{

  //x[0] bjorken x1
  //x[1] Energy Top
  //x[2] Theta Top
  //x[3] Energy Antitop
  //x[4] Theta Antitop
  //x[5] Phi Top
  //x[6] Phi Antitop
  //Integrating out the top and antitop phi
  
  double P1_M = mTop; //Top
  double P2_M = mTop; //Antitop
  double P3_M = mHiggs;  //Higgs   

  TLorentzVector P1; //Top
  SetMomentumFromEThetaPhi(P1_M, x[1], x[2], x[5], &P1);
  TLorentzVector P2; //AntiTop
  SetMomentumFromEThetaPhi(P2_M, x[3], x[4], x[6], &P2);
   
  //Compute x2
  double x1 = x[0];
  double a = -P3_M*P3_M +P1.M2() +P2.M2() +2*P1.Dot(P2) -x1*comEnergy*(P1.E()+P2.E()-P1.Pz()-P2.Pz());
  double b = x1*comEnergy*comEnergy -comEnergy*(P1.E()+P2.E()+P1.Pz()+P2.Pz());
  double x2 = -a/b;

  TLorentzVector Ptot(0, 0, (x1-x2)*comEnergy/2., (x1+x2)*comEnergy/2.);
  TLorentzVector P3 = Ptot -P1 -P2;

  if (!(x2>0 && x2<1)){
    if (verbosity>=2) cout << "x2="<<x2<<"<0, quit"<<endl;
    return 0;
  }
  if (!(P1.E()<Ptot.E() && P2.E()<Ptot.E() && P3.E()<Ptot.E() && P3.E()>P3_M && P3.E()<Ptot.E())){
    if (verbosity>=2) cout << "Energy out of range"<<endl;
    return 0;
  }

  FillTTHPhaseSpacePoint(P1, P2, P3);
  SetInitialPartonMomenta(x1, x2);

  ReadPartonMomenta(pCore, 5);
  int check = 1;
  check =  CheckMomentum(P1, P1_M);
  check *= CheckMomentum(P2, P2_M);
  check *= CheckMomentum(P3, P3_M);
  if (check==0) {
    if (verbosity>=2) cout << "Mass constraint is not respected"<<endl;
    return 0;
  }

  double dfdx2 = TMath::Abs(b);
  double weight_PS1 = 1./dfdx2 /TMath::Power(2*TMath::Pi(), 5);
  double weight_PS2 = P1.Beta()*P1.Vect().Mag()*sin(P1.Theta())/2.;
  double weight_PS3 = P2.Beta()*P2.Vect().Mag()*sin(P2.Theta())/2.;
  double weight_PS = weight_PS1 * weight_PS2 * weight_PS3;

  //weight_PS *= TMath::Power(2*TMath::Pi(), 2); //integration over phi

  return weight_PS; //ME factor
}

double MEPhaseSpace::SetupKinematics2to3_LabFrame_NoBjorken(const double* x) const{

  //Angles integration
  //x[0] Higgs Pz
  //x[1] Energy Top
  //x[2] Theta Top
  //x[3] Phi Top
  //x[4] Energy antitop
  //x[5] Theta antitop
  //x[6] Phi Antitop

  //Momentum integration
  //x[0] Higgs Pz
  //x[1] Px Top
  //x[2] Py Top
  //x[3] Pz Top
  //x[4] Px antitop
  //x[5] Py antitop
  //x[6] Pz Antitop

  double P1_M = mTop;
  double P2_M = mTop;
  double P3_M = mHiggs;

  TLorentzVector P1; //Top
  TLorentzVector P2; //AntiTop
  TLorentzVector P3; //Higgs

  if (iOption==kNoBjorken){
    SetMomentumFromEThetaPhi(P1_M, x[1], x[2], x[3], &P1);
    SetMomentumFromEThetaPhi(P2_M, x[4], x[5], x[6], &P2);
  }
  else if (iOption==kNoBjorkenMomentum){
    double P1_E = sqrt(x[1]*x[1]+x[2]*x[2]+x[3]*x[3]+P1_M*P1_M);
    double P2_E = sqrt(x[4]*x[4]+x[5]*x[5]+x[6]*x[6]+P2_M*P2_M);
    P1.SetPxPyPzE(x[1],x[2],x[3],P1_E);
    P2.SetPxPyPzE(x[4],x[5],x[6],P2_E);
  }

    double P3_Px = -P1.Px() - P2.Px();
    double P3_Py = -P1.Py() - P2.Py();
    double P3_E = sqrt(P3_Px*P3_Px+P3_Py*P3_Py+x[0]*x[0]+P3_M*P3_M);
    P3.SetPxPyPzE(P3_Px, P3_Py, x[0], P3_E);

  TLorentzVector Ptot = P1 + P2 + P3;
  double x1 = (Ptot.Pz()+Ptot.E())/comEnergy;
  double x2 = (-Ptot.Pz()+Ptot.E())/comEnergy;

  if (verbosity>=2) cout << "x1+x2="<<x1+x2<<endl;

  if (!(x2>0 && x2<1)){
    if (verbosity>=2) cout << "x2="<<x2<<" out of range, quit"<<endl;
    return 0;
  }
  if (!(x1>0 && x1<1)){
    if (verbosity>=2) cout << "x1="<<x1<<" out of range, quit"<<endl;
    return 0;
  }

  FillTTHPhaseSpacePoint(P1, P2, P3);
  SetInitialPartonMomenta(x1, x2);

  ReadPartonMomenta(pCore, 5);
  int check = 1;
  check =  CheckMomentum(P1, P1_M);
  check *= CheckMomentum(P2, P2_M);
  check *= CheckMomentum(P3, P3_M);
  if (check==0) {
    if (verbosity>=2) cout << "Mass constraint is not respected"<<endl;
    return 0;
  }

  double weight_PS1 = 2./(TMath::Power(2*TMath::Pi(), 5))/(comEnergy*comEnergy);
  double weight_PS2 = 1/(2*P3.E());
  double weight_PS3 = P1.Beta()*P1.Vect().Mag()*sin(P1.Theta())/2;
  double weight_PS4 = P2.Beta()*P2.Vect().Mag()*sin(P2.Theta())/2;
  double weight_PS = 0;
  if (iOption==kNoBjorken) weight_PS = weight_PS1 * weight_PS2 * weight_PS3 * weight_PS4;
  if (iOption==kNoBjorkenMomentum) weight_PS = 1./TMath::Power(2*TMath::Pi(), 5) /(8*comEnergy*comEnergy*P1.E()*P2.E()*P3.E());

  return weight_PS; //ME factor
}

double MEPhaseSpace::SetupKinematics2to4_LabFrame_NoBjorken(const double* x) const{

  //Angles integration
  //x[0] Lepton1 Pz
  //x[1] Energy Top
  //x[2] Theta Top
  //x[3] Phi Top
  //x[4] Energy antitop
  //x[5] Theta antitop
  //x[6] Phi Antitop
  //x[7] Energy Lep2
  //x[8] Theta Lep2
  //x[9] Phi Lep2

  double P1_M = mTop;
  double P2_M = mTop;
  double P3_M = 0;
  double P4_M = 0;

  TLorentzVector P1; //Top
  TLorentzVector P2; //AntiTop
  TLorentzVector P3; //Lep1
  TLorentzVector P4; //Lep2

  SetMomentumFromEThetaPhi(P1_M, x[1], x[2], x[3], &P1);
  SetMomentumFromEThetaPhi(P2_M, x[4], x[5], x[6], &P2);
  SetMomentumFromEThetaPhi(P4_M, x[7], x[8], x[9], &P4);

  double P3_Px = -P1.Px() - P2.Px() - P4.Px();
  double P3_Py = -P1.Py() - P2.Py() - P4.Py();
  double P3_E = sqrt(P3_Px*P3_Px+P3_Py*P3_Py+x[0]*x[0]+P3_M*P3_M);
  P3.SetPxPyPzE(P3_Px, P3_Py, x[0], P3_E);

  TLorentzVector Ptot = P1 + P2 + P3 + P4;
  double x1 = (Ptot.Pz()+Ptot.E())/comEnergy;
  double x2 = (-Ptot.Pz()+Ptot.E())/comEnergy;

  if (verbosity>=2) cout << "x1+x2="<<x1+x2<<endl;

  if (!(x2>0 && x2<1)){
    if (verbosity>=2) cout << "x2="<<x2<<" out of range, quit"<<endl;
    return 0;
  }
  if (!(x1>0 && x1<1)){
    if (verbosity>=2) cout << "x1="<<x1<<" out of range, quit"<<endl;
    return 0;
  }


  FillTTLLPhaseSpacePoint(P1, P2, P3, P4);
  SetInitialPartonMomenta(x1, x2);

  ReadPartonMomenta(pCore, 6);
  int check = 1;
  check =  CheckMomentum(P1, P1_M);
  check *= CheckMomentum(P2, P2_M);
  check *= CheckMomentum(P3, P3_M);
  check *= CheckMomentum(P4, P4_M);
  if (check==0) {
    if (verbosity>=2) cout << "Mass constraint is not respected"<<endl;
    return 0;
  }

  double weight_PS1 = 2./(TMath::Power(2*TMath::Pi(), 8))/(comEnergy*comEnergy);
  double weight_PS2 = 1/(2*P3.E());
  double weight_PS3 = P1.Beta()*P1.Vect().Mag()*sin(P1.Theta())/2;
  double weight_PS4 = P2.Beta()*P2.Vect().Mag()*sin(P2.Theta())/2;
  double weight_PS5 = P4.Beta()*P4.Vect().Mag()*sin(P4.Theta())/2;
  double weight_PS = weight_PS1 * weight_PS2 * weight_PS3 * weight_PS4 * weight_PS5;

  return weight_PS * 4; //ME factor
}

double MEPhaseSpace::SetupKinematicsTTH_NoBjorken_TopHadDecay(const double* x) const{

  //Not decaying TTH variables
  //x[0] Higgs Pz
  //x[1] Energy antitop
  //x[2] Theta antitop
  //x[3] Phi Antitop

  //Top decay phase space
  //x[4] Energy Bjet
  //x[5] Theta Bjet
  //x[6] Phi Bjet
  //x[7] Energy Jet1
  //x[8] Theta Jet1
  //x[9] Phi Jet1
  //x[10] Theta Jet2
  //x[11] Phi Jet2 

  //First reconstruct the top from its decay
  y[0] = x[4];
  y[1] = x[5];
  y[2] = x[6];
  y[3] = x[7];
  y[4] = x[8];
  y[5] = x[9];
  y[6] = x[10];
  y[7] = x[11];

  TLorentzVector P1;
  double weight_PS_Top = SetupKinematics_TopHadDecay_WithTopPhaseSpace(y, &P1, kTop, kTopHadDecay);
  if (weight_PS_Top==0) return 0;

  double P1_M = mTop;
  double P2_M = mTop;
  double P3_M = mHiggs;

  TLorentzVector P2; //AntiTop
  TLorentzVector P3; //Higgs

  SetMomentumFromEThetaPhi(P2_M, x[1], x[2], x[3], &P2);
  double P3_Px = -P1.Px() - P2.Px();
  double P3_Py = -P1.Py() - P2.Py();
  double P3_E = sqrt(P3_Px*P3_Px+P3_Py*P3_Py+x[0]*x[0]+P3_M*P3_M);
  P3.SetPxPyPzE(P3_Px, P3_Py, x[0], P3_E);

  TLorentzVector Ptot = P1 + P2 + P3;
  double x1 = (Ptot.Pz()+Ptot.E())/comEnergy;
  double x2 = (-Ptot.Pz()+Ptot.E())/comEnergy;

  if (verbosity>=2) cout << "x1+x2="<<x1+x2<<endl;
  if (!(x2>0 && x2<1)){
    if (verbosity>=2) cout << "x2="<<x2<<" out of range, quit"<<endl;
    return 0;
  }
  if (!(x1>0 && x1<1)){
    if (verbosity>=2) cout << "x1="<<x1<<" out of range, quit"<<endl;
    return 0;
  }

  FillTTHPhaseSpacePoint(P1, P2, P3);
  SetInitialPartonMomenta(x1, x2);

  ReadPartonMomenta(pCore, 5);
  int check = 1;
  check =  CheckMomentum(P1, P1_M);
  check *= CheckMomentum(P2, P2_M);
  check *= CheckMomentum(P3, P3_M);
  if (check==0) {
    if (verbosity>=2) cout << "Mass constraint is not respected"<<endl;
    return 0;
  }

  double weight_PS_Init = 2./(TMath::Power(2*TMath::Pi(), 5))/(comEnergy*comEnergy);
  double weight_PS_Higgs = 1/(2*P3.E());
  //double weight_PS_Top;
  double weight_PS_Antitop = P2.Beta()*P2.Vect().Mag()*sin(P2.Theta())/2;
  double weight_PS = weight_PS_Init * weight_PS_Higgs * weight_PS_Top * weight_PS_Antitop;

  return weight_PS; //ME factor
}

double MEPhaseSpace::SetupKinematics_TopHadDecay_WithTopPhaseSpace(const double* x, TLorentzVector* Ptop, int TopType, int DecayType) const{

  //Top decay phase space
  //x[0] Energy Bjet
  //x[1] Theta Bjet
  //x[2] Phi Bjet
  //x[3] Energy Jet1
  //x[4] Theta Jet1
  //x[5] Phi Jet1
  //x[6] Theta Jet2
  //x[7] Phi Jet2 

  if (verbosity>=2) cout << "Top/Antitop momentum reconstruction"<<endl;


  TLorentzVector Bjet;
  TLorentzVector Jet1;
  TLorentzVector Jet2;

  SetMomentumFromEThetaPhi(mB, x[0], x[1], x[2], &Bjet);
  SetMomentumFromEThetaPhi(0, x[3], x[4], x[5], &Jet1);

  //cosAlpha
  TVector3 Jet2unit;
  Jet2unit.SetPtThetaPhi(1, x[6], x[7]);
  Jet2unit.SetMag(1);
  double cosAlpha = (((Bjet+Jet1).Vect()).Unit()).Dot(Jet2unit);
  if (verbosity>=2) cout << "cosAlpha="<<cosAlpha << endl;
  //Find Jet2 Energy
  double a = -mTop*mTop + mB*mB + 2*Bjet.Dot(Jet1);
  double b = -2*(Bjet+Jet1).Vect().Mag()*cosAlpha;
  double c = 2*(Bjet.E()+Jet1.E());
  double Jet2_Mag = -a / (b+c);
  if (Jet2_Mag<0) {
    if (verbosity>=2) cout << "Jet2 E<0, quit" << endl;
    return 0;
  }
  SetMomentumFromEThetaPhi(0, Jet2_Mag, x[6], x[7], &Jet2);

  if (DecayType==kTopLepDecay) Computed_mETvect += Jet2;

  (*Ptop) = Bjet + Jet1 + Jet2;

  FillTopDecayPhaseSpacePoint(*Ptop, Bjet, Jet1, Jet2, TopType);
  if (TopType==kTop) ReadPartonMomenta(pTop, 4);
  if (TopType==kAntitop) ReadPartonMomenta(pAntitop, 4);

  int check = 1;
  check =  CheckMomentum(*Ptop, mTop);
  check *= CheckMomentum(Bjet, mB);
  check *= CheckMomentum(Jet1, 0);
  check *= CheckMomentum(Jet2, 0);
  if (check==0) {
    if (verbosity>=2) cout << "Mass constraint is not respected"<<endl;
    return 0;
  }

  // if (verbosity>=2) cout << "Phase-space computation starts"<<endl;
 
  double dfdPJet2 = TMath::Abs(b + c);
  double weight_PS_Bjet = Bjet.Beta()*Bjet.Vect().Mag()*sin(x[1])/2.;
  double weight_PS_Jet1 = Jet1.Beta()*Jet1.Vect().Mag()*sin(x[4])/2.;
  double weight_PS_Pi = 1/(TMath::Power(2*TMath::Pi(), 8));
  double weight_PS_Jet2 = Jet2.Beta()*Jet2.Vect().Mag2()*sin(x[6])/(2*Jet2.E())/dfdPJet2;
  double weightPS = weight_PS_Pi*weight_PS_Bjet*weight_PS_Jet1*weight_PS_Jet2;

  weightPS /= gammaTop;

  return weightPS * 2;
}

double MEPhaseSpace::SetupKinematics_Higgs2l2nuDecay_WithHiggsPhaseSpace(const double* x, TLorentzVector* Phiggs) const{

  //Higgs decay phase space
  //x[0] Energy Lep1
  //x[1] Theta Lep1
  //x[2] Phi Lep1
  //x[3] Energy Lep2
  //x[4] Theta Lep2
  //x[5] Phi Lep2
  //x[6] Energy Neut1
  //x[7] Theta Neut1
  //x[8] Phi Neut1
  //x[9] Theta Neut2
  //x[10] Phi Neut2 

  if (verbosity>=2) cout << "Higgs momentum reconstruction"<<endl;

  TLorentzVector Lep1;
  TLorentzVector Lep2;
  TLorentzVector Neut1;
  TLorentzVector Neut2;

  SetMomentumFromEThetaPhi(0, x[0], x[1], x[2], &Lep1);
  SetMomentumFromEThetaPhi(0, x[3], x[4], x[5], &Lep2);
  SetMomentumFromEThetaPhi(0, x[6], x[7], x[8], &Neut1);

  //cosAlpha
  TVector3 Neut2unit;
  Neut2unit.SetPtThetaPhi(1, x[9], x[10]);
  Neut2unit.SetMag(1);
  double cosAlpha = (((Lep1+Lep2+Neut1).Vect()).Unit()).Dot(Neut2unit);
  if (verbosity>=2) cout << "cosAlpha="<<cosAlpha << endl;
  if (verbosity>=2) cout << "x[9]="<<x[9]<<" x[10]="<<x[10]<<endl;

  double a = -mHiggs*mHiggs + (Lep1+Lep2+Neut1).Mag2();
  double b = -2*(Lep1+Lep2+Neut1).Vect().Mag()*cosAlpha;
  double c = 2*(Lep1.E()+Lep2.E()+Neut1.E());
  double Neut2_Mag = -a / (b+c);
  if (Neut2_Mag<0) {
    if (verbosity>=2) cout << "a="<<a<<" b="<<b<<" c="<<c<<endl;
    if (verbosity>=2) cout << "Neut2 E<0, quit" << endl;
    return 0;
  }
  SetMomentumFromEThetaPhi(0, Neut2_Mag, x[9], x[10], &Neut2);

  if (iMode==kMEM_TTH_TopAntitopHiggsDecay) Computed_mETvect += (Neut1+Neut2); 
  if (iMode==kMEM_TTH_TopAntitopHiggsSemiLepDecay) Computed_mETvect += Neut2;

  (*Phiggs) = Lep1 + Lep2 + Neut1 + Neut2;

  FillHiggsDecayPhaseSpacePoint(*Phiggs, Lep1, Neut1, Lep2, Neut2);
  ReadPartonMomenta(pHiggs, 5);

  int check = 1;
  check =  CheckMomentum(*Phiggs, mHiggs);
  check *= CheckMomentum(Lep1, 0);
  check *= CheckMomentum(Neut1, 0);
  check *= CheckMomentum(Lep2, 0);
  check *= CheckMomentum(Neut2, 0);
  if (check==0) {
    if (verbosity>=2) cout << "Mass constraint is not respected"<<endl;
    return 0;
  }

  double dfdPNeut2 = TMath::Abs(b + c);
  double weight_PS_Lep1 = Lep1.Beta()*Lep1.Vect().Mag()*sin(x[1])/2.;
  double weight_PS_Lep2 = Lep2.Beta()*Lep2.Vect().Mag()*sin(x[4])/2.;
  double weight_PS_Neut1 = Neut1.Beta()*Neut1.Vect().Mag()*sin(x[7])/2.;
  double weight_PS_Pi = 1/(TMath::Power(2*TMath::Pi(), 11));
  double weight_PS_Neut2 = Neut2.Beta()*Neut2.Vect().Mag2()*sin(x[9])/(2*Neut2.E())/dfdPNeut2;
  double weightPS = weight_PS_Pi*weight_PS_Lep1*weight_PS_Lep2*weight_PS_Neut1*weight_PS_Neut2;

  weightPS /= gammaHiggs;


  return weightPS * 2;
}

double MEPhaseSpace::SetupKinematicsTTH_NoBjorken_TopAntitopHiggsDecay(const double* x) const{


  if (iMode==kMEM_TTH_TopAntitopHiggsDecay && MEMFix_TopHad.TopSign==-1) isTopAntitop = 0;
  else if (iMode==kMEM_TTH_TopAntitopHiggsDecay && MEMFix_TopHad.TopSign==1) isTopAntitop = 1;
  else if (iMode==kMEM_TTH_TopAntitopHiggsSemiLepDecay && MEMFix_TopLep.TopSign==-1) isTopAntitop = 0;
  else if (iMode==kMEM_TTH_TopAntitopHiggsSemiLepDecay && MEMFix_TopLep.TopSign==1) isTopAntitop = 1;

  int TopDecay1=0, TopDecay2=0; 
  if (iMode==kMEM_TTH_TopAntitopHiggsDecay) {
    TopDecay1 = kTopHadDecay;
    TopDecay2 = kTopLepDecay;
  }
  if (iMode==kMEM_TTH_TopAntitopHiggsSemiLepDecay) {
    TopDecay1 = kTopLepDecay;
    TopDecay2 = kTopLepDecay;
  }


  //Top had decay phase space
  //x[0] Energy Bjet
  //x[1] Theta Bjet
  //x[2] Phi Bjet
  //x[3] Energy Jet1
  //x[4] Theta Jet1
  //x[5] Phi Jet1
  //x[6] Theta Jet2
  //x[7] Phi Jet2 

  y[0] = x[0];
  y[1] = x[1];
  y[2] = x[2];
  y[3] = x[3];
  y[4] = x[4];
  y[5] = x[5];
  y[6] = x[6];
  y[7] = x[7];

  TLorentzVector P1;
  double weight_PS_Top1 = 0;
  if (isTopAntitop==0) weight_PS_Top1 = SetupKinematics_TopHadDecay_WithTopPhaseSpace(y, &P1, kTop, TopDecay1);
  else if (isTopAntitop==1) weight_PS_Top1 = SetupKinematics_TopHadDecay_WithTopPhaseSpace(y, &P1, kAntitop, TopDecay1);
  if (weight_PS_Top1==0) return 0;


  //Top lep decay phase space
  //x[8] Energy Bjet
  //x[9] Theta Bjet
  //x[10] Phi Bjet
  //x[11] Energy Jet1
  //x[12] Theta Jet1
  //x[13] Phi Jet1
  //x[14] Theta Jet2
  //x[15] Phi Jet2 
  
  y[0] = x[8];
  y[1] = x[9];
  y[2] = x[10];
  y[3] = x[11];
  y[4] = x[12];
  y[5] = x[13];
  y[6] = x[14];
  y[7] = x[15];

  TLorentzVector P2;
  double weight_PS_Top2 = 0;
  if (isTopAntitop==0) weight_PS_Top2 = SetupKinematics_TopHadDecay_WithTopPhaseSpace(y, &P2, kAntitop, TopDecay2);
  else if (isTopAntitop==1) weight_PS_Top2 = SetupKinematics_TopHadDecay_WithTopPhaseSpace(y, &P2, kTop, TopDecay2);

  //Higgs decay phase space
  //x[16] Energy Lep1
  //x[17] Theta Lep1
  //x[18] Phi Lep1
  //x[19] Energy Lep2
  //x[20] Theta Lep2
  //x[21] Phi Lep2
  //x[22] Energy Neut1
  //x[23] Theta Neut1
  //x[24] Phi Neut1
  //x[25] Theta Neut2
  //x[26] Phi Neut2 

  y[0] = x[16];
  y[1] = x[17];
  y[2] = x[18];
  y[3] = x[19];
  y[4] = x[20];
  y[5] = x[21];
  y[6] = x[22];
  y[7] = x[23];
  y[8] = x[24];
  y[9] = x[25];
  y[10] = x[26];

  TLorentzVector P3;
  double weight_PS_Higgs = SetupKinematics_Higgs2l2nuDecay_WithHiggsPhaseSpace(y, &P3);
  if (weight_PS_Higgs==0) return 0;

  TLorentzVector Ptot = P1 + P2 + P3;
  double x1 = (Ptot.Pz()+Ptot.E())/comEnergy;
  double x2 = (-Ptot.Pz()+Ptot.E())/comEnergy;

  if (verbosity>=2) cout << "x1+x2="<<x1+x2<<endl;
  if (!(x2>0 && x2<1)){
    if (verbosity>=2) cout << "x2="<<x2<<" out of range, quit"<<endl;
    return 0;
  }
  if (!(x1>0 && x1<1)){
    if (verbosity>=2) cout << "x1="<<x1<<" out of range, quit"<<endl;
    return 0;
  }

  if (verbosity>=2) cout << "Ptot Pt="<<Ptot.Pt()<<endl;

    if (isTopAntitop==0) FillTTHPhaseSpacePoint(P1, P2, P3);
    else FillTTHPhaseSpacePoint(P2, P1, P3);
  SetInitialPartonMomenta(x1, x2);

  ReadPartonMomenta(pCore, 5);
  int check = 1;
  check =  CheckMomentum(P1, mTop);
  check *= CheckMomentum(P2, mTop);
  check *= CheckMomentum(P3, mHiggs);
  if (check==0) {
    if (verbosity>=2) cout << "Mass constraint is not respected"<<endl;
    return 0;
  }

  double weight_PS_Init = 2./(comEnergy*comEnergy);
  double weight_PS = weight_PS_Init * weight_PS_Top1 * weight_PS_Top2 * weight_PS_Higgs;

  return weight_PS;
}

double MEPhaseSpace::SetupKinematicsTTW_NoBjorken_TopAntitopDecay(const double* x) const{

  if (MEMFix_TopLep.TopSign==-1) isTopAntitop = 0;
  else if (MEMFix_TopLep.TopSign==1) isTopAntitop = 1;

  int TopDecay1 = kTopLepDecay;
  int TopDecay2 = kTopLepDecay;

  y[0] = x[0];
  y[1] = x[1];
  y[2] = x[2];
  y[3] = x[3];
  y[4] = x[4];
  y[5] = x[5];
  y[6] = x[6];
  y[7] = x[7];

  TLorentzVector P1;
  double weight_PS_Top1 = 0;
  if (isTopAntitop==0) weight_PS_Top1  = SetupKinematics_TopHadDecay_WithTopPhaseSpace(y, &P1, kTop, TopDecay1);
  if (isTopAntitop==1) weight_PS_Top1  = SetupKinematics_TopHadDecay_WithTopPhaseSpace(y, &P1, kAntitop, TopDecay1);
  if (weight_PS_Top1==0) return 0;

  y[0] = x[8];
  y[1] = x[9];
  y[2] = x[10];
  y[3] = x[11];
  y[4] = x[12];
  y[5] = x[13];
  y[6] = x[14];
  y[7] = x[15];

  TLorentzVector P2;
  double weight_PS_Top2 = 0;
  if (isTopAntitop==0) weight_PS_Top2 = SetupKinematics_TopHadDecay_WithTopPhaseSpace(y, &P2, kAntitop, TopDecay2);
  if (isTopAntitop==1) weight_PS_Top2 = SetupKinematics_TopHadDecay_WithTopPhaseSpace(y, &P2, kTop, TopDecay2);
  if (weight_PS_Top2==0) return 0;

  TLorentzVector P3;
  TLorentzVector P4;
  SetMomentumFromEThetaPhi(0, x[16], x[17], x[18], &P3);  //lepton
  SetMomentumFromEThetaPhi(0, x[19], x[20], x[21], &P4);  //neutrino

  Computed_mETvect += P4;

/*
  TVector3 P3u = P3.Vect().Unit();

  //mW
  double tW = x[19];  // belongs to [-arctan(mW/gammaW) ; arctan(mW/gammaW)]
  double W_M = sqrt(mW*mW + mW*gammaW*tan(tW)); 
  if (verbosity>=2) cout << "tW="<<tW<<" W_M="<<W_M<<endl; 
 
  //Neutrino theta ?
  double Neut_E = x[20];
  double Neut_Phi = x[21];
  double cosAlpha = 1 - W_M*W_M/(2*Neut_E*P3.E());
  double a = (P3u.X()*cos(Neut_Phi)+P3u.Y()*sin(Neut_Phi))*(P3u.X()*cos(Neut_Phi)+P3u.Y()*sin(Neut_Phi)) + P3u.Z()*P3u.Z();
  double b = -2*(P3u.X()*cos(Neut_Phi)+P3u.Y()*sin(Neut_Phi)) * cosAlpha;
  double c = cosAlpha*cosAlpha - P3u.Z()*P3u.Z();
  double delta = b*b-4*a*c;
  if (verbosity>=2) cout << "Neutrino delta="<<delta<<endl;
  if (delta<0) return 0;
  double X1 = (-b+sqrt(delta))/(2*a);
  double X2 = (-b-sqrt(delta))/(2*a);
  if (verbosity>=2) cout << "Neutrino sinTheta X1="<<X1<<" X2="<<X2<<endl;
  double Neut_Theta;
  int foundsol = 0;
  if (X1>-1 && X1<1){
    if (X1>0) Neut_Theta = TMath::ASin(X1); 
    if (X1<0) Neut_Theta = TMath::Abs(TMath::ASin(X1)) + TMath::Pi(); 
   SetMomentumFromEThetaPhi(0, Neut_E, Neut_Theta, Neut_Phi, &P4);
    double W_M_check = (P3+P4).M(); 
    if (verbosity>=2) cout << "Check Wmass "<<W_M<<" redone: "<<W_M_check<<" (theta="<<Neut_Theta <<")"<<endl;
    int check = CheckMomentum(P4, 0); 
    if (check && W_M==W_M_check) foundsol++;
  }
  if (X2>-1 && X2<1){
    if (X2>0) Neut_Theta = TMath::ASin(X2);
    if (X2<0) Neut_Theta = -TMath::Abs(TMath::ASin(X2))  + TMath::Pi();
    SetMomentumFromEThetaPhi(0, Neut_E, Neut_Theta, Neut_Phi, &P4);  
    double W_M_check = (P3+P4).M();   
    if (verbosity>=2) cout << "Check Wmass "<<W_M<<" redone: "<<W_M_check<<" (theta="<<Neut_Theta <<")"<<endl;
    int check = CheckMomentum(P4, 0);
    if (check && W_M==W_M_check) foundsol++;
  }
  if (foundsol==0) return 0;
  if (foundsol==2) {
    if (verbosity>=2) cout << "Found two viable neutrino Theta" << endl;
  }
*/

  TLorentzVector Ptot = P1 + P2 + P3 + P4;
  double x1 = (Ptot.Pz()+Ptot.E())/comEnergy;
  double x2 = (-Ptot.Pz()+Ptot.E())/comEnergy;

  if (verbosity>=2) cout << "x1+x2="<<x1+x2<<endl;
  if (!(x2>0 && x2<1)){
    if (verbosity>=2) cout << "x2="<<x2<<" out of range, quit"<<endl;
    return 0;
  }
  if (!(x1>0 && x1<1)){
    if (verbosity>=2) cout << "x1="<<x1<<" out of range, quit"<<endl;
    return 0;
  }
  if (verbosity>=2) cout << "Ptot Pt="<<Ptot.Pt()<<endl;


  if (isTopAntitop==0) FillTTLLPhaseSpacePoint(P1, P2, P3, P4);
  else FillTTLLPhaseSpacePoint(P2, P1, P3, P4);
  SetInitialPartonMomenta(x1, x2);

  ReadPartonMomenta(pCore, 6);
  int check = 1;
  check =  CheckMomentum(P1, mTop);
  check *= CheckMomentum(P2, mTop);
  check *= CheckMomentum(P3, 0);
  check *= CheckMomentum(P4, 0);
  if (check==0) {
    if (verbosity>=2) cout << "Mass constraint is not respected"<<endl;
    return 0;
  }

  double weight_PS_Init = 2./(comEnergy*comEnergy);
  double weight_PS_W = P3.Beta()*P3.Vect().Mag()*sin(x[17])/2. * P4.Beta()*P4.Vect().Mag()*sin(x[20])/2. /(TMath::Power(2*TMath::Pi(),6));
  double weight_PS = weight_PS_Init * weight_PS_Top1 * weight_PS_Top2 * weight_PS_W;

  return weight_PS;
}


double MEPhaseSpace::SetupKinematicsTTWJJ_NoBjorken_TopAntitopDecay(const double* x) const{

  if (MEMFix_TopLep.TopSign==-1) isTopAntitop = 0;
  else if (MEMFix_TopLep.TopSign==1) isTopAntitop = 1;

  int TopDecay1 = kTopLepDecay;
  int TopDecay2 = kTopLepDecay;

  y[0] = x[0];
  y[1] = x[1];
  y[2] = x[2];
  y[3] = x[3];
  y[4] = x[4];
  y[5] = x[5];
  y[6] = x[6];
  y[7] = x[7];

  TLorentzVector P1;
  double weight_PS_Top1 = 0;
  if (isTopAntitop==0) weight_PS_Top1  = SetupKinematics_TopHadDecay_WithTopPhaseSpace(y, &P1, kTop, TopDecay1);
  if (isTopAntitop==1) weight_PS_Top1  = SetupKinematics_TopHadDecay_WithTopPhaseSpace(y, &P1, kAntitop, TopDecay1);
  if (weight_PS_Top1==0) return 0;

  y[0] = x[8];
  y[1] = x[9];
  y[2] = x[10];
  y[3] = x[11];
  y[4] = x[12];
  y[5] = x[13];
  y[6] = x[14];
  y[7] = x[15];

  TLorentzVector P2;
  double weight_PS_Top2 = 0;
  if (isTopAntitop==0) weight_PS_Top2 = SetupKinematics_TopHadDecay_WithTopPhaseSpace(y, &P2, kAntitop, TopDecay2);
  if (isTopAntitop==1) weight_PS_Top2 = SetupKinematics_TopHadDecay_WithTopPhaseSpace(y, &P2, kTop, TopDecay2);
  if (weight_PS_Top2==0) return 0;

  TLorentzVector P3;
  TLorentzVector P4;
  SetMomentumFromEThetaPhi(0, x[16], x[17], x[18], &P3);  //lepton
  SetMomentumFromEThetaPhi(0, x[19], x[20], x[21], &P4);  //neutrino

  Computed_mETvect += P4;

  TLorentzVector P5;
  TLorentzVector P6;
  SetMomentumFromEThetaPhi(0, x[22], x[23], x[24], &P5);  //jet1
  SetMomentumFromEThetaPhi(0, x[25], x[26], x[27], &P6);  //jet2

  TLorentzVector Ptot = P1 + P2 + P3 + P4 + P5 + P6;
  double x1 = (Ptot.Pz()+Ptot.E())/comEnergy;
  double x2 = (-Ptot.Pz()+Ptot.E())/comEnergy;

  if (verbosity>=2) cout << "x1+x2="<<x1+x2<<endl;
  if (!(x2>0 && x2<1)){
    if (verbosity>=2) cout << "x2="<<x2<<" out of range, quit"<<endl;
    return 0;
  }
  if (!(x1>0 && x1<1)){
    if (verbosity>=2) cout << "x1="<<x1<<" out of range, quit"<<endl;
    return 0;
  }
  if (verbosity>=2) cout << "Ptot Pt="<<Ptot.Pt()<<endl;


  if (isTopAntitop==0) FillTTLNuJJPhaseSpacePoint(P1, P2, P3, P4, P5, P6);
  else FillTTLNuJJPhaseSpacePoint(P2, P1, P3, P4, P5, P6);
  SetInitialPartonMomenta(x1, x2);

  ReadPartonMomenta(pCore, 8);
  int check = 1;
  check =  CheckMomentum(P1, mTop);
  check *= CheckMomentum(P2, mTop);
  check *= CheckMomentum(P3, 0);
  check *= CheckMomentum(P4, 0);
  check *= CheckMomentum(P5, 0);
  check *= CheckMomentum(P6, 0);
  if (check==0) {
    if (verbosity>=2) cout << "Mass constraint is not respected"<<endl;
    return 0;
  }

  double weight_PS_Init = 2./(comEnergy*comEnergy);
  double weight_PS_W = P3.Beta()*P3.Vect().Mag()*sin(x[17])/2. * P4.Beta()*P4.Vect().Mag()*sin(x[20])/2. /(TMath::Power(2*TMath::Pi(),6));
  double weight_PS_JJ = P5.Beta()*P5.Vect().Mag()*sin(x[23])/2. * P6.Beta()*P6.Vect().Mag()*sin(x[26])/2. /(TMath::Power(2*TMath::Pi(),6));
  double weight_PS = weight_PS_Init * weight_PS_Top1 * weight_PS_Top2 * weight_PS_W * weight_PS_JJ;

  return weight_PS;

}


double MEPhaseSpace::SetupKinematicsTTLL_NoBjorken_TopAntitopDecay(const double* x) const{

  if (MEMFix_TopHad.TopSign==-1) isTopAntitop = 0;
  else if (MEMFix_TopHad.TopSign==1) isTopAntitop = 1;

  int TopDecay1 = kTopHadDecay;
  int TopDecay2 = kTopLepDecay;

  y[0] = x[0];
  y[1] = x[1];
  y[2] = x[2];
  y[3] = x[3];
  y[4] = x[4];
  y[5] = x[5];
  y[6] = x[6];
  y[7] = x[7];

  TLorentzVector P1;
  double weight_PS_Top1 = 0;
  if (isTopAntitop==0) weight_PS_Top1  = SetupKinematics_TopHadDecay_WithTopPhaseSpace(y, &P1, kTop, TopDecay1);
  if (isTopAntitop==1) weight_PS_Top1  = SetupKinematics_TopHadDecay_WithTopPhaseSpace(y, &P1, kAntitop, TopDecay1);
  if (weight_PS_Top1==0) return 0;

  y[0] = x[8];
  y[1] = x[9];
  y[2] = x[10];
  y[3] = x[11];
  y[4] = x[12];
  y[5] = x[13];
  y[6] = x[14];
  y[7] = x[15];

  TLorentzVector P2;
  double weight_PS_Top2 = 0;
  if (isTopAntitop==0) weight_PS_Top2 = SetupKinematics_TopHadDecay_WithTopPhaseSpace(y, &P2, kAntitop, TopDecay2);
  if (isTopAntitop==1) weight_PS_Top2 = SetupKinematics_TopHadDecay_WithTopPhaseSpace(y, &P2, kTop, TopDecay2);
  if (weight_PS_Top2==0) return 0;

  TLorentzVector P3;
  TLorentzVector P4;
  SetMomentumFromEThetaPhi(0, x[16], x[17], x[18], &P3);
  SetMomentumFromEThetaPhi(0, x[19], x[20], x[21], &P4);

  TLorentzVector Ptot = P1 + P2 + P3 + P4;
  double x1 = (Ptot.Pz()+Ptot.E())/comEnergy;
  double x2 = (-Ptot.Pz()+Ptot.E())/comEnergy;

  if (verbosity>=2) cout << "x1+x2="<<x1+x2<<endl;
  if (!(x2>0 && x2<1)){
    if (verbosity>=2) cout << "x2="<<x2<<" out of range, quit"<<endl;
    return 0;
  }
  if (!(x1>0 && x1<1)){
    if (verbosity>=2) cout << "x1="<<x1<<" out of range, quit"<<endl;
    return 0;
  }

  if (verbosity>=2) cout << "Ptot Pt="<<Ptot.Pt()<<endl;

  if (isTopAntitop==0) FillTTLLPhaseSpacePoint(P1, P2, P3, P4);
  else FillTTLLPhaseSpacePoint(P2, P1, P3, P4);
  SetInitialPartonMomenta(x1, x2);

  ReadPartonMomenta(pCore, 6);
  int check = 1;
  check =  CheckMomentum(P1, mTop);
  check *= CheckMomentum(P2, mTop);
  check *= CheckMomentum(P3, 0);
  check *= CheckMomentum(P4, 0);
  if (check==0) {
    if (verbosity>=2) cout << "Mass constraint is not respected"<<endl;
    return 0;
  }

  double weight_PS_Init = 2./(comEnergy*comEnergy);
  double weight_PS_LL = P3.Beta()*P3.Vect().Mag()*sin(x[17])/2. * P4.Beta()*P4.Vect().Mag()*sin(x[20])/2. /(TMath::Power(2*TMath::Pi(),6));
  double weight_PS = weight_PS_Init * weight_PS_Top1 * weight_PS_Top2 * weight_PS_LL;

  return weight_PS;
}

double MEPhaseSpace::ComputeDecayMomenta(TLorentzVector& P, double M1, double M2, double theta, double phi, TLorentzVector* Decay1, TLorentzVector* Decay2) const{

  //Energies of decay particles in cm
  TLorentzVector Pcm(0,0,0,P.M());
  double Kallen = KallenFunction(Pcm.E(), M1, M2);
  if (Kallen<0){
   if (verbosity>=2) cout <<"Kallen("<<P.M()<<","<<M1<<","<<M2<<")="<< Kallen<<"<0, Kinematics forbidden, quit"<<endl;
    return 0;
  }
  double Decay_P = sqrt(Kallen) / (2*P.M());
  double Decay1_E = sqrt(M1*M1+Decay_P*Decay_P);
  //double Decay2_E = sqrt(M2*M2+Decay_P*Decay_P);

  //Apply the opening angles
    double Decay1_Phi = phi;
    double Decay1_Theta = theta;
    double Decay1_Eta = -log(tan(Decay1_Theta/2));
    if (Decay1_E<M1) {
      if (verbosity>=2) cout << "Decay1 E<M, quit"<<endl;
      return 0;
    }
    double Decay1_Pt = sqrt(Decay1_E*Decay1_E-M1*M1) / TMath::CosH(Decay1_Eta);
    (*Decay1).SetPtEtaPhiE(Decay1_Pt, Decay1_Eta, Decay1_Phi, Decay1_E);
    (*Decay2) = Pcm - (*Decay1);

  //Boost
  (*Decay1).Boost(P.BoostVector());
  (*Decay2).Boost(P.BoostVector());

  //Return phase space
  double weightPS;
   //if (P.Vect().Mag()==0)  
     weightPS = Decay_P * sin(theta) / (4*TMath::Power(2*TMath::Pi(), 2)) / P.M();
   //else 
     //weightPS = (*Decay1).Vect().Mag() * sin((*Decay1).Theta()) / (4*TMath::Power(2*TMath::Pi(), 2)) / P.E();

  return weightPS;
}

int MEPhaseSpace::SetMomentumFromEThetaPhi(double M, double E, double Theta, double Phi, TLorentzVector* P2) const{

  double P2_M = M;
  double P2_Theta = Theta;
  double P2_Phi = Phi;
  double P2_E = E;
  if (P2_E<P2_M){
    if (verbosity>=2) cout << "Energy out of range, quit"<<endl;
    return 0;
  }
  double P2_Eta = -log(tan(P2_Theta/2.));
  double P2_Pt = sqrt(-P2_M*P2_M+P2_E*P2_E) / TMath::CosH(P2_Eta);
  (*P2).SetPtEtaPhiE(P2_Pt, P2_Eta, P2_Phi, P2_E);

  return 1;
}

double MEPhaseSpace::ComputeSubMatrixElement(int iProc, int ip1, int ip2) const {

  double weight = 0;
  const double* matrix_elements = 0;
    if (iGen==kMadgraph){
      if (iProc==kTTH){
        process->setMomenta(*pCore);
        process->sigmaKin();
        matrix_elements = process->getMatrixElements();
        weight = matrix_elements[0];
        if (verbosity>=2) cout << "TTH ME = "<<weight<<endl;
      }
      if (iProc==kTTLL){
        process_ggttll->setMomenta(*pCore);
        process_ggttll->sigmaKin();
        matrix_elements = process_ggttll->getMatrixElements();
        weight = matrix_elements[0];
        if (verbosity>=2) cout << "TTLL ME = "<<weight<<endl;
      }
      if (iProc==kTTW){
        if (MEMFix_HiggsSemiLep.LepSign>0){
	  if ((ip1==-1 && ip2==4) || (ip1==4 && ip2==-1)){
            process_qqttlpvl_cdx->setInitial(ip1, ip2);
            process_qqttlpvl_cdx->setMomenta(*pCore);
            process_qqttlpvl_cdx->sigmaKin();
            weight = process_qqttlpvl_cdx->sigmaHat();
	  }
	  if ((ip1==-3 && ip2==4) || (ip1==4 && ip2==-3) || (ip1==2 && ip2==-1) || (ip1==-1 && ip2==2)){
            process_qqttlpvl_udx->setInitial(ip1, ip2);
            process_qqttlpvl_udx->setMomenta(*pCore);
            process_qqttlpvl_udx->sigmaKin();
            weight = process_qqttlpvl_udx->sigmaHat();
	  }
	  if ((ip1==-3 && ip2==2) || (ip1==2 && ip2==-3)){
            process_qqttlpvl_usx->setInitial(ip1, ip2);
            process_qqttlpvl_usx->setMomenta(*pCore);
            process_qqttlpvl_usx->sigmaKin();
            weight = process_qqttlpvl_usx->sigmaHat();
	  }
        }
        else if (MEMFix_HiggsSemiLep.LepSign<0){
	  if ((ip1==-4 && ip2==1) || (ip1==1 && ip2==-4)){
            process_qqttlmvl_dcx->setInitial(ip1, ip2);
            process_qqttlmvl_dcx->setMomenta(*pCore);
            process_qqttlmvl_dcx->sigmaKin();
            weight = process_qqttlmvl_dcx->sigmaHat();
	  }
          if ((ip1==-4 && ip2==3) || (ip1==3 && ip2==-4) || (ip1==1 && ip2==-2) || (ip1==-2 && ip2==1)){
            process_qqttlmvl_dux->setInitial(ip1, ip2);
            process_qqttlmvl_dux->setMomenta(*pCore);
            process_qqttlmvl_dux->sigmaKin();
            weight = process_qqttlmvl_dux->sigmaHat();
          }
	  if ((ip1==-2 && ip2==3) || (ip1==3 && ip2==-2)){
            process_qqttlmvl_sux->setInitial(ip1, ip2);
            process_qqttlmvl_sux->setMomenta(*pCore);
            process_qqttlmvl_sux->sigmaKin();
            weight = process_qqttlmvl_sux->sigmaHat();
          }
        }
        if (verbosity>=2) cout << "TTW ip1="<<ip1<<" ip2="<<ip2<<" ME = "<<weight<<endl;
      }
      if (iProc==kTTWJJ){
        if (MEMFix_HiggsSemiLep.LepSign<0){
          if ((ip1==21 && ip2==-4) || (ip1==-4 && ip2==21)){
	    //process_P0_Sigma_sm_ckm_gcx_ttxemvexgdx->setInitial(ip1, ip2);
            //process_P0_Sigma_sm_ckm_gcx_ttxemvexgdx->setMomenta(*pCore);
            //process_P0_Sigma_sm_ckm_gcx_ttxemvexgdx->sigmaKin();
            //weight = process_P0_Sigma_sm_ckm_gcx_ttxemvexgdx->sigmaHat();
	    process_P0_Sigma_sm_ckm_gux_ttxemvexgdx->setInitial(ip1, ip2);
            process_P0_Sigma_sm_ckm_gux_ttxemvexgdx->setMomenta(*pCore);
            process_P0_Sigma_sm_ckm_gux_ttxemvexgdx->sigmaKin();
            weight = process_P0_Sigma_sm_ckm_gux_ttxemvexgdx->sigmaHat();
          }
          if ((ip1==21 && ip2==1) || (ip1==1 && ip2==21)){
	    process_P0_Sigma_sm_ckm_gd_ttxemvexgc->setInitial(ip1, ip2);
            process_P0_Sigma_sm_ckm_gd_ttxemvexgc->setMomenta(*pCore);
            process_P0_Sigma_sm_ckm_gd_ttxemvexgc->sigmaKin();
            weight = process_P0_Sigma_sm_ckm_gd_ttxemvexgc->sigmaHat();
	    process_P0_Sigma_sm_ckm_gd_ttxemvexgu->setInitial(ip1, ip2);
            process_P0_Sigma_sm_ckm_gd_ttxemvexgu->setMomenta(*pCore);
            process_P0_Sigma_sm_ckm_gd_ttxemvexgu->sigmaKin();
            weight += process_P0_Sigma_sm_ckm_gd_ttxemvexgu->sigmaHat();
          }
          if ((ip1==21 && ip2==3) || (ip1==3 && ip2==21)){
	    process_P0_Sigma_sm_ckm_gd_ttxemvexgu->setInitial(ip1, ip2);
            process_P0_Sigma_sm_ckm_gd_ttxemvexgu->setMomenta(*pCore);
            process_P0_Sigma_sm_ckm_gd_ttxemvexgu->sigmaKin();
            weight = process_P0_Sigma_sm_ckm_gd_ttxemvexgu->sigmaHat();
	    process_P0_Sigma_sm_ckm_gs_ttxemvexgu->setInitial(ip1, ip2);
            process_P0_Sigma_sm_ckm_gs_ttxemvexgu->setMomenta(*pCore);
            process_P0_Sigma_sm_ckm_gs_ttxemvexgu->sigmaKin();
            weight += process_P0_Sigma_sm_ckm_gs_ttxemvexgu->sigmaHat();
          }
          if ((ip1==21 && ip2==-2) || (ip1==-2 && ip2==21)){
            process_P0_Sigma_sm_ckm_gux_ttxemvexgdx->setInitial(ip1, ip2);
            process_P0_Sigma_sm_ckm_gux_ttxemvexgdx->setMomenta(*pCore);
            process_P0_Sigma_sm_ckm_gux_ttxemvexgdx->sigmaKin();
            weight = process_P0_Sigma_sm_ckm_gux_ttxemvexgdx->sigmaHat();
	    process_P0_Sigma_sm_ckm_gux_ttxemvexgsx->setInitial(ip1, ip2);
            process_P0_Sigma_sm_ckm_gux_ttxemvexgsx->setMomenta(*pCore);
            process_P0_Sigma_sm_ckm_gux_ttxemvexgsx->sigmaKin();
            weight += process_P0_Sigma_sm_ckm_gux_ttxemvexgsx->sigmaHat();
          }
        }
	if (MEMFix_HiggsSemiLep.LepSign>0){
	  if ((ip1==21 && ip2==4) || (ip1==4 && ip2==21)){
	    process_P0_Sigma_sm_ckm_gc_ttxepvegd->setInitial(ip1, ip2);
            process_P0_Sigma_sm_ckm_gc_ttxepvegd->setMomenta(*pCore);
            process_P0_Sigma_sm_ckm_gc_ttxepvegd->sigmaKin();
	    weight = process_P0_Sigma_sm_ckm_gc_ttxepvegd->sigmaHat();
	    process_P0_Sigma_sm_ckm_gu_ttxepvegd->setInitial(ip1, ip2);
            process_P0_Sigma_sm_ckm_gu_ttxepvegd->setMomenta(*pCore);
            process_P0_Sigma_sm_ckm_gu_ttxepvegd->sigmaKin();
            weight += process_P0_Sigma_sm_ckm_gu_ttxepvegd->sigmaHat();
	  }
          if ((ip1==21 && ip2==-1) || (ip1==-1 && ip2==21)){
	    process_P0_Sigma_sm_ckm_gdx_ttxepvegcx->setInitial(ip1, ip2);
            process_P0_Sigma_sm_ckm_gdx_ttxepvegcx->setMomenta(*pCore);
            process_P0_Sigma_sm_ckm_gdx_ttxepvegcx->sigmaKin();
            weight = process_P0_Sigma_sm_ckm_gdx_ttxepvegcx->sigmaHat();
            process_P0_Sigma_sm_ckm_gdx_ttxepvegux->setInitial(ip1, ip2);
            process_P0_Sigma_sm_ckm_gdx_ttxepvegux->setMomenta(*pCore);
            process_P0_Sigma_sm_ckm_gdx_ttxepvegux->sigmaKin();
            weight += process_P0_Sigma_sm_ckm_gdx_ttxepvegux->sigmaHat();
	  }
          if ((ip1==21 && ip2==-3) || (ip1==-3 && ip2==21)){
            process_P0_Sigma_sm_ckm_gdx_ttxepvegux->setInitial(ip1, ip2);
            process_P0_Sigma_sm_ckm_gdx_ttxepvegux->setMomenta(*pCore);
            process_P0_Sigma_sm_ckm_gdx_ttxepvegux->sigmaKin();
            weight = process_P0_Sigma_sm_ckm_gdx_ttxepvegux->sigmaHat();
	    process_P0_Sigma_sm_ckm_gsx_ttxepvegux->setInitial(ip1, ip2);
            process_P0_Sigma_sm_ckm_gsx_ttxepvegux->setMomenta(*pCore);
            process_P0_Sigma_sm_ckm_gsx_ttxepvegux->sigmaKin();
            weight += process_P0_Sigma_sm_ckm_gsx_ttxepvegux->sigmaHat();
          }
          if ((ip1==21 && ip2==2) || (ip1==2 && ip2==21)){
            process_P0_Sigma_sm_ckm_gu_ttxepvegd->setInitial(ip1, ip2);
            process_P0_Sigma_sm_ckm_gu_ttxepvegd->setMomenta(*pCore);
            process_P0_Sigma_sm_ckm_gu_ttxepvegd->sigmaKin();
            weight = process_P0_Sigma_sm_ckm_gu_ttxepvegd->sigmaHat();
	    //process_P0_Sigma_sm_ckm_gu_ttxepvegs->setInitial(ip1, ip2);
            //process_P0_Sigma_sm_ckm_gu_ttxepvegs->setMomenta(*pCore);
            //process_P0_Sigma_sm_ckm_gu_ttxepvegs->sigmaKin();
            //weight += process_P0_Sigma_sm_ckm_gu_ttxepvegs->sigmaHat();
          }
	}
      if (verbosity>=2) cout << "TTWJJ ip1="<<ip1<<" ip2="<<ip2<<" ME = "<<weight<<endl;
      }

      if (iProc==kTop){
        process_tbwjj->setMomenta(*pTop);
        process_tbwjj->sigmaKin();
        matrix_elements = process_tbwjj->getMatrixElements();
        weight = matrix_elements[0];
	if (verbosity>=2) cout << "Top ME = "<<weight<<endl;
      }
      if (iProc==kAntitop){
        process_antitbwjj->setMomenta(*pAntitop);
        process_antitbwjj->sigmaKin();
        matrix_elements = process_antitbwjj->getMatrixElements();
        weight = matrix_elements[0];
	if (verbosity>=2) cout << "Antitop ME = "<<weight<<endl;
      }
    }

  return weight;
}


double MEPhaseSpace::ComputeMatrixElement() const {

  double weight = 0;
  const double* matrix_elements = 0;
  if (iMode==kMEM_TTH_TopAntitopHiggsDecay || iMode==kMEM_TTLL_TopAntitopDecay || iMode==kMEM_TTH_TopAntitopHiggsSemiLepDecay || iMode==kMEM_TTW_TopAntitopDecay){
    if (iGen==kMadgraph){
      if (iCore==kTTH){
        process->setMomenta(*pCore);
        process->sigmaKin();
        matrix_elements = process->getMatrixElements();
        weight = matrix_elements[0];
      }
      if (iCore==kTTLL){
        process_ggttll->setMomenta(*pCore);
        process_ggttll->sigmaKin();
        matrix_elements = process_ggttll->getMatrixElements();
        weight = matrix_elements[0];
      }
      /*
      if (iCore==kTTW){
	if (MEMFix_HiggsSemiLep.LepSign>0){
          process_qqttlpvl->setMomenta(*pCore);
          process_qqttlpvl->sigmaKin();
          matrix_elements = process_qqttlpvl->getMatrixElements();
          weight = matrix_elements[0];
	}
        else {
          process_qqttlmvl->setMomenta(*pCore);
          process_qqttlmvl->sigmaKin();
          matrix_elements = process_qqttlmvl->getMatrixElements();
          weight = matrix_elements[0];
        }
      } 
      */     
      if (verbosity>=2) cout << "Core process ME = "<<weight<<endl;
      process_tbwjj->setMomenta(*pTop);
      process_tbwjj->sigmaKin();
      matrix_elements = process_tbwjj->getMatrixElements();
      weight *= matrix_elements[0];
      if (verbosity>=2) cout << "Top ME = "<<matrix_elements[0]<<endl;
      process_antitbwjj->setMomenta(*pAntitop);
      process_antitbwjj->sigmaKin();
      matrix_elements = process_antitbwjj->getMatrixElements();
      weight *= matrix_elements[0];
      if (verbosity>=2) cout << "Antitop ME = "<<matrix_elements[0]<<endl;
    }
  }
  if (iMode==kAllPartonsTTH_TopDecay){
    if (iGen==kMadgraph){
      process->setMomenta(*pCore);
      process->sigmaKin();
      matrix_elements = process->getMatrixElements();
      weight = matrix_elements[0];

      process_tbwjj->setMomenta(*pTop);
      process_tbwjj->sigmaKin();
      matrix_elements = process_tbwjj->getMatrixElements();
      weight *= matrix_elements[0];
    }
  }
  if (iMode==kAllPartonsTTLL){
      process_ggttll->setMomenta(*pCore);
      process_ggttll->sigmaKin();
      matrix_elements = process_ggttll->getMatrixElements();
      weight = matrix_elements[0];
  }
  if (iMode==kAllPartonsTTH || iMode==kInitialPartons){
    if (iGen==kMadgraph){ 
      process->setMomenta(*pCore);
      process->sigmaKin();
      matrix_elements = process->getMatrixElements();
      weight = matrix_elements[0];
    }
    /*
    else if (iGen==kGosam){
      double* p_momenta1 = new double[25];
      double mmu = 1;
      double param(1.);
      double p_result[4];

      //glu
      p_momenta1[0] = pCore->at(0)[0];
      p_momenta1[1] = pCore->at(0)[1];
      p_momenta1[2] = pCore->at(0)[2];
      p_momenta1[3] = pCore->at(0)[3];
      p_momenta1[4] = 0;

      //glu
      p_momenta1[5] = pCore->at(1)[0];
      p_momenta1[6] = pCore->at(1)[1];
      p_momenta1[7] = pCore->at(1)[2];
      p_momenta1[8] = pCore->at(1)[3];
      p_momenta1[9] = 0;

      //h
      p_momenta1[10] = pCore->at(4)[0];
      p_momenta1[11] = pCore->at(4)[1];
      p_momenta1[12] = pCore->at(4)[2];
      p_momenta1[13] = pCore->at(4)[3];
      p_momenta1[14] = mHiggs;

      //top
      p_momenta1[15] = pCore->at(2)[0];
      p_momenta1[16] = pCore->at(2)[1];
      p_momenta1[17] = pCore->at(2)[2];
      p_momenta1[18] = pCore->at(2)[3];
      p_momenta1[19] = mTop;

      //top
      p_momenta1[20] = pCore->at(3)[0];
      p_momenta1[21] = pCore->at(3)[1];
      p_momenta1[22] = pCore->at(3)[2];
      p_momenta1[23] = pCore->at(3)[3];
      p_momenta1[24] = mTop;

      OLP_EvalSubProcess(2,p_momenta1,mmu,&param,p_result);

      matrix_elements = p_result;
      weight = matrix_elements[3]* (4 * 3.1416 * 0.118 * 4 * 3.1416 * 0.118);
    }
    */
  }
  if (iMode==kAllPartonsTopHad || iMode==kAllPartonsTopHad_FiniteWidth || iMode==kAllPartonsTopHad_FixedTopM_FiniteWwidth){
    process_tbwjj->setMomenta(*pTop);
    process_tbwjj->sigmaKin();
    matrix_elements = process_tbwjj->getMatrixElements();
      weight = matrix_elements[0];
  }
  if (iMode==kAllPartonsAntiTopHad_FixedTopM_FiniteWwidth){
    process_antitbwjj->setMomenta(*pAntitop);
    process_antitbwjj->sigmaKin();
    matrix_elements = process_antitbwjj->getMatrixElements();
      weight = matrix_elements[0];
  }
  if (iMode==kAllPartonsTopLep || iMode==kAllPartonsTopLep_FixedTopM_FiniteWwidth){
    process_tbwlnu->setMomenta(*pTop);
    process_tbwlnu->sigmaKin();
    matrix_elements = process_tbwlnu->getMatrixElements();
      weight = matrix_elements[0];
  }
  if (iMode==kAllPartonsHiggsWWLep_FixedHiggsM_FiniteWwidth){
    process_hw2l2nu->setMomenta(*pHiggs);
    process_hw2l2nu->sigmaKin();
    matrix_elements = process_hw2l2nu->getMatrixElements();
      weight = matrix_elements[0];
  }

  if (verbosity>=2 && iGen==kMadgraph) {
 for(int i=0; i<1;i++)
    cout << " Matrix element = "
         << setiosflags(ios::fixed) //<< setprecision(17)
         << matrix_elements[i]
         << " GeV^" << -(2*nExternals-8) << endl;
  }
/*
  if (verbosity>=2 && iGen==kGosam) cout << " Matrix element = "
         << setiosflags(ios::fixed) //<< setprecision(17)
         << matrix_elements[3] * (4 * 3.1416 * 0.118 * 4 * 3.1416 * 0.118)
         << " GeV^" << -(2*nExternals-8) << endl;
*/
  return weight;
}

double MEPhaseSpace::ComputePDF(double x1, double x2, double mu) const {

  if (verbosity>=2) cout << "x1="<<x1<<" x2="<<x2<<" mu="<<mu<<endl;
  if(x1<0 || x1>0.99 || x2<0 || x2>0.99) return 0.;
  
  int pid = 0; 
  double xf1 = 0;
  double xf2 = 0;
  double Pdf = 0;

  //gluon only
  if (iMode==kMEM_TTH_TopAntitopHiggsDecay || iMode==kMEM_TTH_TopAntitopHiggsSemiLepDecay || iMode==kMEM_TTLL_TopAntitopDecay) {
    pid = 21;
    xf1 = pdf->xfxQ2(pid, x1, mu*mu) / x1;
    xf2 = pdf->xfxQ2(pid, x2, mu*mu) / x2;
    Pdf = 2 * xf1 * xf2;
  }

  double weightPDF = Pdf /x1 /x2 /4 /(comEnergy*comEnergy);
  if (verbosity>=2) {
    cout << "x1="<<x1<<" x2="<<x2<<" Pdf="<<Pdf<<" weightPDF="<<weightPDF<<endl;
  }

  return weightPDF;
}

double MEPhaseSpace::ConvolvePdfCrossSection(double x1, double x2, double mu) const{

  double weight = 0;
  if (iMode==kMEM_TTW_TopAntitopDecay && MEMFix_HiggsSemiLep.LepSign==-1){
    weight = pdf->xfxQ2(1, x1, mu*mu) * pdf->xfxQ2(-2, x2, mu*mu) * ComputeSubMatrixElement(kTTW, 1, -2) 
           + pdf->xfxQ2(-2, x1, mu*mu) * pdf->xfxQ2(1, x2, mu*mu) * ComputeSubMatrixElement(kTTW, -2, 1)
  	   + pdf->xfxQ2(3, x1, mu*mu) * pdf->xfxQ2(-4, x2, mu*mu) * ComputeSubMatrixElement(kTTW, 3, -4)
	   + pdf->xfxQ2(-4, x1, mu*mu) * pdf->xfxQ2(3, x2, mu*mu) * ComputeSubMatrixElement(kTTW, -4, 3)
	   + pdf->xfxQ2(1, x1, mu*mu) * pdf->xfxQ2(-4, x2, mu*mu) * ComputeSubMatrixElement(kTTW, 1, -4)
	   + pdf->xfxQ2(-4, x1, mu*mu) * pdf->xfxQ2(1, x2, mu*mu) * ComputeSubMatrixElement(kTTW, -4, 1)
	   + pdf->xfxQ2(-2, x1, mu*mu) * pdf->xfxQ2(3, x2, mu*mu) * ComputeSubMatrixElement(kTTW, -2, 3)
	   + pdf->xfxQ2(3, x1, mu*mu) * pdf->xfxQ2(-2, x2, mu*mu) * ComputeSubMatrixElement(kTTW, 3, -2);
  }
  if (iMode==kMEM_TTW_TopAntitopDecay && MEMFix_HiggsSemiLep.LepSign==1){
    weight = pdf->xfxQ2(-1, x1, mu*mu) * pdf->xfxQ2(2, x2, mu*mu) * ComputeSubMatrixElement(kTTW, -1, 2)       
           + pdf->xfxQ2(2, x1, mu*mu) * pdf->xfxQ2(-1, x2, mu*mu) * ComputeSubMatrixElement(kTTW, 2, -1)
           + pdf->xfxQ2(-3, x1, mu*mu) * pdf->xfxQ2(4, x2, mu*mu) * ComputeSubMatrixElement(kTTW, -3, 4)
           + pdf->xfxQ2(4, x1, mu*mu) * pdf->xfxQ2(-3, x2, mu*mu) * ComputeSubMatrixElement(kTTW, 4, -3)
           + pdf->xfxQ2(-1, x1, mu*mu) * pdf->xfxQ2(4, x2, mu*mu) * ComputeSubMatrixElement(kTTW, -1, 4)
           + pdf->xfxQ2(4, x1, mu*mu) * pdf->xfxQ2(-1, x2, mu*mu) * ComputeSubMatrixElement(kTTW, 4, -1)
           + pdf->xfxQ2(2, x1, mu*mu) * pdf->xfxQ2(-3, x2, mu*mu) * ComputeSubMatrixElement(kTTW, 2, -3)
           + pdf->xfxQ2(-3, x1, mu*mu) * pdf->xfxQ2(2, x2, mu*mu) * ComputeSubMatrixElement(kTTW, -3, 2);
  }

  if (iMode==kMEM_TTWJJ_TopAntitopDecay && MEMFix_HiggsSemiLep.LepSign==-1){
    weight =// pdf->xfxQ2(-4, x1, mu*mu) * pdf->xfxQ2(21, x2, mu*mu) * ComputeSubMatrixElement(kTTWJJ, -4, 21)
	   //+ pdf->xfxQ2(21, x1, mu*mu) * pdf->xfxQ2(-4, x2, mu*mu) * ComputeSubMatrixElement(kTTWJJ, 21, -4)
	    pdf->xfxQ2(21, x1, mu*mu) * pdf->xfxQ2(1, x2, mu*mu) * ComputeSubMatrixElement(kTTWJJ, 21, 1)
           + pdf->xfxQ2(1, x1, mu*mu) * pdf->xfxQ2(21, x2, mu*mu) * ComputeSubMatrixElement(kTTWJJ, 1, 21);
	   //+ pdf->xfxQ2(-2, x1, mu*mu) * pdf->xfxQ2(21, x2, mu*mu) * ComputeSubMatrixElement(kTTWJJ, -2, 21)
           //+ pdf->xfxQ2(21, x1, mu*mu) * pdf->xfxQ2(-2, x2, mu*mu) * ComputeSubMatrixElement(kTTWJJ, 21, -2)
           //+ pdf->xfxQ2(21, x1, mu*mu) * pdf->xfxQ2(3, x2, mu*mu) * ComputeSubMatrixElement(kTTWJJ, 21, 3)
           //+ pdf->xfxQ2(3, x1, mu*mu) * pdf->xfxQ2(21, x2, mu*mu) * ComputeSubMatrixElement(kTTWJJ, 3, 21);
  }
  if (iMode==kMEM_TTWJJ_TopAntitopDecay && MEMFix_HiggsSemiLep.LepSign==1){
    weight = //pdf->xfxQ2(4, x1, mu*mu) * pdf->xfxQ2(21, x2, mu*mu) * ComputeSubMatrixElement(kTTWJJ, 4, 21)
           //+ pdf->xfxQ2(21, x1, mu*mu) * pdf->xfxQ2(4, x2, mu*mu) * ComputeSubMatrixElement(kTTWJJ, 21, 4)
           //+ pdf->xfxQ2(21, x1, mu*mu) * pdf->xfxQ2(-1, x2, mu*mu) * ComputeSubMatrixElement(kTTWJJ, 21, -1)
           //+ pdf->xfxQ2(-1, x1, mu*mu) * pdf->xfxQ2(21, x2, mu*mu) * ComputeSubMatrixElement(kTTWJJ, -1, 21)
            pdf->xfxQ2(2, x1, mu*mu) * pdf->xfxQ2(21, x2, mu*mu) * ComputeSubMatrixElement(kTTWJJ, 2, 21)
           + pdf->xfxQ2(21, x1, mu*mu) * pdf->xfxQ2(2, x2, mu*mu) * ComputeSubMatrixElement(kTTWJJ, 21, 2);
           //+ pdf->xfxQ2(21, x1, mu*mu) * pdf->xfxQ2(-3, x2, mu*mu) * ComputeSubMatrixElement(kTTWJJ, 21, -3)
           //+ pdf->xfxQ2(-3, x1, mu*mu) * pdf->xfxQ2(21, x2, mu*mu) * ComputeSubMatrixElement(kTTWJJ, -3, 21);
  }
  if (iMode==kMEM_TTW_TopAntitopDecay || iMode==kMEM_TTWJJ_TopAntitopDecay){
    weight *= ComputeSubMatrixElement(kTop, 0, 0);
    weight *= ComputeSubMatrixElement(kAntitop, 0, 0);
  }
 
  weight = weight /x1/x1 /x2/x2 /4 /(comEnergy*comEnergy);

  if (verbosity>=2) cout << "Convolve weightMEPDF="<<weight<<endl; 

  return weight;
}

double MEPhaseSpace::ComputeTFProduct(double Gen_Bjet1_E, double Gen_Bjet2_E, double Gen_Jet_E, double Gen_Jet2_E, double Gen_Recoil_Px, double Gen_Recoil_Py) const {

  double weightTF_Bjet1 = 1;
  double weightTF_Bjet2 = 1;
  double weightTF_Jet = 1;
  double weightTF_Jet2 = 1;
  double weightTF_Recoil_Px = 1;
  double weightTF_Recoil_Py = 1;

  if (iTF == kTFGaussian){
    weightTF_Bjet1 = ComputeSingleTF_Gaus(Gen_Bjet1_E, MeasuredVarForTF.Bjet1_E, 0.2);
    weightTF_Bjet2 = ComputeSingleTF_Gaus(Gen_Bjet2_E, MeasuredVarForTF.Bjet2_E, 0.2);
    if (iMode!=kMEM_TTW_TopAntitopDecay) weightTF_Jet = ComputeSingleTF_Gaus(Gen_Jet_E, MeasuredVarForTF.Jet1_E, 0.2);
    if (iMode==kMEM_TTH_TopAntitopHiggsSemiLepDecay || iMode==kMEM_TTWJJ_TopAntitopDecay) weightTF_Jet2 = ComputeSingleTF_Gaus(Gen_Jet2_E, MeasuredVarForTF.Jet2_E, 0.2);
    //weightTF_Jet2 = ComputeSingleTF_Gaus(Gen_Jet2_E, MeasuredVarForTF.Jet2_E, 0.2);
    if (iTFOption==kTFRecoilPtot){
      weightTF_Recoil_Px = ComputeSingleTF_Gaus(Gen_Recoil_Px, MeasuredVarForTF.Recoil_Px, 0.2);
      weightTF_Recoil_Py = ComputeSingleTF_Gaus(Gen_Recoil_Py, MeasuredVarForTF.Recoil_Py, 0.2);
    }
    else if (iTFOption==kTFRecoilmET){
      weightTF_Recoil_Px = ComputeSingleTF_Gaus(ComputedVarForTF.mET_Px, MeasuredVarForTF.mET_Px, 0.2);
      weightTF_Recoil_Py = ComputeSingleTF_Gaus(ComputedVarForTF.mET_Py, MeasuredVarForTF.mET_Py, 0.2);
    }
  }
  if (iTF == kTFHistoBnonB_GausRecoil){
    weightTF_Bjet1 = ComputeSingleTF_Histo("Bjet", Gen_Bjet1_E, MeasuredVarForTF.Bjet1_E, MeasuredVarForTF.Bjet1_Eta);
    weightTF_Bjet2 = ComputeSingleTF_Histo("Bjet", Gen_Bjet2_E, MeasuredVarForTF.Bjet2_E, MeasuredVarForTF.Bjet2_Eta);
    if (iMode!=kMEM_TTW_TopAntitopDecay) weightTF_Jet = ComputeSingleTF_Histo("Jet", Gen_Jet_E, MeasuredVarForTF.Jet1_E, MeasuredVarForTF.Jet1_Eta);
    if (iMode==kMEM_TTH_TopAntitopHiggsSemiLepDecay || iMode==kMEM_TTWJJ_TopAntitopDecay) weightTF_Jet2 = ComputeSingleTF_Histo("Jet", Gen_Jet2_E, MeasuredVarForTF.Jet2_E, MeasuredVarForTF.Jet2_Eta); 
    if (iTFOption==kTFRecoilPtot){
      weightTF_Recoil_Px = ComputeSingleTF_Gaus(Gen_Recoil_Px, MeasuredVarForTF.Recoil_Px, 0.2);
      weightTF_Recoil_Py = ComputeSingleTF_Gaus(Gen_Recoil_Py, MeasuredVarForTF.Recoil_Py, 0.2);
    }
    else if (iTFOption==kTFRecoilmET){
      weightTF_Recoil_Px = ComputeSingleTF_Gaus(ComputedVarForTF.mET_Px, MeasuredVarForTF.mET_Px, 0.2);      
      weightTF_Recoil_Py = ComputeSingleTF_Gaus(ComputedVarForTF.mET_Py, MeasuredVarForTF.mET_Py, 0.2);
    }
  }
  if (iTF == kTFHistoBnonBmET){
    weightTF_Bjet1 = ComputeSingleTF_Histo("Bjet", Gen_Bjet1_E, MeasuredVarForTF.Bjet1_E, MeasuredVarForTF.Bjet1_Eta);
    weightTF_Bjet2 = ComputeSingleTF_Histo("Bjet", Gen_Bjet2_E, MeasuredVarForTF.Bjet2_E, MeasuredVarForTF.Bjet2_Eta);
    if (iMode!=kMEM_TTW_TopAntitopDecay) weightTF_Jet = ComputeSingleTF_Histo("Jet", Gen_Jet_E, MeasuredVarForTF.Jet1_E, MeasuredVarForTF.Jet1_Eta);
    if (iMode==kMEM_TTH_TopAntitopHiggsSemiLepDecay || iMode==kMEM_TTWJJ_TopAntitopDecay) weightTF_Jet2 = ComputeSingleTF_Histo("Jet", Gen_Jet2_E, MeasuredVarForTF.Jet2_E, MeasuredVarForTF.Jet2_Eta);
    if (iTFOption==kTFRecoilmET){
      weightTF_Recoil_Px = ComputeSingleTF_Histo("mET_Px", ComputedVarForTF.mET_Px, MeasuredVarForTF.mET_Px, 0);
      weightTF_Recoil_Py = ComputeSingleTF_Histo("mET_Py", ComputedVarForTF.mET_Py, MeasuredVarForTF.mET_Py, 0);
    }
  }

  double weightTF = weightTF_Bjet1*weightTF_Bjet2*weightTF_Jet*weightTF_Jet2*weightTF_Recoil_Px*weightTF_Recoil_Py;

  if (verbosity>=2) cout << "TFProduct weightTF="<<weightTF<<endl; 

  return weightTF;
}

double MEPhaseSpace::ComputeSingleTF_Gaus(double Gen_E, double Reco_E, double Res) const {

  double weightTF = TMath::Gaus(Gen_E-Reco_E, 0, TMath::Abs(Res*Gen_E), kTRUE);
  if (verbosity>=2) cout << "TF Gaus Gen="<<Gen_E<<" Reco="<<Reco_E<<" weightTF="<<weightTF<<endl;

  return weightTF;
}

double MEPhaseSpace::ComputeSingleTF_Histo(string Part, double Gen_E, double Reco_E, double Reco_Eta) const {

  double weightTF = 1;

  float Jet_Pt = Reco_E / TMath::CosH(Reco_Eta);
  if ((Part=="Bjet" || Part=="Jet") && Jet_Pt<25) {
    if (verbosity>=2) cout << "TF Histo JetPt="<<Jet_Pt<<"<20, set weightTF=1"<<endl;
    return 1;
  }

  int iEta = -1;
  for (int iBin=0; iBin<3; iBin++){
    if (fabs(Reco_Eta) > EtaRange[iBin] && fabs(Reco_Eta) < EtaRange[iBin+1]) iEta=iBin;
  }
  int iEnergy = -1;
  for (int iBin=0; iBin<6; iBin++){
    if (Reco_E > EnergyRange[iBin] && Reco_E < EnergyRange[iBin+1]) iEnergy=iBin;
  }
  //if (Reco_E > EnergyRange[6]) iEnergy = 6;
  if ((Part=="Bjet" || Part=="Jet") && (iEta==-1 || iEnergy==-1)) {
    if (verbosity>=2) cout << "TF Histo Jet Energy="<<Reco_E<<" |Eta|="<<fabs(Reco_Eta)<<" out of range, set weightTF=1"<<endl;
    return 1;
  }
  //cout << "ComputeSingleTF_Histo iEta="<<iEta<<" iEnergy="<<iEnergy<<endl;

  int iMet = -1;
  for (int iBin=0; iBin<3; iBin++){
    if (fabs(Reco_E) > MetRange[iBin] && fabs(Reco_E)<MetRange[iBin+1]) iMet = iBin;
  }
  if ((Part=="mET_Px" || Part=="mET_Py") && iMet==-1){
    if (verbosity>=2) cout << "TF Histo mET component = "<<Reco_E<<" out of range, set weightTF=1"<<endl;
    return 1;
  }

  int iTFbin = -1;
  if (Part=="Bjet") {
    iTFbin = TF_B[iEta][iEnergy]->FindBin(Gen_E-Reco_E);
    if (iTFbin >= TF_B[iEta][iEnergy]->GetNbinsX() || iTFbin==0) weightTF = 0;
    else weightTF = TF_B[iEta][iEnergy]->GetBinContent(iTFbin);
  }
  if (Part=="Jet") {
    iTFbin = TF_nonB[iEta][iEnergy]->FindBin(Gen_E-Reco_E);
    if (iTFbin >= TF_nonB[iEta][iEnergy]->GetNbinsX() || iTFbin==0) weightTF = 0;
    weightTF = TF_nonB[iEta][iEnergy]->GetBinContent(iTFbin);
  }
  if (Part=="mET_Px") {
    iTFbin = TF_MetPx[iMet]->FindBin(Gen_E-Reco_E);
    if (iTFbin >= TF_MetPx[iMet]->GetNbinsX() || iTFbin==0) weightTF = 0;
    weightTF = TF_MetPx[iMet]->GetBinContent(iTFbin);
  }
  if (Part=="mET_Py") { 
    iTFbin = TF_MetPy[iMet]->FindBin(Gen_E-Reco_E);
    if (iTFbin >= TF_MetPy[iMet]->GetNbinsX() || iTFbin==0) weightTF = 0;
    weightTF = TF_MetPy[iMet]->GetBinContent(iTFbin);
  }

  if (verbosity>=2) cout << "TF Histo "<<Part<< " Gen="<<Gen_E<<" Reco="<<Reco_E<<" weightTF="<<weightTF<<endl;

  return weightTF;
}

void MEPhaseSpace::LoadTFfromHisto(string FileName) const {

  TFile* fTF = new TFile(FileName.c_str(),"READ"); 
  //fTF->ls();
  //fTF->SetDirectory(gROOT);
  //fTF->GetListOfKeys()->Print();
 
  gDirectory->pwd();

  TH1F* Histo_tmp;
  string histoname;

  //Loading histos
  for (int iEta=0; iEta<3; iEta++){
    for (int iEnergy=0; iEnergy<6; iEnergy++){

      histoname = "TF_nonB_" + EtaRangeLabel[iEta] + "_" + EnergyRangeLabel[iEnergy];
      cout << histoname << endl;
      //Histo_tmp = (TH1F*) fTF->Get(histoname.c_str());
      //Histo_tmp = (TH1F*) gDirectory->GetList()->FindObject(histoname.c_str());
      //Histo_tmp = (TH1F*) gDirectory->Get(histoname.c_str());
      fTF->GetObject(histoname.c_str(), Histo_tmp);
      //cout << "Nbins="<<Histo_tmp->GetNbinsX()<<endl;
      //for (int iBin=1; iBin<100;  iBin++) cout << "Bin "<<iBin<<" val="<<Histo_tmp->GetBinContent(iBin)<<endl;


      TF_nonB[iEta][iEnergy] = (TH1F*) Histo_tmp->Clone();

      //for (int iBin=0; iBin<TF_nonB[iEta][iEnergy]->GetNbinsX(); iBin++) cout << "Bin "<<iBin<<" val="<<TF_nonB[iEta][iEnergy]->GetBinContent(iBin)<<endl;

      histoname = "TF_B_" + EtaRangeLabel[iEta] + "_" + EnergyRangeLabel[iEnergy];
      Histo_tmp = (TH1F*) fTF->Get(histoname.c_str());
      TF_B[iEta][iEnergy] = (TH1F*) Histo_tmp->Clone();

      //for (int iBin=0; iBin<TF_B[iEta][iEnergy]->GetNbinsX(); iBin++) cout << "Bin "<<iBin<<" val="<<TF_B[iEta][iEnergy]->GetBinContent(iBin)<<endl;

    }
  }

  for (int iMet=0; iMet<3; iMet++){
    histoname = "TF_mET_Px_" + MetRangeLabel[iMet];
    Histo_tmp = (TH1F*) fTF->Get(histoname.c_str());
    TF_MetPx[iMet] = (TH1F*) Histo_tmp->Clone();

    histoname = "TF_mET_Py_" + MetRangeLabel[iMet];
    Histo_tmp = (TH1F*) fTF->Get(histoname.c_str());
    TF_MetPy[iMet] = (TH1F*) Histo_tmp->Clone();
  }


  //fTF->Close();

  return;
}

vector<double*> MEPhaseSpace::GetPhaseSpacePoint(){

  double effectiveLumi = (float)(rand() % 10000) / 10000. * comEnergy;
    if (verbosity>=1) cout << "effectiveLumi="<<effectiveLumi<<endl;

  double weight;
  vector<double*> part = get_momenta(process->ninitial, effectiveLumi,
                                 process->getMasses(), weight);

  return part;
}

void MEPhaseSpace::CheckMatrixElement(){
  
  vector<double*> part = GetPhaseSpacePoint();

  bool CheckDirectly = false;

  if (CheckDirectly){

  process->setMomenta(part);

  process->sigmaKin();

  const double* matrix_elements = process->getMatrixElements();

  if (verbosity>=1) {
  cout << "Momenta:" << endl;
  for(int i=0;i < process->nexternal; i++)
    cout << setprecision(10) << i+1
         << setiosflags(ios::scientific) << setw(14) << part[i][0]
         << setiosflags(ios::scientific) << setw(14) << part[i][1]
         << setiosflags(ios::scientific) << setw(14) << part[i][2]
         << setiosflags(ios::scientific) << setw(14) << part[i][3] << endl;

  for(int i=0; i<process->nprocesses;i++)
    cout << " Matrix element = "
         << setiosflags(ios::fixed) << setprecision(17)
         << matrix_elements[i]
         << " GeV^" << -(2*process->nexternal-8) << endl;

  }
  }
  else {

  (*pCore)[2] = part[2];
  (*pCore)[3] = part[3];
 
  double x[2];
  x[0] = 1;
  x[1] = 1;

  Eval(x);

    if (verbosity>=1) {
  cout << "Momenta:" << endl;
  for(int i=0;i < process->nexternal; i++)
    cout << setprecision(10) << i+1
         << setiosflags(ios::scientific) << setw(14) << (*pCore)[i][0]
         << setiosflags(ios::scientific) << setw(14) << (*pCore)[i][1]
         << setiosflags(ios::scientific) << setw(14) << (*pCore)[i][2]
         << setiosflags(ios::scientific) << setw(14) << (*pCore)[i][3] << endl;
  }


  }
}

void MEPhaseSpace::SetVerbosity(int v=1){
  verbosity = v;
  return;
}

#endif
