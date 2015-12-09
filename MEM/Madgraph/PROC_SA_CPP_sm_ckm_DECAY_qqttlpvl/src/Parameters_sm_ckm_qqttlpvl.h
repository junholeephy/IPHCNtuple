//==========================================================================
// This file has been automatically generated for C++
// MadGraph5_aMC@NLO v. 2.3.0, 2015-07-01
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef Parameters_sm_ckm_qqttlpvl_H
#define Parameters_sm_ckm_qqttlpvl_H

#include <complex> 

#include "../../PROC_SA_CPP_sm_4/src/read_slha.h"
using namespace std; 

class Parameters_sm_ckm_qqttlpvl
{
  public:

    static Parameters_sm_ckm_qqttlpvl * getInstance(); 

    // Define "zero"
    double zero, ZERO; 
    // Model parameters independent of aS
    double mdl_WH, mdl_WW, mdl_WZ, mdl_WT, mdl_ymtau, mdl_ymt, mdl_ymb,
        mdl_etaWS, mdl_rhoWS, mdl_AWS, mdl_lamWS, aS, mdl_Gf, aEWM1, mdl_MH,
        mdl_MZ, mdl_MTA, mdl_MT, mdl_MB, mdl_CKM3x3, mdl_conjg__CKM3x3,
        mdl_lamWS__exp__2, mdl_lamWS__exp__3, mdl_MZ__exp__2, mdl_MZ__exp__4,
        mdl_sqrt__2, mdl_MH__exp__2, mdl_aEW, mdl_MW, mdl_sqrt__aEW, mdl_ee,
        mdl_MW__exp__2, mdl_sw2, mdl_cw, mdl_sqrt__sw2, mdl_sw, mdl_g1, mdl_gw,
        mdl_vev, mdl_vev__exp__2, mdl_lam, mdl_yb, mdl_yt, mdl_ytau, mdl_muH,
        mdl_ee__exp__2, mdl_sw__exp__2, mdl_cw__exp__2;
    std::complex<double> mdl_CKM1x1, mdl_CKM1x2, mdl_complexi, mdl_CKM1x3,
        mdl_CKM2x1, mdl_CKM2x2, mdl_CKM2x3, mdl_CKM3x1, mdl_CKM3x2,
        mdl_conjg__CKM1x3, mdl_conjg__CKM2x3, mdl_conjg__CKM2x1,
        mdl_conjg__CKM3x1, mdl_conjg__CKM2x2, mdl_conjg__CKM3x2,
        mdl_conjg__CKM1x1, mdl_conjg__CKM1x2, mdl_I1x31, mdl_I1x32, mdl_I1x33,
        mdl_I2x13, mdl_I2x23, mdl_I2x33, mdl_I3x31, mdl_I3x32, mdl_I3x33,
        mdl_I4x13, mdl_I4x23, mdl_I4x33;
    // Model parameters dependent on aS
    double mdl_sqrt__aS, G, mdl_G__exp__2; 
    // Model couplings independent of aS
    std::complex<double> GC_44, GC_100, GC_101, GC_108; 
    // Model couplings dependent on aS
    std::complex<double> GC_11; 

    // Set parameters that are unchanged during the run
    void setIndependentParameters(SLHAReader& slha); 
    // Set couplings that are unchanged during the run
    void setIndependentCouplings(); 
    // Set parameters that are changed event by event
    void setDependentParameters(); 
    // Set couplings that are changed event by event
    void setDependentCouplings(); 

    // Print parameters that are unchanged during the run
    void printIndependentParameters(); 
    // Print couplings that are unchanged during the run
    void printIndependentCouplings(); 
    // Print parameters that are changed event by event
    void printDependentParameters(); 
    // Print couplings that are changed event by event
    void printDependentCouplings(); 


  private:
    static Parameters_sm_ckm_qqttlpvl * instance; 
}; 

#endif  // Parameters_sm_ckm_qqttlpvl_H

