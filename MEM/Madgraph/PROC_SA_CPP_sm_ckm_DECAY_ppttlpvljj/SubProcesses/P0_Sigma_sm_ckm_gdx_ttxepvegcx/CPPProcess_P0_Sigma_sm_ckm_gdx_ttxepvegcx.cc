//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.3.0, 2015-07-01
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "CPPProcess_P0_Sigma_sm_ckm_gdx_ttxepvegcx.h"
#include "HelAmps_sm_ckm_ppttlpvljj.h"

using namespace MG5_sm_ckm_ppttlpvljj; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g d~ > t t~ e+ ve g c~ WEIGHTED=8
// Process: g d~ > t t~ mu+ vm g c~ WEIGHTED=8

//--------------------------------------------------------------------------
// Initialize process.

void CPPProcess_P0_Sigma_sm_ckm_gdx_ttxepvegcx::initProc(string param_card_name) 
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_sm_ckm_ppttlpvljj::getInstance(); 
  SLHAReader slha(param_card_name); 
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings(); 
  pars->printIndependentParameters(); 
  pars->printIndependentCouplings(); 
  // Set external particle masses for this matrix element
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->mdl_MT); 
  mME.push_back(pars->mdl_MT); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  jamp2[0] = new double[12]; 
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void CPPProcess_P0_Sigma_sm_ckm_gdx_ttxepvegcx::sigmaKin() 
{
  // Set the parameters which change event by event
  pars->setDependentParameters(); 
  pars->setDependentCouplings(); 
  static bool firsttime = true; 
  if (firsttime)
  {
    pars->printDependentParameters(); 
    pars->printDependentCouplings(); 
    firsttime = false; 
  }

  // Reset color flows
  for(int i = 0; i < 12; i++ )
    jamp2[0][i] = 0.; 

  // Local variables and constants
  const int ncomb = 256; 
  static bool goodhel[ncomb] = {ncomb * false}; 
  static int ntry = 0, sum_hel = 0, ngood = 0; 
  static int igood[ncomb]; 
  static int jhel; 
  std::complex<double> * * wfs; 
  double t[nprocesses]; 
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1, -1, -1, -1,
      -1}, {-1, -1, -1, -1, -1, -1, -1, 1}, {-1, -1, -1, -1, -1, -1, 1, -1},
      {-1, -1, -1, -1, -1, -1, 1, 1}, {-1, -1, -1, -1, -1, 1, -1, -1}, {-1, -1,
      -1, -1, -1, 1, -1, 1}, {-1, -1, -1, -1, -1, 1, 1, -1}, {-1, -1, -1, -1,
      -1, 1, 1, 1}, {-1, -1, -1, -1, 1, -1, -1, -1}, {-1, -1, -1, -1, 1, -1,
      -1, 1}, {-1, -1, -1, -1, 1, -1, 1, -1}, {-1, -1, -1, -1, 1, -1, 1, 1},
      {-1, -1, -1, -1, 1, 1, -1, -1}, {-1, -1, -1, -1, 1, 1, -1, 1}, {-1, -1,
      -1, -1, 1, 1, 1, -1}, {-1, -1, -1, -1, 1, 1, 1, 1}, {-1, -1, -1, 1, -1,
      -1, -1, -1}, {-1, -1, -1, 1, -1, -1, -1, 1}, {-1, -1, -1, 1, -1, -1, 1,
      -1}, {-1, -1, -1, 1, -1, -1, 1, 1}, {-1, -1, -1, 1, -1, 1, -1, -1}, {-1,
      -1, -1, 1, -1, 1, -1, 1}, {-1, -1, -1, 1, -1, 1, 1, -1}, {-1, -1, -1, 1,
      -1, 1, 1, 1}, {-1, -1, -1, 1, 1, -1, -1, -1}, {-1, -1, -1, 1, 1, -1, -1,
      1}, {-1, -1, -1, 1, 1, -1, 1, -1}, {-1, -1, -1, 1, 1, -1, 1, 1}, {-1, -1,
      -1, 1, 1, 1, -1, -1}, {-1, -1, -1, 1, 1, 1, -1, 1}, {-1, -1, -1, 1, 1, 1,
      1, -1}, {-1, -1, -1, 1, 1, 1, 1, 1}, {-1, -1, 1, -1, -1, -1, -1, -1},
      {-1, -1, 1, -1, -1, -1, -1, 1}, {-1, -1, 1, -1, -1, -1, 1, -1}, {-1, -1,
      1, -1, -1, -1, 1, 1}, {-1, -1, 1, -1, -1, 1, -1, -1}, {-1, -1, 1, -1, -1,
      1, -1, 1}, {-1, -1, 1, -1, -1, 1, 1, -1}, {-1, -1, 1, -1, -1, 1, 1, 1},
      {-1, -1, 1, -1, 1, -1, -1, -1}, {-1, -1, 1, -1, 1, -1, -1, 1}, {-1, -1,
      1, -1, 1, -1, 1, -1}, {-1, -1, 1, -1, 1, -1, 1, 1}, {-1, -1, 1, -1, 1, 1,
      -1, -1}, {-1, -1, 1, -1, 1, 1, -1, 1}, {-1, -1, 1, -1, 1, 1, 1, -1}, {-1,
      -1, 1, -1, 1, 1, 1, 1}, {-1, -1, 1, 1, -1, -1, -1, -1}, {-1, -1, 1, 1,
      -1, -1, -1, 1}, {-1, -1, 1, 1, -1, -1, 1, -1}, {-1, -1, 1, 1, -1, -1, 1,
      1}, {-1, -1, 1, 1, -1, 1, -1, -1}, {-1, -1, 1, 1, -1, 1, -1, 1}, {-1, -1,
      1, 1, -1, 1, 1, -1}, {-1, -1, 1, 1, -1, 1, 1, 1}, {-1, -1, 1, 1, 1, -1,
      -1, -1}, {-1, -1, 1, 1, 1, -1, -1, 1}, {-1, -1, 1, 1, 1, -1, 1, -1}, {-1,
      -1, 1, 1, 1, -1, 1, 1}, {-1, -1, 1, 1, 1, 1, -1, -1}, {-1, -1, 1, 1, 1,
      1, -1, 1}, {-1, -1, 1, 1, 1, 1, 1, -1}, {-1, -1, 1, 1, 1, 1, 1, 1}, {-1,
      1, -1, -1, -1, -1, -1, -1}, {-1, 1, -1, -1, -1, -1, -1, 1}, {-1, 1, -1,
      -1, -1, -1, 1, -1}, {-1, 1, -1, -1, -1, -1, 1, 1}, {-1, 1, -1, -1, -1, 1,
      -1, -1}, {-1, 1, -1, -1, -1, 1, -1, 1}, {-1, 1, -1, -1, -1, 1, 1, -1},
      {-1, 1, -1, -1, -1, 1, 1, 1}, {-1, 1, -1, -1, 1, -1, -1, -1}, {-1, 1, -1,
      -1, 1, -1, -1, 1}, {-1, 1, -1, -1, 1, -1, 1, -1}, {-1, 1, -1, -1, 1, -1,
      1, 1}, {-1, 1, -1, -1, 1, 1, -1, -1}, {-1, 1, -1, -1, 1, 1, -1, 1}, {-1,
      1, -1, -1, 1, 1, 1, -1}, {-1, 1, -1, -1, 1, 1, 1, 1}, {-1, 1, -1, 1, -1,
      -1, -1, -1}, {-1, 1, -1, 1, -1, -1, -1, 1}, {-1, 1, -1, 1, -1, -1, 1,
      -1}, {-1, 1, -1, 1, -1, -1, 1, 1}, {-1, 1, -1, 1, -1, 1, -1, -1}, {-1, 1,
      -1, 1, -1, 1, -1, 1}, {-1, 1, -1, 1, -1, 1, 1, -1}, {-1, 1, -1, 1, -1, 1,
      1, 1}, {-1, 1, -1, 1, 1, -1, -1, -1}, {-1, 1, -1, 1, 1, -1, -1, 1}, {-1,
      1, -1, 1, 1, -1, 1, -1}, {-1, 1, -1, 1, 1, -1, 1, 1}, {-1, 1, -1, 1, 1,
      1, -1, -1}, {-1, 1, -1, 1, 1, 1, -1, 1}, {-1, 1, -1, 1, 1, 1, 1, -1},
      {-1, 1, -1, 1, 1, 1, 1, 1}, {-1, 1, 1, -1, -1, -1, -1, -1}, {-1, 1, 1,
      -1, -1, -1, -1, 1}, {-1, 1, 1, -1, -1, -1, 1, -1}, {-1, 1, 1, -1, -1, -1,
      1, 1}, {-1, 1, 1, -1, -1, 1, -1, -1}, {-1, 1, 1, -1, -1, 1, -1, 1}, {-1,
      1, 1, -1, -1, 1, 1, -1}, {-1, 1, 1, -1, -1, 1, 1, 1}, {-1, 1, 1, -1, 1,
      -1, -1, -1}, {-1, 1, 1, -1, 1, -1, -1, 1}, {-1, 1, 1, -1, 1, -1, 1, -1},
      {-1, 1, 1, -1, 1, -1, 1, 1}, {-1, 1, 1, -1, 1, 1, -1, -1}, {-1, 1, 1, -1,
      1, 1, -1, 1}, {-1, 1, 1, -1, 1, 1, 1, -1}, {-1, 1, 1, -1, 1, 1, 1, 1},
      {-1, 1, 1, 1, -1, -1, -1, -1}, {-1, 1, 1, 1, -1, -1, -1, 1}, {-1, 1, 1,
      1, -1, -1, 1, -1}, {-1, 1, 1, 1, -1, -1, 1, 1}, {-1, 1, 1, 1, -1, 1, -1,
      -1}, {-1, 1, 1, 1, -1, 1, -1, 1}, {-1, 1, 1, 1, -1, 1, 1, -1}, {-1, 1, 1,
      1, -1, 1, 1, 1}, {-1, 1, 1, 1, 1, -1, -1, -1}, {-1, 1, 1, 1, 1, -1, -1,
      1}, {-1, 1, 1, 1, 1, -1, 1, -1}, {-1, 1, 1, 1, 1, -1, 1, 1}, {-1, 1, 1,
      1, 1, 1, -1, -1}, {-1, 1, 1, 1, 1, 1, -1, 1}, {-1, 1, 1, 1, 1, 1, 1, -1},
      {-1, 1, 1, 1, 1, 1, 1, 1}, {1, -1, -1, -1, -1, -1, -1, -1}, {1, -1, -1,
      -1, -1, -1, -1, 1}, {1, -1, -1, -1, -1, -1, 1, -1}, {1, -1, -1, -1, -1,
      -1, 1, 1}, {1, -1, -1, -1, -1, 1, -1, -1}, {1, -1, -1, -1, -1, 1, -1, 1},
      {1, -1, -1, -1, -1, 1, 1, -1}, {1, -1, -1, -1, -1, 1, 1, 1}, {1, -1, -1,
      -1, 1, -1, -1, -1}, {1, -1, -1, -1, 1, -1, -1, 1}, {1, -1, -1, -1, 1, -1,
      1, -1}, {1, -1, -1, -1, 1, -1, 1, 1}, {1, -1, -1, -1, 1, 1, -1, -1}, {1,
      -1, -1, -1, 1, 1, -1, 1}, {1, -1, -1, -1, 1, 1, 1, -1}, {1, -1, -1, -1,
      1, 1, 1, 1}, {1, -1, -1, 1, -1, -1, -1, -1}, {1, -1, -1, 1, -1, -1, -1,
      1}, {1, -1, -1, 1, -1, -1, 1, -1}, {1, -1, -1, 1, -1, -1, 1, 1}, {1, -1,
      -1, 1, -1, 1, -1, -1}, {1, -1, -1, 1, -1, 1, -1, 1}, {1, -1, -1, 1, -1,
      1, 1, -1}, {1, -1, -1, 1, -1, 1, 1, 1}, {1, -1, -1, 1, 1, -1, -1, -1},
      {1, -1, -1, 1, 1, -1, -1, 1}, {1, -1, -1, 1, 1, -1, 1, -1}, {1, -1, -1,
      1, 1, -1, 1, 1}, {1, -1, -1, 1, 1, 1, -1, -1}, {1, -1, -1, 1, 1, 1, -1,
      1}, {1, -1, -1, 1, 1, 1, 1, -1}, {1, -1, -1, 1, 1, 1, 1, 1}, {1, -1, 1,
      -1, -1, -1, -1, -1}, {1, -1, 1, -1, -1, -1, -1, 1}, {1, -1, 1, -1, -1,
      -1, 1, -1}, {1, -1, 1, -1, -1, -1, 1, 1}, {1, -1, 1, -1, -1, 1, -1, -1},
      {1, -1, 1, -1, -1, 1, -1, 1}, {1, -1, 1, -1, -1, 1, 1, -1}, {1, -1, 1,
      -1, -1, 1, 1, 1}, {1, -1, 1, -1, 1, -1, -1, -1}, {1, -1, 1, -1, 1, -1,
      -1, 1}, {1, -1, 1, -1, 1, -1, 1, -1}, {1, -1, 1, -1, 1, -1, 1, 1}, {1,
      -1, 1, -1, 1, 1, -1, -1}, {1, -1, 1, -1, 1, 1, -1, 1}, {1, -1, 1, -1, 1,
      1, 1, -1}, {1, -1, 1, -1, 1, 1, 1, 1}, {1, -1, 1, 1, -1, -1, -1, -1}, {1,
      -1, 1, 1, -1, -1, -1, 1}, {1, -1, 1, 1, -1, -1, 1, -1}, {1, -1, 1, 1, -1,
      -1, 1, 1}, {1, -1, 1, 1, -1, 1, -1, -1}, {1, -1, 1, 1, -1, 1, -1, 1}, {1,
      -1, 1, 1, -1, 1, 1, -1}, {1, -1, 1, 1, -1, 1, 1, 1}, {1, -1, 1, 1, 1, -1,
      -1, -1}, {1, -1, 1, 1, 1, -1, -1, 1}, {1, -1, 1, 1, 1, -1, 1, -1}, {1,
      -1, 1, 1, 1, -1, 1, 1}, {1, -1, 1, 1, 1, 1, -1, -1}, {1, -1, 1, 1, 1, 1,
      -1, 1}, {1, -1, 1, 1, 1, 1, 1, -1}, {1, -1, 1, 1, 1, 1, 1, 1}, {1, 1, -1,
      -1, -1, -1, -1, -1}, {1, 1, -1, -1, -1, -1, -1, 1}, {1, 1, -1, -1, -1,
      -1, 1, -1}, {1, 1, -1, -1, -1, -1, 1, 1}, {1, 1, -1, -1, -1, 1, -1, -1},
      {1, 1, -1, -1, -1, 1, -1, 1}, {1, 1, -1, -1, -1, 1, 1, -1}, {1, 1, -1,
      -1, -1, 1, 1, 1}, {1, 1, -1, -1, 1, -1, -1, -1}, {1, 1, -1, -1, 1, -1,
      -1, 1}, {1, 1, -1, -1, 1, -1, 1, -1}, {1, 1, -1, -1, 1, -1, 1, 1}, {1, 1,
      -1, -1, 1, 1, -1, -1}, {1, 1, -1, -1, 1, 1, -1, 1}, {1, 1, -1, -1, 1, 1,
      1, -1}, {1, 1, -1, -1, 1, 1, 1, 1}, {1, 1, -1, 1, -1, -1, -1, -1}, {1, 1,
      -1, 1, -1, -1, -1, 1}, {1, 1, -1, 1, -1, -1, 1, -1}, {1, 1, -1, 1, -1,
      -1, 1, 1}, {1, 1, -1, 1, -1, 1, -1, -1}, {1, 1, -1, 1, -1, 1, -1, 1}, {1,
      1, -1, 1, -1, 1, 1, -1}, {1, 1, -1, 1, -1, 1, 1, 1}, {1, 1, -1, 1, 1, -1,
      -1, -1}, {1, 1, -1, 1, 1, -1, -1, 1}, {1, 1, -1, 1, 1, -1, 1, -1}, {1, 1,
      -1, 1, 1, -1, 1, 1}, {1, 1, -1, 1, 1, 1, -1, -1}, {1, 1, -1, 1, 1, 1, -1,
      1}, {1, 1, -1, 1, 1, 1, 1, -1}, {1, 1, -1, 1, 1, 1, 1, 1}, {1, 1, 1, -1,
      -1, -1, -1, -1}, {1, 1, 1, -1, -1, -1, -1, 1}, {1, 1, 1, -1, -1, -1, 1,
      -1}, {1, 1, 1, -1, -1, -1, 1, 1}, {1, 1, 1, -1, -1, 1, -1, -1}, {1, 1, 1,
      -1, -1, 1, -1, 1}, {1, 1, 1, -1, -1, 1, 1, -1}, {1, 1, 1, -1, -1, 1, 1,
      1}, {1, 1, 1, -1, 1, -1, -1, -1}, {1, 1, 1, -1, 1, -1, -1, 1}, {1, 1, 1,
      -1, 1, -1, 1, -1}, {1, 1, 1, -1, 1, -1, 1, 1}, {1, 1, 1, -1, 1, 1, -1,
      -1}, {1, 1, 1, -1, 1, 1, -1, 1}, {1, 1, 1, -1, 1, 1, 1, -1}, {1, 1, 1,
      -1, 1, 1, 1, 1}, {1, 1, 1, 1, -1, -1, -1, -1}, {1, 1, 1, 1, -1, -1, -1,
      1}, {1, 1, 1, 1, -1, -1, 1, -1}, {1, 1, 1, 1, -1, -1, 1, 1}, {1, 1, 1, 1,
      -1, 1, -1, -1}, {1, 1, 1, 1, -1, 1, -1, 1}, {1, 1, 1, 1, -1, 1, 1, -1},
      {1, 1, 1, 1, -1, 1, 1, 1}, {1, 1, 1, 1, 1, -1, -1, -1}, {1, 1, 1, 1, 1,
      -1, -1, 1}, {1, 1, 1, 1, 1, -1, 1, -1}, {1, 1, 1, 1, 1, -1, 1, 1}, {1, 1,
      1, 1, 1, 1, -1, -1}, {1, 1, 1, 1, 1, 1, -1, 1}, {1, 1, 1, 1, 1, 1, 1,
      -1}, {1, 1, 1, 1, 1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {96, 96}; 

  ntry = ntry + 1; 

  // Reset the matrix elements
  for(int i = 0; i < nprocesses; i++ )
  {
    matrix_element[i] = 0.; 
  }
  // Define permutation
  int perm[nexternal]; 
  for(int i = 0; i < nexternal; i++ )
  {
    perm[i] = i; 
  }

  if (sum_hel == 0 || ntry < 10)
  {
    // Calculate the matrix element for all helicities
    for(int ihel = 0; ihel < ncomb; ihel++ )
    {
      if (goodhel[ihel] || ntry < 2)
      {
        calculate_wavefunctions(perm, helicities[ihel]); 
        t[0] = matrix_gdx_ttxepvegcx(); 
        // Mirror initial state momenta for mirror process
        perm[0] = 1; 
        perm[1] = 0; 
        // Calculate wavefunctions
        calculate_wavefunctions(perm, helicities[ihel]); 
        // Mirror back
        perm[0] = 0; 
        perm[1] = 1; 
        // Calculate matrix elements
        t[1] = matrix_gdx_ttxepvegcx(); 
        double tsum = 0; 
        for(int iproc = 0; iproc < nprocesses; iproc++ )
        {
          matrix_element[iproc] += t[iproc]; 
          tsum += t[iproc]; 
        }
        // Store which helicities give non-zero result
        if (tsum != 0. && !goodhel[ihel])
        {
          goodhel[ihel] = true; 
          ngood++; 
          igood[ngood] = ihel; 
        }
      }
    }
    jhel = 0; 
    sum_hel = min(sum_hel, ngood); 
  }
  else
  {
    // Only use the "good" helicities
    for(int j = 0; j < sum_hel; j++ )
    {
      jhel++; 
      if (jhel >= ngood)
        jhel = 0; 
      double hwgt = double(ngood)/double(sum_hel); 
      int ihel = igood[jhel]; 
      calculate_wavefunctions(perm, helicities[ihel]); 
      t[0] = matrix_gdx_ttxepvegcx(); 
      // Mirror initial state momenta for mirror process
      perm[0] = 1; 
      perm[1] = 0; 
      // Calculate wavefunctions
      calculate_wavefunctions(perm, helicities[ihel]); 
      // Mirror back
      perm[0] = 0; 
      perm[1] = 1; 
      // Calculate matrix elements
      t[1] = matrix_gdx_ttxepvegcx(); 
      for(int iproc = 0; iproc < nprocesses; iproc++ )
      {
        matrix_element[iproc] += t[iproc] * hwgt; 
      }
    }
  }

  for (int i = 0; i < nprocesses; i++ )
    matrix_element[i] /= denominators[i]; 



}

//--------------------------------------------------------------------------
// Evaluate |M|^2, including incoming flavour dependence.

double CPPProcess_P0_Sigma_sm_ckm_gdx_ttxepvegcx::sigmaHat() 
{
  // Select between the different processes
  if(id1 == 21 && id2 == -1)
  {
    // Add matrix elements for processes with beams (21, -1)
    return matrix_element[0] * 2; 
  }
  else if(id1 == -1 && id2 == 21)
  {
    // Add matrix elements for processes with beams (-1, 21)
    return matrix_element[1] * 2; 
  }
  else
  {
    // Return 0 if not correct initial state assignment
    return 0.; 
  }
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void CPPProcess_P0_Sigma_sm_ckm_gdx_ttxepvegcx::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
  int i, j; 

  // Calculate all wavefunctions
  vxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]); 
  oxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
  oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]); 
  ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]); 
  ixxxxx(p[perm[4]], mME[4], hel[4], -1, w[4]); 
  oxxxxx(p[perm[5]], mME[5], hel[5], +1, w[5]); 
  vxxxxx(p[perm[6]], mME[6], hel[6], +1, w[6]); 
  ixxxxx(p[perm[7]], mME[7], hel[7], -1, w[7]); 
  VVV1P0_1(w[0], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[8]); 
  FFV1P0_3(w[3], w[2], pars->GC_11, pars->ZERO, pars->ZERO, w[9]); 
  FFV2_3(w[4], w[5], pars->GC_108, pars->mdl_MW, pars->mdl_WW, w[10]); 
  VVV1P0_1(w[8], w[9], pars->GC_10, pars->ZERO, pars->ZERO, w[11]); 
  FFV2_1(w[1], w[10], pars->GC_44, pars->ZERO, pars->ZERO, w[12]); 
  FFV2_2(w[7], w[10], pars->GC_44, pars->ZERO, pars->ZERO, w[13]); 
  FFV1_1(w[1], w[8], pars->GC_11, pars->ZERO, pars->ZERO, w[14]); 
  FFV1_2(w[7], w[9], pars->GC_11, pars->ZERO, pars->ZERO, w[15]); 
  FFV1_2(w[7], w[8], pars->GC_11, pars->ZERO, pars->ZERO, w[16]); 
  FFV1_1(w[1], w[9], pars->GC_11, pars->ZERO, pars->ZERO, w[17]); 
  FFV1_1(w[2], w[8], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[18]); 
  FFV1P0_3(w[3], w[18], pars->GC_11, pars->ZERO, pars->ZERO, w[19]); 
  FFV1_2(w[3], w[8], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[20]); 
  FFV1P0_3(w[20], w[2], pars->GC_11, pars->ZERO, pars->ZERO, w[21]); 
  FFV1_1(w[2], w[0], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[22]); 
  FFV1_2(w[3], w[6], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[23]); 
  FFV1P0_3(w[23], w[22], pars->GC_11, pars->ZERO, pars->ZERO, w[24]); 
  FFV1_1(w[1], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[25]); 
  FFV1P0_3(w[3], w[22], pars->GC_11, pars->ZERO, pars->ZERO, w[26]); 
  FFV2_1(w[25], w[10], pars->GC_44, pars->ZERO, pars->ZERO, w[27]); 
  FFV1_2(w[7], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[28]); 
  FFV2_2(w[28], w[10], pars->GC_44, pars->ZERO, pars->ZERO, w[29]); 
  FFV1_1(w[22], w[6], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[30]); 
  FFV1P0_3(w[3], w[30], pars->GC_11, pars->ZERO, pars->ZERO, w[31]); 
  VVV1P0_1(w[26], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[32]); 
  FFV1_2(w[7], w[26], pars->GC_11, pars->ZERO, pars->ZERO, w[33]); 
  FFV1_1(w[1], w[26], pars->GC_11, pars->ZERO, pars->ZERO, w[34]); 
  FFV1_2(w[3], w[0], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[35]); 
  FFV1_1(w[2], w[6], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[36]); 
  FFV1P0_3(w[35], w[36], pars->GC_11, pars->ZERO, pars->ZERO, w[37]); 
  FFV1P0_3(w[35], w[2], pars->GC_11, pars->ZERO, pars->ZERO, w[38]); 
  FFV1_2(w[35], w[6], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[39]); 
  FFV1P0_3(w[39], w[2], pars->GC_11, pars->ZERO, pars->ZERO, w[40]); 
  VVV1P0_1(w[38], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[41]); 
  FFV1_2(w[7], w[38], pars->GC_11, pars->ZERO, pars->ZERO, w[42]); 
  FFV1_1(w[1], w[38], pars->GC_11, pars->ZERO, pars->ZERO, w[43]); 
  FFV1_1(w[1], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[44]); 
  FFV2_1(w[44], w[10], pars->GC_44, pars->ZERO, pars->ZERO, w[45]); 
  FFV1P0_3(w[3], w[36], pars->GC_11, pars->ZERO, pars->ZERO, w[46]); 
  FFV1P0_3(w[23], w[2], pars->GC_11, pars->ZERO, pars->ZERO, w[47]); 
  FFV1_1(w[44], w[9], pars->GC_11, pars->ZERO, pars->ZERO, w[48]); 
  FFV1_1(w[44], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[49]); 
  VVV1P0_1(w[6], w[9], pars->GC_10, pars->ZERO, pars->ZERO, w[50]); 
  FFV1_2(w[7], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[51]); 
  FFV2_2(w[51], w[10], pars->GC_44, pars->ZERO, pars->ZERO, w[52]); 
  FFV1_2(w[51], w[9], pars->GC_11, pars->ZERO, pars->ZERO, w[53]); 
  FFV1_2(w[51], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[54]); 
  FFV1_1(w[36], w[0], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[55]); 
  FFV1P0_3(w[3], w[55], pars->GC_11, pars->ZERO, pars->ZERO, w[56]); 
  VVV1P0_1(w[0], w[46], pars->GC_10, pars->ZERO, pars->ZERO, w[57]); 
  FFV1_1(w[12], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[58]); 
  FFV1_2(w[13], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[59]); 
  FFV1_2(w[23], w[0], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[60]); 
  FFV1P0_3(w[60], w[2], pars->GC_11, pars->ZERO, pars->ZERO, w[61]); 
  VVV1P0_1(w[0], w[47], pars->GC_10, pars->ZERO, pars->ZERO, w[62]); 
  FFV1_1(w[25], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[63]); 
  VVV1P0_1(w[0], w[9], pars->GC_10, pars->ZERO, pars->ZERO, w[64]); 
  FFV1_1(w[25], w[9], pars->GC_11, pars->ZERO, pars->ZERO, w[65]); 
  FFV1_2(w[28], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[66]); 
  FFV1_2(w[28], w[9], pars->GC_11, pars->ZERO, pars->ZERO, w[67]); 
  VVV1P0_1(w[64], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[68]); 
  FFV1_2(w[7], w[64], pars->GC_11, pars->ZERO, pars->ZERO, w[69]); 
  FFV1_1(w[1], w[64], pars->GC_11, pars->ZERO, pars->ZERO, w[70]); 
  VVV1P0_1(w[0], w[50], pars->GC_10, pars->ZERO, pars->ZERO, w[71]); 
  FFV1_1(w[17], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[72]); 
  FFV1_2(w[15], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[73]); 
  VVVV1P0_1(w[0], w[6], w[9], pars->GC_12, pars->ZERO, pars->ZERO, w[74]); 
  VVVV3P0_1(w[0], w[6], w[9], pars->GC_12, pars->ZERO, pars->ZERO, w[75]); 
  VVVV4P0_1(w[0], w[6], w[9], pars->GC_12, pars->ZERO, pars->ZERO, w[76]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV1_0(w[7], w[12], w[11], pars->GC_11, amp[0]); 
  FFV1_0(w[13], w[1], w[11], pars->GC_11, amp[1]); 
  FFV2_0(w[15], w[14], w[10], pars->GC_44, amp[2]); 
  FFV1_0(w[13], w[14], w[9], pars->GC_11, amp[3]); 
  FFV2_0(w[16], w[17], w[10], pars->GC_44, amp[4]); 
  FFV1_0(w[16], w[12], w[9], pars->GC_11, amp[5]); 
  FFV1_0(w[13], w[17], w[8], pars->GC_11, amp[6]); 
  FFV1_0(w[15], w[12], w[8], pars->GC_11, amp[7]); 
  FFV1_0(w[7], w[12], w[19], pars->GC_11, amp[8]); 
  FFV1_0(w[13], w[1], w[19], pars->GC_11, amp[9]); 
  FFV1_0(w[7], w[12], w[21], pars->GC_11, amp[10]); 
  FFV1_0(w[13], w[1], w[21], pars->GC_11, amp[11]); 
  FFV1_0(w[7], w[12], w[24], pars->GC_11, amp[12]); 
  FFV1_0(w[13], w[1], w[24], pars->GC_11, amp[13]); 
  FFV1_0(w[7], w[27], w[26], pars->GC_11, amp[14]); 
  FFV1_0(w[13], w[25], w[26], pars->GC_11, amp[15]); 
  FFV1_0(w[29], w[1], w[26], pars->GC_11, amp[16]); 
  FFV1_0(w[28], w[12], w[26], pars->GC_11, amp[17]); 
  FFV1_0(w[7], w[12], w[31], pars->GC_11, amp[18]); 
  FFV1_0(w[13], w[1], w[31], pars->GC_11, amp[19]); 
  FFV1_0(w[7], w[12], w[32], pars->GC_11, amp[20]); 
  FFV1_0(w[33], w[12], w[6], pars->GC_11, amp[21]); 
  FFV1_0(w[13], w[1], w[32], pars->GC_11, amp[22]); 
  FFV1_0(w[13], w[34], w[6], pars->GC_11, amp[23]); 
  FFV1_0(w[7], w[12], w[37], pars->GC_11, amp[24]); 
  FFV1_0(w[13], w[1], w[37], pars->GC_11, amp[25]); 
  FFV1_0(w[7], w[27], w[38], pars->GC_11, amp[26]); 
  FFV1_0(w[13], w[25], w[38], pars->GC_11, amp[27]); 
  FFV1_0(w[29], w[1], w[38], pars->GC_11, amp[28]); 
  FFV1_0(w[28], w[12], w[38], pars->GC_11, amp[29]); 
  FFV1_0(w[7], w[12], w[40], pars->GC_11, amp[30]); 
  FFV1_0(w[13], w[1], w[40], pars->GC_11, amp[31]); 
  FFV1_0(w[7], w[12], w[41], pars->GC_11, amp[32]); 
  FFV1_0(w[42], w[12], w[6], pars->GC_11, amp[33]); 
  FFV1_0(w[13], w[1], w[41], pars->GC_11, amp[34]); 
  FFV1_0(w[13], w[43], w[6], pars->GC_11, amp[35]); 
  FFV1_0(w[7], w[45], w[46], pars->GC_11, amp[36]); 
  FFV1_0(w[13], w[44], w[46], pars->GC_11, amp[37]); 
  FFV1_0(w[7], w[45], w[47], pars->GC_11, amp[38]); 
  FFV1_0(w[13], w[44], w[47], pars->GC_11, amp[39]); 
  FFV2_0(w[28], w[48], w[10], pars->GC_44, amp[40]); 
  FFV1_0(w[28], w[45], w[9], pars->GC_11, amp[41]); 
  FFV2_0(w[15], w[49], w[10], pars->GC_44, amp[42]); 
  FFV1_0(w[13], w[49], w[9], pars->GC_11, amp[43]); 
  FFV1_0(w[13], w[48], w[6], pars->GC_11, amp[44]); 
  FFV1_0(w[7], w[45], w[50], pars->GC_11, amp[45]); 
  FFV1_0(w[15], w[45], w[6], pars->GC_11, amp[46]); 
  FFV1_0(w[13], w[44], w[50], pars->GC_11, amp[47]); 
  FFV1_0(w[52], w[1], w[46], pars->GC_11, amp[48]); 
  FFV1_0(w[51], w[12], w[46], pars->GC_11, amp[49]); 
  FFV1_0(w[52], w[1], w[47], pars->GC_11, amp[50]); 
  FFV1_0(w[51], w[12], w[47], pars->GC_11, amp[51]); 
  FFV2_0(w[53], w[25], w[10], pars->GC_44, amp[52]); 
  FFV1_0(w[52], w[25], w[9], pars->GC_11, amp[53]); 
  FFV2_0(w[54], w[17], w[10], pars->GC_44, amp[54]); 
  FFV1_0(w[54], w[12], w[9], pars->GC_11, amp[55]); 
  FFV1_0(w[53], w[12], w[6], pars->GC_11, amp[56]); 
  FFV1_0(w[52], w[1], w[50], pars->GC_11, amp[57]); 
  FFV1_0(w[52], w[17], w[6], pars->GC_11, amp[58]); 
  FFV1_0(w[51], w[12], w[50], pars->GC_11, amp[59]); 
  FFV1_0(w[7], w[12], w[56], pars->GC_11, amp[60]); 
  FFV1_0(w[13], w[1], w[56], pars->GC_11, amp[61]); 
  FFV1_0(w[7], w[12], w[57], pars->GC_11, amp[62]); 
  FFV1_0(w[7], w[58], w[46], pars->GC_11, amp[63]); 
  FFV1_0(w[13], w[1], w[57], pars->GC_11, amp[64]); 
  FFV1_0(w[59], w[1], w[46], pars->GC_11, amp[65]); 
  FFV1_0(w[7], w[12], w[61], pars->GC_11, amp[66]); 
  FFV1_0(w[13], w[1], w[61], pars->GC_11, amp[67]); 
  FFV1_0(w[7], w[12], w[62], pars->GC_11, amp[68]); 
  FFV1_0(w[7], w[58], w[47], pars->GC_11, amp[69]); 
  FFV1_0(w[13], w[1], w[62], pars->GC_11, amp[70]); 
  FFV1_0(w[59], w[1], w[47], pars->GC_11, amp[71]); 
  FFV2_0(w[15], w[63], w[10], pars->GC_44, amp[72]); 
  FFV1_0(w[13], w[63], w[9], pars->GC_11, amp[73]); 
  FFV1_0(w[7], w[27], w[64], pars->GC_11, amp[74]); 
  FFV1_0(w[13], w[25], w[64], pars->GC_11, amp[75]); 
  FFV1_0(w[13], w[65], w[0], pars->GC_11, amp[76]); 
  FFV1_0(w[15], w[27], w[0], pars->GC_11, amp[77]); 
  FFV2_0(w[66], w[17], w[10], pars->GC_44, amp[78]); 
  FFV1_0(w[66], w[12], w[9], pars->GC_11, amp[79]); 
  FFV1_0(w[29], w[1], w[64], pars->GC_11, amp[80]); 
  FFV1_0(w[28], w[12], w[64], pars->GC_11, amp[81]); 
  FFV1_0(w[67], w[12], w[0], pars->GC_11, amp[82]); 
  FFV1_0(w[29], w[17], w[0], pars->GC_11, amp[83]); 
  FFV1_0(w[7], w[12], w[68], pars->GC_11, amp[84]); 
  FFV1_0(w[69], w[12], w[6], pars->GC_11, amp[85]); 
  FFV1_0(w[13], w[1], w[68], pars->GC_11, amp[86]); 
  FFV1_0(w[13], w[70], w[6], pars->GC_11, amp[87]); 
  FFV1_0(w[7], w[12], w[71], pars->GC_11, amp[88]); 
  FFV1_0(w[7], w[58], w[50], pars->GC_11, amp[89]); 
  FFV1_0(w[13], w[1], w[71], pars->GC_11, amp[90]); 
  FFV1_0(w[59], w[1], w[50], pars->GC_11, amp[91]); 
  FFV1_0(w[13], w[72], w[6], pars->GC_11, amp[92]); 
  FFV1_0(w[59], w[17], w[6], pars->GC_11, amp[93]); 
  FFV1_0(w[73], w[12], w[6], pars->GC_11, amp[94]); 
  FFV1_0(w[15], w[58], w[6], pars->GC_11, amp[95]); 
  FFV1_0(w[7], w[12], w[74], pars->GC_11, amp[96]); 
  FFV1_0(w[7], w[12], w[75], pars->GC_11, amp[97]); 
  FFV1_0(w[7], w[12], w[76], pars->GC_11, amp[98]); 
  FFV1_0(w[13], w[1], w[74], pars->GC_11, amp[99]); 
  FFV1_0(w[13], w[1], w[75], pars->GC_11, amp[100]); 
  FFV1_0(w[13], w[1], w[76], pars->GC_11, amp[101]); 

}
double CPPProcess_P0_Sigma_sm_ckm_gdx_ttxepvegcx::matrix_gdx_ttxepvegcx() 
{
  int i, j; 
  // Local variables
  const int ngraphs = 102; 
  const int ncolor = 12; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}; 
  static const double cf[ncolor][ncolor] = {{48, 16, 16, 6, 0, 16, -2, 0, 0,
      16, -2, 0}, {16, 48, 6, 16, 16, 0, 0, -2, 16, 0, 0, -2}, {16, 6, 48, 16,
      -2, 0, 0, 16, -2, 0, 0, 16}, {6, 16, 16, 48, 0, -2, 16, 0, 0, -2, 16, 0},
      {0, 16, -2, 0, 48, 16, 16, 6, 6, -2, -2, -6}, {16, 0, 0, -2, 16, 48, 6,
      16, -2, 6, -6, -2}, {-2, 0, 0, 16, 16, 6, 48, 16, -2, -6, 6, -2}, {0, -2,
      16, 0, 6, 16, 16, 48, -6, -2, -2, 6}, {0, 16, -2, 0, 6, -2, -2, -6, 48,
      16, 16, 6}, {16, 0, 0, -2, -2, 6, -6, -2, 16, 48, 6, 16}, {-2, 0, 0, 16,
      -2, -6, 6, -2, 16, 6, 48, 16}, {0, -2, 16, 0, -6, -2, -2, 6, 6, 16, 16,
      48}};

  // Calculate color flows
  jamp[0] = +1./2. * (+amp[24] + amp[25] + amp[28] + amp[29] -
      std::complex<double> (0, 1) * amp[32] - std::complex<double> (0, 1) *
      amp[34] + amp[35] + amp[36] + amp[37] + amp[40] + amp[41] + amp[44] +
      std::complex<double> (0, 1) * amp[45] + std::complex<double> (0, 1) *
      amp[47] - std::complex<double> (0, 1) * amp[62] + amp[63] -
      std::complex<double> (0, 1) * amp[64] - std::complex<double> (0, 1) *
      amp[80] - std::complex<double> (0, 1) * amp[81] + amp[82] - amp[84] -
      amp[86] - std::complex<double> (0, 1) * amp[87] + amp[88] +
      std::complex<double> (0, 1) * amp[89] + amp[90] + amp[96] + amp[97] +
      amp[99] + amp[100]);
  jamp[1] = +1./2. * (-1./3. * amp[36] - 1./3. * amp[37] - 1./3. * amp[38] -
      1./3. * amp[39] - 1./3. * amp[48] - 1./3. * amp[49] - 1./3. * amp[50] -
      1./3. * amp[51] - 1./3. * amp[63] - 1./3. * amp[65] - 1./3. * amp[69] -
      1./3. * amp[71]);
  jamp[2] = +1./2. * (-1./3. * amp[14] - 1./3. * amp[15] - 1./3. * amp[16] -
      1./3. * amp[17] - 1./3. * amp[21] - 1./3. * amp[23] - 1./3. * amp[26] -
      1./3. * amp[27] - 1./3. * amp[28] - 1./3. * amp[29] - 1./3. * amp[33] -
      1./3. * amp[35]);
  jamp[3] = +1./2. * (+amp[12] + amp[13] + amp[14] + amp[15] +
      std::complex<double> (0, 1) * amp[20] + amp[21] + std::complex<double>
      (0, 1) * amp[22] + amp[50] + amp[51] + amp[52] + amp[53] + amp[56] -
      std::complex<double> (0, 1) * amp[57] - std::complex<double> (0, 1) *
      amp[59] + std::complex<double> (0, 1) * amp[68] + std::complex<double>
      (0, 1) * amp[70] + amp[71] + std::complex<double> (0, 1) * amp[74] +
      std::complex<double> (0, 1) * amp[75] + amp[76] - amp[84] +
      std::complex<double> (0, 1) * amp[85] - amp[86] + amp[88] + amp[90] -
      std::complex<double> (0, 1) * amp[91] + amp[96] + amp[97] + amp[99] +
      amp[100]);
  jamp[4] = +1./2. * (-amp[0] - amp[1] - std::complex<double> (0, 1) * amp[2] -
      std::complex<double> (0, 1) * amp[3] - std::complex<double> (0, 1) *
      amp[7] - std::complex<double> (0, 1) * amp[10] - std::complex<double> (0,
      1) * amp[11] + amp[38] + amp[39] + amp[42] + amp[43] -
      std::complex<double> (0, 1) * amp[45] + amp[46] - std::complex<double>
      (0, 1) * amp[47] + amp[66] + amp[67] - std::complex<double> (0, 1) *
      amp[68] + amp[69] - std::complex<double> (0, 1) * amp[70] - amp[88] -
      std::complex<double> (0, 1) * amp[89] - amp[90] + amp[95] + amp[98] -
      amp[96] + amp[101] - amp[99]);
  jamp[5] = +1./2. * (+1./3. * std::complex<double> (0, 1) * amp[2] + 1./3. *
      std::complex<double> (0, 1) * amp[3] + 1./3. * std::complex<double> (0,
      1) * amp[4] + 1./3. * std::complex<double> (0, 1) * amp[5] + 1./3. *
      std::complex<double> (0, 1) * amp[6] + 1./3. * std::complex<double> (0,
      1) * amp[7] - 1./3. * amp[40] - 1./3. * amp[41] - 1./3. * amp[42] - 1./3.
      * amp[43] - 1./3. * amp[44] - 1./3. * amp[46] - 1./3. * amp[78] - 1./3. *
      amp[79] - 1./3. * amp[82] - 1./3. * amp[83] - 1./3. * amp[92] - 1./3. *
      amp[95]);
  jamp[6] = +1./2. * (+1./3. * std::complex<double> (0, 1) * amp[8] + 1./3. *
      std::complex<double> (0, 1) * amp[9] + 1./3. * std::complex<double> (0,
      1) * amp[10] + 1./3. * std::complex<double> (0, 1) * amp[11] - 1./3. *
      amp[12] - 1./3. * amp[13] - 1./3. * amp[18] - 1./3. * amp[19] - 1./3. *
      amp[66] - 1./3. * amp[67]);
  jamp[7] = +1./2. * (+amp[0] + amp[1] - std::complex<double> (0, 1) * amp[4] -
      std::complex<double> (0, 1) * amp[5] - std::complex<double> (0, 1) *
      amp[6] - std::complex<double> (0, 1) * amp[8] - std::complex<double> (0,
      1) * amp[9] + amp[16] + amp[17] + amp[18] + amp[19] -
      std::complex<double> (0, 1) * amp[20] - std::complex<double> (0, 1) *
      amp[22] + amp[23] + amp[78] + amp[79] + std::complex<double> (0, 1) *
      amp[80] + std::complex<double> (0, 1) * amp[81] + amp[83] + amp[84] +
      amp[86] + std::complex<double> (0, 1) * amp[87] + amp[92] - amp[98] -
      amp[97] - amp[101] - amp[100]);
  jamp[8] = +1./2. * (-amp[0] - amp[1] + std::complex<double> (0, 1) * amp[4] +
      std::complex<double> (0, 1) * amp[5] + std::complex<double> (0, 1) *
      amp[6] + std::complex<double> (0, 1) * amp[8] + std::complex<double> (0,
      1) * amp[9] + amp[48] + amp[49] + amp[54] + amp[55] +
      std::complex<double> (0, 1) * amp[57] + amp[58] + std::complex<double>
      (0, 1) * amp[59] + amp[60] + amp[61] + std::complex<double> (0, 1) *
      amp[62] + std::complex<double> (0, 1) * amp[64] + amp[65] - amp[88] -
      amp[90] + std::complex<double> (0, 1) * amp[91] + amp[93] + amp[98] -
      amp[96] + amp[101] - amp[99]);
  jamp[9] = +1./2. * (-1./3. * std::complex<double> (0, 1) * amp[8] - 1./3. *
      std::complex<double> (0, 1) * amp[9] - 1./3. * std::complex<double> (0,
      1) * amp[10] - 1./3. * std::complex<double> (0, 1) * amp[11] - 1./3. *
      amp[24] - 1./3. * amp[25] - 1./3. * amp[30] - 1./3. * amp[31] - 1./3. *
      amp[60] - 1./3. * amp[61]);
  jamp[10] = +1./2. * (-1./3. * std::complex<double> (0, 1) * amp[2] - 1./3. *
      std::complex<double> (0, 1) * amp[3] - 1./3. * std::complex<double> (0,
      1) * amp[4] - 1./3. * std::complex<double> (0, 1) * amp[5] - 1./3. *
      std::complex<double> (0, 1) * amp[6] - 1./3. * std::complex<double> (0,
      1) * amp[7] - 1./3. * amp[52] - 1./3. * amp[53] - 1./3. * amp[54] - 1./3.
      * amp[55] - 1./3. * amp[56] - 1./3. * amp[58] - 1./3. * amp[72] - 1./3. *
      amp[73] - 1./3. * amp[76] - 1./3. * amp[77] - 1./3. * amp[93] - 1./3. *
      amp[94]);
  jamp[11] = +1./2. * (+amp[0] + amp[1] + std::complex<double> (0, 1) * amp[2]
      + std::complex<double> (0, 1) * amp[3] + std::complex<double> (0, 1) *
      amp[7] + std::complex<double> (0, 1) * amp[10] + std::complex<double> (0,
      1) * amp[11] + amp[26] + amp[27] + amp[30] + amp[31] +
      std::complex<double> (0, 1) * amp[32] + amp[33] + std::complex<double>
      (0, 1) * amp[34] + amp[72] + amp[73] - std::complex<double> (0, 1) *
      amp[74] - std::complex<double> (0, 1) * amp[75] + amp[77] + amp[84] -
      std::complex<double> (0, 1) * amp[85] + amp[86] + amp[94] - amp[98] -
      amp[97] - amp[101] - amp[100]);

  // Sum and square the color flows to get the matrix element
  double matrix = 0; 
  for(i = 0; i < ncolor; i++ )
  {
    ztemp = 0.; 
    for(j = 0; j < ncolor; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j]; 
    matrix = matrix + real(ztemp * conj(jamp[i]))/denom[i]; 
  }

  // Store the leading color flows for choice of color
  for(i = 0; i < ncolor; i++ )
    jamp2[0][i] += real(jamp[i] * conj(jamp[i])); 

  return matrix; 
}



