
#ifndef SETUPINTEGRATOR_H
#define SETUPINTEGRATOR_H

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

#include <ctime>


struct IntegrationResult {
      double weight;
      double time;
      double err;
      double chi2;
} ;


class HypIntegrator
{
  public:
  HypIntegrator();
  ~HypIntegrator();

  MEPhaseSpace* meIntegrator;

  ROOT::Math::Functor* toIntegrateTTHhyp;
  ROOT::Math::Functor* toIntegrateTTLLhyp;
  ROOT::Math::Functor* toIntegrateTTWhyp;
  ROOT::Math::Functor* toIntegrateTTWJJhyp;

  ROOT::Math::GSLMCIntegrator* ig2;
  ROOT::Math::VegasParameters* param;

  int intPoints;

  void InitializeIntegrator(double , int , int , int, int);
  void SetNCalls(int);
  void ResetCounters();
  void SetupIntegrationHypothesis(int, int, int);
  IntegrationResult DoIntegration(double* , double*);
  void FillErrHist(TH1F**);

  private:  
};


HypIntegrator::HypIntegrator(){

  meIntegrator = new MEPhaseSpace();

  toIntegrateTTHhyp = new ROOT::Math::Functor(meIntegrator, &MEPhaseSpace::Eval, 10);
  toIntegrateTTLLhyp = new ROOT::Math::Functor(meIntegrator, &MEPhaseSpace::Eval, 5);
  toIntegrateTTWhyp = new ROOT::Math::Functor(meIntegrator, &MEPhaseSpace::Eval, 9);
  toIntegrateTTWJJhyp = new ROOT::Math::Functor(meIntegrator, &MEPhaseSpace::Eval, 11);
 
  intPoints = 10000; 
  ig2 = new ROOT::Math::GSLMCIntegrator( ROOT::Math::IntegrationMultiDim::kVEGAS, 1.e-12, 1.e-5, intPoints);
 
  param = new ROOT::Math::VegasParameters( *(ig2->ExtraOptions()) );

}

HypIntegrator::~HypIntegrator(){


}

void HypIntegrator::InitializeIntegrator(double comEnergy, int kGenerator, int kTFChoice, int kTFOption, int nPoints){

  meIntegrator->SetComEnergy(comEnergy);
  meIntegrator->SetGenerator(kGenerator);
  meIntegrator->SetTFChoice(kTFChoice);
  meIntegrator->SetTFOption(kTFOption);
  SetNCalls(nPoints);

 return;
}

void HypIntegrator::SetNCalls(int nPoints)
{
  intPoints = nPoints;
  ROOT::Math::IntegratorMultiDimOptions opts = ig2->Options();
  opts.SetNCalls( intPoints  );
  ig2->SetOptions( opts ) ;
  return;
}

void HypIntegrator::SetupIntegrationHypothesis(int kMode, int stageValue, int nPoints){

  meIntegrator->SetIntegrationMode(kMode);

  ROOT::Math::Functor* FunctorHyp;
  if (kMode==kMEM_TTLL_TopAntitopDecay) FunctorHyp = toIntegrateTTLLhyp;
  if (kMode==kMEM_TTH_TopAntitopHiggsDecay || kMode==kMEM_TTH_TopAntitopHiggsSemiLepDecay) FunctorHyp = toIntegrateTTHhyp;
  if (kMode==kMEM_TTW_TopAntitopDecay) FunctorHyp = toIntegrateTTWhyp;
  if (kMode==kMEM_TTWJJ_TopAntitopDecay) FunctorHyp = toIntegrateTTWJJhyp;

  ig2->SetFunction(*FunctorHyp);

  param->stage      = stageValue;
  ig2->SetParameters(*param);

  SetNCalls(nPoints);
  ResetCounters();

  return;
}

void HypIntegrator::ResetCounters(){

  meIntegrator->iIteration = 0;
  meIntegrator->iCall = 0;
  for (int i=0; i<20; i++) meIntegrator->errorCounter[i]=0;

  return;
}

IntegrationResult HypIntegrator::DoIntegration(double* xL, double* xU)
{

  ResetCounters();
  IntegrationResult res;

  double click1 = std::clock(); 
  
  res.weight = ig2->Integral(xL, xU);
  res.time = ( std::clock() - click1 ) / (double) CLOCKS_PER_SEC;
  res.err = ig2->Error();
  res.chi2 = ig2->ChiSqr();

  return res;
}

void HypIntegrator::FillErrHist(TH1F** h){
 
  cout << "nIteration="<<meIntegrator->iIteration<<endl;

  for (int i=0; i<5; i++) {
    cout << "Err "<<i<<", sum "<<meIntegrator->errorCounter[i]<<endl;
    cout<<"proba "<<((double)meIntegrator->errorCounter[i])/((double)meIntegrator->iIteration)<<endl;
    (*h)->Fill(i+0.5, ((double)meIntegrator->errorCounter[i])/((double)meIntegrator->iIteration));
  }

  return;
}


#endif
