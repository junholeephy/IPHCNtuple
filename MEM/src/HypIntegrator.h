
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
#include "ConfigParser.h"
#include "MultiLepton.h"

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

/*
  ROOT::Math::Functor* toIntegrateTTHhyp;
  ROOT::Math::Functor* toIntegrateTTLLhyp;
  ROOT::Math::Functor* toIntegrateTTWhyp;
  ROOT::Math::Functor* toIntegrateTTWJJhyp;
  ROOT::Math::Functor* toIntegrateTTbarflhyp;
  ROOT::Math::Functor* toIntegrateTTbarslhyp;
  ROOT::Math::Functor* toIntegrateTTHhyp_miss1j;
  ROOT::Math::Functor* toIntegrateTTLLhyp_miss1j;
  ROOT::Math::Functor* toIntegrateTTbarslhyp_miss1j;
  ROOT::Math::Functor* toIntegrateTTWhyp_miss1j;
  ROOT::Math::Functor* toIntegrateTTbarflhyp_miss1j;
  ROOT::Math::Functor* toIntegrateTTHhyp_miss2j;
  ROOT::Math::Functor* toIntegrateTTLLhyp_miss2j;
  ROOT::Math::Functor* toIntegrateTTbarslhyp_miss2j;
*/
  ROOT::Math::Functor** toIntegrate;


  ROOT::Math::GSLMCIntegrator* ig2;
  ROOT::Math::VegasParameters* param;

  int intPoints;

  void InitializeIntegrator(double , int , int , int, int, ConfigParser*);
  void SetNCalls(int);
  void ResetCounters();
  void SetupIntegrationHypothesis(int, int, int, int);
  IntegrationResult DoIntegration(double* , double*);
  void FillErrHist(TH1F**);

  private:  
};


HypIntegrator::HypIntegrator(){

  meIntegrator = new MEPhaseSpace();

/*
  toIntegrateTTHhyp = new ROOT::Math::Functor(meIntegrator, &MEPhaseSpace::Eval, 10);
  toIntegrateTTLLhyp = new ROOT::Math::Functor(meIntegrator, &MEPhaseSpace::Eval, 5);
  toIntegrateTTWhyp = new ROOT::Math::Functor(meIntegrator, &MEPhaseSpace::Eval, 9);
  toIntegrateTTWJJhyp = new ROOT::Math::Functor(meIntegrator, &MEPhaseSpace::Eval, 11);
  toIntegrateTTbarslhyp = new ROOT::Math::Functor(meIntegrator, &MEPhaseSpace::Eval, 5); 
  toIntegrateTTbarflhyp = new ROOT::Math::Functor(meIntegrator, &MEPhaseSpace::Eval, 6);

  toIntegrateTTHhyp_miss1j = new ROOT::Math::Functor(meIntegrator, &MEPhaseSpace::Eval, 12);
  toIntegrateTTWhyp_miss1j = new ROOT::Math::Functor(meIntegrator, &MEPhaseSpace::Eval, 11);
  toIntegrateTTLLhyp_miss1j = new ROOT::Math::Functor(meIntegrator, &MEPhaseSpace::Eval, 7);
  toIntegrateTTbarslhyp_miss1j = new ROOT::Math::Functor(meIntegrator, &MEPhaseSpace::Eval, 7);
  toIntegrateTTbarflhyp_miss1j = new ROOT::Math::Functor(meIntegrator, &MEPhaseSpace::Eval, 8);

  toIntegrateTTHhyp_miss2j = new ROOT::Math::Functor(meIntegrator, &MEPhaseSpace::Eval, 14);
  toIntegrateTTLLhyp_miss2j = new ROOT::Math::Functor(meIntegrator, &MEPhaseSpace::Eval, 9);
  toIntegrateTTbarslhyp_miss2j = new ROOT::Math::Functor(meIntegrator, &MEPhaseSpace::Eval, 9);
*/

  toIntegrate = new ROOT::Math::Functor*[15];
  for (int i=0; i<15; i++){
    toIntegrate[i] = new ROOT::Math::Functor(meIntegrator, &MEPhaseSpace::Eval, i);
  }

  intPoints = 10000; 
  ig2 = new ROOT::Math::GSLMCIntegrator( ROOT::Math::IntegrationMultiDim::kVEGAS, 1.e-12, 1.e-5, intPoints);
 
  param = new ROOT::Math::VegasParameters( *(ig2->ExtraOptions()) );

}

HypIntegrator::~HypIntegrator(){


}

void HypIntegrator::InitializeIntegrator(double comEnergy, int kGenerator, int kTFChoice, int kTFOption, int nPoints, ConfigParser* cfgParser){

  meIntegrator->SetComEnergy(comEnergy);
  meIntegrator->SetGenerator(kGenerator);
  meIntegrator->SetTFChoice(kTFChoice);

  meIntegrator->SetTFOption(cfgParser->valTFOption);
  meIntegrator->SetOptimization(cfgParser->valOptim);
  meIntegrator->SetOptimization(cfgParser->valOptimTopHad, cfgParser->valOptimTopLep, cfgParser->valOptimHiggs, cfgParser->valOptimW);

  meIntegrator->InitializeMadgraphProcesses(cfgParser->valMadgraphDir);
  meIntegrator->LoadTFfromHisto(cfgParser->valTFfile);

  meIntegrator->SetVerbosity(cfgParser->valVerbosity);

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

void HypIntegrator::SetupIntegrationHypothesis(int kMode, int kCat, int stageValue, int nPoints){

  meIntegrator->SetIntegrationMode(kMode);

  int nparam = meIntegrator->GetNumberIntegrationVar(kMode, kCat);

  ROOT::Math::Functor* FunctorHyp = NULL;
  FunctorHyp = toIntegrate[nparam];
/*
  //2b_2j
  if (kCat==kCat_3l_2b_2j && kMode==kMEM_TTLL_TopAntitopDecay) FunctorHyp = toIntegrateTTLLhyp;
  if (kCat==kCat_3l_2b_2j && (kMode==kMEM_TTH_TopAntitopHiggsDecay || kMode==kMEM_TTH_TopAntitopHiggsSemiLepDecay)) FunctorHyp = toIntegrateTTHhyp;
  if (kMode==kMEM_TTW_TopAntitopDecay) FunctorHyp = toIntegrateTTWhyp;
  if (kCat==kCat_3l_2b_2j && kMode==kMEM_TTWJJ_TopAntitopDecay) FunctorHyp = toIntegrateTTWJJhyp;
  if (kMode==kMEM_TTbar_TopAntitopFullyLepDecay) FunctorHyp = toIntegrateTTbarflhyp;
  if (kCat==kCat_3l_2b_2j && kMode==kMEM_TTbar_TopAntitopSemiLepDecay) FunctorHyp = toIntegrateTTbarslhyp;
  //miss1j
  if ((kCat==kCat_3l_1b_2j || kCat==kCat_3l_2b_1j) && kMode==kMEM_TTLL_TopAntitopDecay) FunctorHyp = toIntegrateTTLLhyp_miss1j;
  if ((kCat==kCat_3l_1b_2j || kCat==kCat_3l_2b_1j) && (kMode==kMEM_TTH_TopAntitopHiggsDecay || kMode==kMEM_TTH_TopAntitopHiggsSemiLepDecay)) FunctorHyp = toIntegrateTTHhyp_miss1j;
  if ((kCat==kCat_3l_1b_2j || kCat==kCat_3l_2b_1j) && kMode==kMEM_TTbar_TopAntitopSemiLepDecay) FunctorHyp = toIntegrateTTbarslhyp_miss1j;
  if ((kCat==kCat_3l_1b_2j || kCat==kCat_3l_1b_1j) && kMode==kMEM_TTW_TopAntitopDecay) FunctorHyp = toIntegrateTTWhyp_miss1j;
  if ((kCat==kCat_3l_1b_2j || kCat==kCat_3l_1b_1j) && kMode==kMEM_TTbar_TopAntitopFullyLepDecay) FunctorHyp = toIntegrateTTbarflhyp_miss1j;
  //miss2j
  if ((kCat==kCat_3l_1b_1j || kCat==kCat_3l_2b_0j) && kMode==kMEM_TTLL_TopAntitopDecay) FunctorHyp = toIntegrateTTLLhyp_miss2j;
  if ((kCat==kCat_3l_1b_1j || kCat==kCat_3l_2b_0j) && (kMode==kMEM_TTH_TopAntitopHiggsDecay || kMode==kMEM_TTH_TopAntitopHiggsSemiLepDecay)) FunctorHyp = toIntegrateTTHhyp_miss2j;
  if ((kCat==kCat_3l_1b_1j || kCat==kCat_3l_2b_0j) && kMode==kMEM_TTbar_TopAntitopSemiLepDecay) FunctorHyp = toIntegrateTTbarslhyp_miss2j;
*/

  ig2->SetFunction(*FunctorHyp);

  param->stage      = stageValue;
  ig2->SetParameters(*param);

  SetNCalls(nPoints);
  if (meIntegrator->iNleptons==3){
    if (kMode!=kMEM_TTW_TopAntitopDecay && kMode!=kMEM_TTbar_TopAntitopFullyLepDecay){
      if (kCat==kCat_3l_1b_2j || kCat==kCat_3l_2b_1j) SetNCalls(nPoints*10);
      if (kCat==kCat_3l_1b_1j || kCat==kCat_3l_2b_0j) SetNCalls(nPoints*50);
    }
    if (kMode==kMEM_TTW_TopAntitopDecay || kMode==kMEM_TTbar_TopAntitopFullyLepDecay){
      if (kCat==kCat_3l_1b_2j || kCat==kCat_3l_1b_1j) SetNCalls(nPoints*10);
    }
  }
  else if (meIntegrator->iNleptons==4){
    if (kCat==kCat_4l_2b) SetNCalls(nPoints*3);   
    if (kCat==kCat_4l_1b) SetNCalls(nPoints*30); 
  }
  else if (meIntegrator->iNleptons==2){
    if (kMode!=kMEM_TTW_TopAntitopDecay && kMode!=kMEM_TTbar_TopAntitopSemiLepDecay){
      if (kCat==kCat_2lss_2b_3j || kCat==kCat_2lss_1b_4j) SetNCalls(nPoints*10);
      if (kCat==kCat_2lss_1b_3j || kCat==kCat_2lss_2b_2j) SetNCalls(nPoints*50);
    }
    if (kMode==kMEM_TTW_TopAntitopDecay || kMode==kMEM_TTbar_TopAntitopSemiLepDecay){
      if (kCat==kCat_2lss_1b_4j || kCat==kCat_2lss_1b_3j) SetNCalls(nPoints*10);
    }
  }
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
 
  //cout << "FillErrHist nIteration="<<meIntegrator->iIteration<<endl;

  for (int i=0; i<5; i++) {
    //cout << "Err "<<i<<", sum "<<meIntegrator->errorCounter[i]<<endl;
    //cout<<"proba "<<((double)meIntegrator->errorCounter[i])/((double)meIntegrator->iIteration)<<endl;
    (*h)->Fill(i+0.5, ((double)meIntegrator->errorCounter[i])/((double)meIntegrator->iIteration));
  }

  return;
}
/*
int HypIntegrator::GetNumberIntegrationVar(int kMode, int kCatJet){

  int nparam=5;
  if (kMode==kMEM_TTLL_TopAntitopDecay) nparam = 5;
  if (kMode==kMEM_TTH_TopAntitopHiggsDecay || kMode==kMEM_TTH_TopAntitopHiggsSemiLepDecay) nparam = 10;
  if (kMode==kMEM_TTW_TopAntitopDecay) nparam = 9;
  if (kMode==kMEM_TTWJJ_TopAntitopDecay) nparam = 11;
  if (kMode==kMEM_TTbar_TopAntitopSemiLepDecay) nparam = 5;
  if (kMode==kMEM_TTbar_TopAntitopFullyLepDecay) nparam = 6;

  if (kMode!=kMEM_TTW_TopAntitopDecay && kMode!=kMEM_TTbar_TopAntitopFullyLepDecay){
    if (kCatJets==kCat_3l_2b_1j || kCatJets==kCat_3l_1b_2j) nparam += 2;
    else if (kCatJets==kCat_3l_2b_0j || kCatJets==kCat_3l_1b_1j) nparam += 4;
  }
  if (kMode==kMEM_TTW_TopAntitopDecay || kMode==kMEM_TTbar_TopAntitopFullyLepDecay){
    if (kCatJets==kCat_3l_1b_1j || kCatJets==kCat_3l_1b_2j) nparam += 2;
  }

  return nparam;
}
*/
#endif
