#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TInterpreter.h"

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
#include "ReadGenFlatTree.h"
#include "ConfigParser.h"

#include <iostream>
#include <ctime>

using namespace std;

const float mZ = 91.18;
const float mW = 80.38;

int main(int argc, char *argv[])
{

  gInterpreter->EnableAutoLoading();

  string InputFileName = string(argv[1]);
  ReadGenFlatTree tree;

  int evMax = atoi(argv[3]);
  cout << "Will run on "<<argv[1]<<" from event "<<argv[2]<<" to "<<evMax<<" with option "<< argv[4]<<endl;

  int nhyp;
  string* shyp;
  int* hyp;
  int* nPointsHyp;

  ConfigParser cfgParser;
  //cfgParser.GetConfigFromFile("/afs/cern.ch/work/c/chanon/MEM/src/config.cfg");
  cfgParser.GetConfigFromFile(string(argv[5]));
  cfgParser.LoadHypotheses(&nhyp, &shyp, &hyp, &nPointsHyp);

  MultiLepton multiLepton;
  cfgParser.LoadIntegrationRange(&multiLepton.JetTFfracmin, &multiLepton.JetTFfracmax, &multiLepton.NeutMaxE);


  if (strcmp(argv[4],"--dryRun")==0){

    tree.InitializeDryRun(InputFileName);
    if (atoi(argv[3])==-1 || evMax>tree.tInput->GetEntries()) evMax = tree.tInput->GetEntries();

    cout << "Start running on events with --dryRun option"<<endl;
    for (Long64_t iEvent=atoi(argv[2]); iEvent<evMax; iEvent++){

      if (iEvent % 1000 == 0) cout << "Event "<<iEvent<<endl;

      tree.mc_event = iEvent;
      tree.FillGenMultilepton(iEvent, &multiLepton);

      bool sel = tree.ApplyGenSelection(iEvent, &multiLepton);
      if (sel==0) continue;

      tree.WriteMultilepton(&multiLepton);
   
      tree.tOutput->Fill();
    }

    tree.tOutput->Write();
    tree.fOutput->Close();
  }
  else if (strcmp(argv[4],"--MEMRun")==0){ 
 
    cout << "Test MEM weight computation" <<endl;

    int index[7];
    for (int ih=0; ih<nhyp; ih++){
      if (shyp[ih]=="TTLL") index[ih] = 0;
      if (shyp[ih]=="TTHfl") index[ih] = 1;
      if (shyp[ih]=="TTHsl") index[ih] = 2;
      if (shyp[ih]=="TTW") index[ih] = 3;
      if (shyp[ih]=="TTWJJ") index[ih] = 4;
      if (shyp[ih]=="TTbarfl") index[ih] = 5;
      if (shyp[ih]=="TTbarsl") index[ih] = 6;
    }
 
    int doOptim;
    cfgParser.LoadOptim(&doOptim);
    int doOptimTopHad, doOptimTopLep, doOptimHiggs, doOptimW;
    cfgParser.LoadOptim(&doOptimTopHad, &doOptimTopLep, &doOptimHiggs, &doOptimW);

    HypIntegrator hypIntegrator;
    //hypIntegrator.InitializeIntegrator(13000, kMadgraph, kTFGaussian, kTFRecoilPtot, 10000, kOptimizeNone);
    //hypIntegrator.InitializeIntegrator(13000, kMadgraph, kTFHistoBnonB_GausRecoil, kTFRecoilmET, 10000, kOptimizeNone);
    //hypIntegrator.InitializeIntegrator(13000, kMadgraph, kTFHistoBnonBmET, kTFRecoilmET, 10000, kOptimizeNone);
    //hypIntegrator.InitializeIntegrator(13000, kMadgraph, kTFHistoBnonBmET, kTFRecoilmET, 10000, kOptimizeMw);
    //hypIntegrator.InitializeIntegrator(13000, kMadgraph, kTFHistoBnonBmET, kTFRecoilmET, 10000, doOptim);
    //hypIntegrator.InitializeIntegrator(13000, kMadgraph, kTFHistoBnonBmET, kTFRecoilmET, 10000, &cfgParser);
    hypIntegrator.InitializeIntegrator(&cfgParser);

    int doMinimization = cfgParser.valDoMinimization;

    //hypIntegrator.meIntegrator->SetVerbosity(2);
    //hypIntegrator.meIntegrator->SetVerbosity(1);

    string JetChoice;
    cfgParser.LoadJetChoice(&JetChoice);

    double xsTTH = hypIntegrator.meIntegrator->xsTTH * hypIntegrator.meIntegrator->brTopHad * hypIntegrator.meIntegrator->brTopLep * hypIntegrator.meIntegrator->brHiggsFullLep;
    double xsTTLL = hypIntegrator.meIntegrator->xsTTLL * hypIntegrator.meIntegrator->brTopHad * hypIntegrator.meIntegrator->brTopLep;
    double xsTTW = hypIntegrator.meIntegrator->xsTTW * hypIntegrator.meIntegrator->brTopLep * hypIntegrator.meIntegrator->brTopLep;
    double xsTTbar = hypIntegrator.meIntegrator->xsTTbar * hypIntegrator.meIntegrator->brTopHad * hypIntegrator.meIntegrator->brTopLep;

    bool doPermutationLep = true;
    bool doPermutationJet = true;
    bool doPermutationBjet = true;

    IntegrationResult res;

    tree.InitializeMEMRun(InputFileName);
    if (atoi(argv[3])==-1 || evMax>tree.tInput->GetEntries()) evMax = tree.tInput->GetEntries();

  int nEvents = evMax - atoi(argv[2]);
  int nErrHist = nEvents*12*nhyp;
  tree.hMEPhaseSpace_Error = new TH1F*[nErrHist];
  tree.hMEPhaseSpace_ErrorTot = new TH1F*[nhyp];
  tree.hMEPhaseSpace_ErrorTot_Pass = new TH1F*[nhyp];
  tree.hMEPhaseSpace_ErrorTot_Fail = new TH1F*[nhyp];
  string snum1[12];
  stringstream ss1[12];
  for (int i=0; i<12; i++) {
    ss1[i] << i;
    ss1[i] >> snum1[i];
  }
  int i0=0;
  string snum2[nEvents];
  stringstream ss2[nEvents];
  for (int i=atoi(argv[2]); i<evMax; i++) {
    ss2[i0] << i0;
    ss2[i0] >> snum2[i0];
    i0++;
  }
  string sname="";
  int it=0, ie0=0;
  for (int ie=atoi(argv[2]); ie<evMax; ie++){
    for (int ip=0; ip<12; ip++){
      for (int ih=0; ih<nhyp; ih++){
        sname = "hErrHist_ev" + snum2[ie0] + "_perm" + snum1[ip] + "_hyp" + shyp[ih];
        tree.hMEPhaseSpace_Error[it] = new TH1F(sname.c_str(), sname.c_str(), 10, 0, 10);
        cout << "Creating error histo, event "<<ie<<", permutation "<<ip<<", hyp "<<ih<<" "<<tree.hMEPhaseSpace_Error[it]->GetName() <<endl;
        it++;
      }
    }
    ie0++;
  }
  for (int ih=0; ih<nhyp; ih++){
    sname = "hErrHist_Tot_hyp" + shyp[ih];
    tree.hMEPhaseSpace_ErrorTot[ih] = new TH1F(sname.c_str(), sname.c_str(), 10, 0, 10);
    sname = "hErrHist_Tot_hyp" + shyp[ih] + "_Pass";
    tree.hMEPhaseSpace_ErrorTot_Pass[ih] = new TH1F(sname.c_str(), sname.c_str(), 10, 0, 10);
    sname = "hErrHist_Tot_hyp" + shyp[ih] + "_Fail";
    tree.hMEPhaseSpace_ErrorTot_Fail[ih] = new TH1F(sname.c_str(), sname.c_str(), 10, 0, 10);
  }
  
  string index_CatJets[12];
  index_CatJets[0] = "3l_2b_2j";
  index_CatJets[1] = "3l_1b_2j";
  index_CatJets[2] = "3l_2b_1j";
  index_CatJets[3] = "3l_1b_1j";
  index_CatJets[4] = "3l_2b_0j";
  index_CatJets[5] = "4l_2b";
  index_CatJets[6] = "4l_1b";
  index_CatJets[7] = "2lss_2b_4j";
  index_CatJets[8] = "2lss_1b_4j";
  index_CatJets[9] = "2lss_2b_3j";
  index_CatJets[10] = "2lss_1b_3j";
  index_CatJets[11] = "2lss_2b_2j";

  string index_Hyp[7]; 
  index_Hyp[0] = "TTLL";
  index_Hyp[1] = "TTHfl";
  index_Hyp[2] = "TTHsl";
  index_Hyp[3] = "TTW";
  index_Hyp[4] = "TTWJJ";
  index_Hyp[5] = "TTbarfl";
  index_Hyp[6] = "TTbarsl";

  Double_t index_XS[7];
  index_XS[0] = xsTTLL;
  index_XS[1] = xsTTH;
  index_XS[2] = xsTTH;
  index_XS[3] = xsTTW;
  index_XS[4] = xsTTW;
  index_XS[5] = xsTTbar;
  index_XS[6] = xsTTbar;

  bool sel_tmp=0;
  
  for(int index_hyp=0; index_hyp<7; index_hyp++){
     for(int index_cat=0; index_cat<12; index_cat++){
   	  sname = "hMEAllWeights_" + index_Hyp[index_hyp] + "_" + index_CatJets[index_cat];
   	  tree.hMEAllWeights[index_hyp][index_cat] = new TH1D(sname.c_str(),sname.c_str(),1000,0,1e-9);
   	  sname = "hMEAllWeights_nlog_" + index_Hyp[index_hyp] + "_" + index_CatJets[index_cat];
   	  tree.hMEAllWeights_nlog[index_hyp][index_cat] = new TH1D(sname.c_str(),sname.c_str(),700,0,700);
     }
  }
                                                                                                                     
    int ievent=0;
    for (Long64_t iEvent=atoi(argv[2]); iEvent<evMax; iEvent++){

      cout << "iEvent="<<iEvent<<endl;

      tree.ReadMultilepton(iEvent, &multiLepton);

      tree.multilepton_Bjet1_P4 = *tree.multilepton_Bjet1_P4_ptr;
      tree.multilepton_Bjet2_P4 = *tree.multilepton_Bjet2_P4_ptr;
      tree.multilepton_Lepton1_P4 = *tree.multilepton_Lepton1_P4_ptr;
      tree.multilepton_Lepton2_P4 = *tree.multilepton_Lepton2_P4_ptr;
      tree.multilepton_Lepton3_P4 = *tree.multilepton_Lepton3_P4_ptr;
      tree.multilepton_Lepton4_P4 = *tree.multilepton_Lepton4_P4_ptr;
      tree.multilepton_JetHighestPt1_P4 = *tree.multilepton_JetHighestPt1_P4_ptr;
      tree.multilepton_JetHighestPt2_P4 = *tree.multilepton_JetHighestPt2_P4_ptr;
      tree.multilepton_JetClosestMw1_P4 = *tree.multilepton_JetClosestMw1_P4_ptr;
      tree.multilepton_JetClosestMw2_P4 = *tree.multilepton_JetClosestMw2_P4_ptr;
      tree.multilepton_JetLowestMjj1_P4 = *tree.multilepton_JetLowestMjj1_P4_ptr;
      tree.multilepton_JetLowestMjj2_P4 = *tree.multilepton_JetLowestMjj2_P4_ptr;
      tree.multilepton_JetHighestPt1_2ndPair_P4 = *tree.multilepton_JetHighestPt1_2ndPair_P4_ptr;
      tree.multilepton_JetHighestPt2_2ndPair_P4 = *tree.multilepton_JetHighestPt2_2ndPair_P4_ptr;
      tree.multilepton_JetClosestMw1_2ndPair_P4 = *tree.multilepton_JetClosestMw1_2ndPair_P4_ptr;
      tree.multilepton_JetClosestMw2_2ndPair_P4 = *tree.multilepton_JetClosestMw2_2ndPair_P4_ptr;
      tree.multilepton_JetLowestMjj1_2ndPair_P4 = *tree.multilepton_JetLowestMjj1_2ndPair_P4_ptr;
      tree.multilepton_JetLowestMjj2_2ndPair_P4 = *tree.multilepton_JetLowestMjj2_2ndPair_P4_ptr;

      tree.multilepton_mET = *tree.multilepton_mET_ptr;
      tree.multilepton_Ptot = *tree.multilepton_Ptot_ptr;

      int ncombcheck = 0;
      int combcheck = 1;
      int permreslep = 1;
      int permresjet = 1;
      int permresbjet = 1;
      int nHypAllowed[7];
      int nNullResult[7];
      double weight_max[7];
      for (int ih=0; ih<7; ih++) {
        nHypAllowed[ih] = 0;
        nNullResult[ih] = 0;
        weight_max[ih] = 0;
      }

      tree.mc_mem_ttz_weight_logmean = 0;
      tree.mc_mem_ttz_weight_log = 0;
      tree.mc_mem_ttz_weight = 0;
      tree.mc_mem_ttz_weight_time = 0;
      tree.mc_mem_ttz_weight_err = 0;
      tree.mc_mem_ttz_weight_chi2 = 0;
      tree.mc_mem_ttz_weight_max = 0;
      tree.mc_mem_ttz_weight_avg = 0;

      tree.mc_mem_ttw_weight_logmean = 0;
      tree.mc_mem_ttw_weight_log = 0;
      tree.mc_mem_ttw_weight = 0;
      tree.mc_mem_ttw_weight_time = 0;
      tree.mc_mem_ttw_weight_err = 0;
      tree.mc_mem_ttw_weight_chi2 = 0;
      tree.mc_mem_ttw_weight_max = 0;
      tree.mc_mem_ttw_weight_avg = 0;

      tree.mc_mem_ttwjj_weight_logmean = 0;
      tree.mc_mem_ttwjj_weight_log = 0;
      tree.mc_mem_ttwjj_weight = 0;
      tree.mc_mem_ttwjj_weight_time = 0;
      tree.mc_mem_ttwjj_weight_err = 0;
      tree.mc_mem_ttwjj_weight_chi2 = 0;
      tree.mc_mem_ttwjj_weight_max = 0;
      tree.mc_mem_ttwjj_weight_avg = 0;

      tree.mc_mem_tthfl_weight_logmean = 0;
      tree.mc_mem_tthfl_weight_log = 0;
      tree.mc_mem_tthfl_weight = 0;
      tree.mc_mem_tthfl_weight_time = 0;
      tree.mc_mem_tthfl_weight_err = 0;
      tree.mc_mem_tthfl_weight_chi2 = 0;
      tree.mc_mem_tthfl_weight_max = 0;
      tree.mc_mem_tthfl_weight_avg = 0;   

      tree.mc_mem_tthsl_weight_logmean = 0;
      tree.mc_mem_tthsl_weight_log = 0;
      tree.mc_mem_tthsl_weight = 0;
      tree.mc_mem_tthsl_weight_time = 0;
      tree.mc_mem_tthsl_weight_err = 0;
      tree.mc_mem_tthsl_weight_chi2 = 0;
      tree.mc_mem_tthsl_weight_max = 0;
      tree.mc_mem_tthsl_weight_avg = 0;

      tree.mc_mem_tth_weight_logmean = 0;
      tree.mc_mem_tth_weight_log = 0;
      tree.mc_mem_tth_weight = 0;
      tree.mc_mem_tth_weight_time = 0;
      tree.mc_mem_tth_weight_err = 0;
      tree.mc_mem_tth_weight_chi2 = 0;
      tree.mc_mem_tth_weight_max = 0;
      tree.mc_mem_tth_weight_avg = 0;

      tree.mc_mem_ttbarfl_weight_logmean = 0;
      tree.mc_mem_ttbarfl_weight_log = 0;
      tree.mc_mem_ttbarfl_weight = 0;
      tree.mc_mem_ttbarfl_weight_time = 0;
      tree.mc_mem_ttbarfl_weight_err = 0;
      tree.mc_mem_ttbarfl_weight_chi2 = 0;
      tree.mc_mem_ttbarfl_weight_max = 0;
      tree.mc_mem_ttbarfl_weight_avg = 0;

      tree.mc_mem_ttbarsl_weight_logmean = 0;
      tree.mc_mem_ttbarsl_weight_log = 0;
      tree.mc_mem_ttbarsl_weight = 0;
      tree.mc_mem_ttbarsl_weight_time = 0;
      tree.mc_mem_ttbarsl_weight_err = 0;
      tree.mc_mem_ttbarsl_weight_chi2 = 0;
      tree.mc_mem_ttbarsl_weight_max = 0;
      tree.mc_mem_ttbarsl_weight_avg = 0;

      tree.mc_mem_ttbar_weight_logmean = 0;
      tree.mc_mem_ttbar_weight_log = 0;
      tree.mc_mem_ttbar_weight = 0;
      tree.mc_mem_ttbar_weight_time = 0;
      tree.mc_mem_ttbar_weight_err = 0;
      tree.mc_mem_ttbar_weight_chi2 = 0;
      tree.mc_mem_ttbar_weight_max = 0;
      tree.mc_mem_ttbar_weight_avg = 0;


   //multiLepton.kCatJets = tree.catJets;
   cout << "Nleptons="<< multiLepton.Leptons.size()<<", Category: "<<multiLepton.kCatJets<<endl;

   hypIntegrator.meIntegrator->SetNleptonMode(multiLepton.Leptons.size());

   int iperm=0;
   for (int ih=0; ih<nhyp; ih++){
     hypIntegrator.SetupIntegrationHypothesis(hyp[ih], multiLepton.kCatJets, 0, nPointsHyp[ih]);

     if (doMinimization) hypIntegrator.SetupMinimizerHypothesis(hyp[ih], multiLepton.kCatJets, 0, nPointsHyp[ih]);

     if (multiLepton.Leptons.size()==4 && (hyp[ih]!=kMEM_TTH_TopAntitopHiggsDecay && hyp[ih]!=kMEM_TTLL_TopAntitopDecay && hyp[ih]!=kMEM_TTbar_TopAntitopFullyLepDecay)) continue;
     if (multiLepton.Leptons.size()==2 && (hyp[ih]!=kMEM_TTH_TopAntitopHiggsSemiLepDecay && hyp[ih]!=kMEM_TTW_TopAntitopDecay && hyp[ih]!=kMEM_TTbar_TopAntitopSemiLepDecay)) continue;

     if (multiLepton.kCatJets==kCat_3l_2b_2j || multiLepton.kCatJets==kCat_3l_1b_2j){
       if (JetChoice=="HighestPt") multiLepton.SwitchJetsFromAllJets(kJetPair_HighestPt);
       else if (JetChoice=="MjjLowest") multiLepton.SwitchJetsFromAllJets(kJetPair_MjjLowest);
       else if (JetChoice=="MwClosest") multiLepton.SwitchJetsFromAllJets(kJetPair_MwClosest);
       else if (JetChoice=="Mixed") {
         if (shyp[ih]=="TTHsl") multiLepton.SwitchJetsFromAllJets(kJetPair_MjjLowest);
         else if (shyp[ih]=="TTWJJ") multiLepton.SwitchJetsFromAllJets(kJetPair_HighestPt); 
         else multiLepton.SwitchJetsFromAllJets(kJetPair_MwClosest);
       }
     }
     else if (multiLepton.kCatJets==kCat_3l_2b_1j || multiLepton.kCatJets==kCat_3l_1b_1j) multiLepton.SwitchJetsFromAllJets(kJetSingle);
     else if (multiLepton.kCatJets==kCat_2lss_2b_4j || multiLepton.kCatJets==kCat_2lss_1b_4j){
       if (JetChoice=="Mixed"){
         if (shyp[ih]=="TTHsl") multiLepton.SwitchJetsFromAllJets(kTwoJetPair_MwClosest_2ndMjjLowest);
	 else if (shyp[ih]=="TTW" || shyp[ih]=="TTbarsl") multiLepton.SwitchJetsFromAllJets(kJetPair_MwClosest); 
	 else if (shyp[ih]=="TTWJJ") multiLepton.SwitchJetsFromAllJets(kTwoJetPair_MwClosest_2ndHighestPt);
       }
       else if (JetChoice=="MwClosest") {
	 if (shyp[ih]=="TTW" || shyp[ih]=="TTbarsl") multiLepton.SwitchJetsFromAllJets(kJetPair_MwClosest); 
         else multiLepton.SwitchJetsFromAllJets(kTwoJetPair_MwClosest_2ndMwClosest);
       }
     }
     else if (multiLepton.kCatJets==kCat_2lss_2b_3j || multiLepton.kCatJets==kCat_2lss_1b_3j) {
	 if (shyp[ih]=="TTHsl") multiLepton.SwitchJetsFromAllJets(kThreeJets_MwClosest_HighestPt);
         else if (shyp[ih]=="TTW" || shyp[ih]=="TTbarsl") multiLepton.SwitchJetsFromAllJets(kJetPair_MwClosest);
     }
     else if (multiLepton.kCatJets==kCat_2lss_2b_2j) {
       multiLepton.SwitchJetsFromAllJets(kJetPair_MwClosest);
     }

     multiLepton.DoSort(&multiLepton.Leptons);
     multiLepton.DoSort(&multiLepton.Bjets); 
     multiLepton.DoSort(&multiLepton.Jets);

     ncombcheck = 0;

     if ((multiLepton.Leptons.size()==3 && shyp[ih]=="TTW") || shyp[ih]=="TTbarfl" || multiLepton.kCatJets==kCat_3l_2b_0j || multiLepton.Leptons.size()==4) doPermutationJet=false;
     else doPermutationJet=true;

     do {
       multiLepton.DoSort(&multiLepton.Jets);
       do {
         multiLepton.DoSort(&multiLepton.Leptons);

         do {
             //cout << "Lepton0Pt="<<multiLepton.Leptons.at(0).P4.Pt()<<" Lepton1Pt="<<multiLepton.Leptons.at(1).P4.Pt() << " Lepton2Pt="<<multiLepton.Leptons.at(2).P4.Pt()<<endl;
	     //cout << "Jet0Pt="<<multiLepton.Jets.at(0).P4.Pt() <<" Jet1Pt="<<multiLepton.Jets.at(1).P4.Pt() <<endl;
	     //cout << "Bjet0Pt="<<multiLepton.Bjets.at(0).P4.Pt() <<" Bjet1Pt="<<multiLepton.Bjets.at(1).P4.Pt() <<endl;

             combcheck = multiLepton.CheckPermutationHyp(hyp[ih]);
	     //cout <<"Checking comb with hyp="<<shyp[ih]<<": check="<<combcheck<<" (ncombcheck="<< ncombcheck<<")"<<endl; ncombcheck++;
             if (combcheck) {
	       for (unsigned int il=0; il<multiLepton.Leptons.size(); il++) cout << " Lepton"<< il<<"Id="<<multiLepton.Leptons.at(il).Id; cout<<endl;
               //cout << "Jet0Pt="<<multiLepton.Jets.at(0).P4.Pt() <<" Jet1Pt="<<multiLepton.Jets.at(1).P4.Pt() <<endl;
	       //if (multiLepton.Bjets.size()>=2) cout << "BjetSize="<<multiLepton.Bjets.size() <<" Bjet0Pt="<<multiLepton.Bjets.at(0).P4.Pt() <<" Bjet1Pt="<<multiLepton.Bjets.at(1).P4.Pt() <<endl;
               //else if (multiLepton.Bjets.size()==1) cout << "BjetSize="<<multiLepton.Bjets.size() <<" Bjet0Pt="<<multiLepton.Bjets.at(0).P4.Pt() <<endl;

               multiLepton.FillParticlesHypothesis(hyp[ih], &hypIntegrator.meIntegrator);

	       if (doMinimization){
                 hypIntegrator.meIntegrator->SetMinimization(1);
	         double * var;
	         var = hypIntegrator.FindMinimizationiInitialValues(multiLepton.xL, multiLepton.xU);
	         res = hypIntegrator.DoMinimization(multiLepton.xL, multiLepton.xU, var);
	         cout << "Hyp "<< shyp[ih]<<" Vegas Ncall="<<nPointsHyp[ih] <<" -log(max)=" << res.weight<<" Time(s)="<<res.time<<endl;
	         hypIntegrator.meIntegrator->SetMinimization(0);
	       }
	       else {
                 res = hypIntegrator.DoIntegration(multiLepton.xL, multiLepton.xU); 
                 cout << "Hyp "<< shyp[ih]<<" Vegas Ncall="<<nPointsHyp[ih] <<" Cross section (pb) : " << res.weight<< " +/- "<< res.err<<" chi2/ndof="<< res.chi2<<" Time(s)="<<res.time<<endl;
	       }

	       //cout << "ihyp="<<ih<<" iperm="<<iperm<<" Filling error hist "<<ih+nhyp*(nHypAllowed[index[ih]]+nNullResult[index[ih]])<<endl;
	       //hypIntegrator.FillErrHist(&(tree.hMEPhaseSpace_Error[ih+nhyp*(nHypAllowed[index[ih]]+nNullResult[index[ih]])+12*nhyp*ievent]));
	       hypIntegrator.FillErrHist(&(tree.hMEPhaseSpace_ErrorTot[ih]));

               if (res.weight>0){
                 nHypAllowed[index[ih]]++;
                 hypIntegrator.FillErrHist(&(tree.hMEPhaseSpace_ErrorTot_Pass[ih]));
	       }
	       else {
		 //res.weight = std::numeric_limits< double >::min();
		 res.weight = 0;
		 res.err = 0;
		 res.chi2 = 0;
		 cout << "Setting result weight to almost 0:"<< res.weight << endl;
		 hypIntegrator.FillErrHist(&(tree.hMEPhaseSpace_ErrorTot_Fail[ih]));
		 nNullResult[index[ih]]++;
	       }
               if (res.weight > weight_max[ih]) weight_max[ih] = res.weight;

	       if (shyp[ih]=="TTLL"){ 
  		 tree.mc_mem_ttz_weight += res.weight / xsTTLL;
		 if(res.weight>0)
		    tree.mc_mem_ttz_weight_logmean += log(res.weight / xsTTLL);
  	 	 tree.mc_mem_ttz_weight_time += res.time;
  		 tree.mc_mem_ttz_weight_err += res.err / xsTTLL;
  		 tree.mc_mem_ttz_weight_chi2 += res.chi2;
	       }
               if (shyp[ih]=="TTHfl"){
	         tree.mc_mem_tthfl_weight += res.weight / xsTTH;
	         if(res.weight>0)
	            tree.mc_mem_tthfl_weight_logmean += log(res.weight / xsTTH);
                 tree.mc_mem_tthfl_weight_time += res.time;
                 tree.mc_mem_tthfl_weight_err += res.err / xsTTH;
                 tree.mc_mem_tthfl_weight_chi2 += res.chi2;

                 tree.mc_mem_tth_weight += res.weight / xsTTH;
                 if(res.weight>0)
                    tree.mc_mem_tth_weight_logmean += log(res.weight / xsTTH);
                 tree.mc_mem_tth_weight_time += res.time;
                 tree.mc_mem_tth_weight_err += res.err / xsTTH;
                 tree.mc_mem_tth_weight_chi2 += res.chi2;
               }
               if (shyp[ih]=="TTHsl"){
                 tree.mc_mem_tthsl_weight += res.weight / xsTTH;
                 if(res.weight>0)
                    tree.mc_mem_tthsl_weight_logmean += log(res.weight / xsTTH);
                 tree.mc_mem_tthsl_weight_time += res.time;
                 tree.mc_mem_tthsl_weight_err += res.err / xsTTH;
                 tree.mc_mem_tthsl_weight_chi2 += res.chi2;

                 tree.mc_mem_tth_weight += res.weight / xsTTH;
                 if(res.weight>0)
                    tree.mc_mem_tth_weight_logmean += log(res.weight / xsTTH);
                 tree.mc_mem_tth_weight_time += res.time;
                 tree.mc_mem_tth_weight_err += res.err / xsTTH;
                 tree.mc_mem_tth_weight_chi2 += res.chi2;
               }
               if (shyp[ih]=="TTW"){
                 tree.mc_mem_ttw_weight += res.weight / xsTTW;
                 if(res.weight>0)
                    tree.mc_mem_ttw_weight_logmean += log(res.weight / xsTTW);
                 tree.mc_mem_ttw_weight_time += res.time;
                 tree.mc_mem_ttw_weight_err += res.err / xsTTW;
                 tree.mc_mem_ttw_weight_chi2 += res.chi2;
               }
               if (shyp[ih]=="TTWJJ"){
                 tree.mc_mem_ttwjj_weight += res.weight / xsTTW;
                 if(res.weight>0)
                    tree.mc_mem_ttwjj_weight_logmean += log(res.weight / xsTTW);
                 tree.mc_mem_ttwjj_weight_time += res.time;
                 tree.mc_mem_ttwjj_weight_err += res.err / xsTTW;
                 tree.mc_mem_ttwjj_weight_chi2 += res.chi2;
               }
               if (shyp[ih]=="TTbarfl"){
                 tree.mc_mem_ttbarfl_weight += res.weight / xsTTbar;
                 if(res.weight>0)
                    tree.mc_mem_ttbarfl_weight_logmean += log(res.weight / xsTTbar);
                 tree.mc_mem_ttbarfl_weight_time += res.time;
                 tree.mc_mem_ttbarfl_weight_err += res.err / xsTTbar;
                 tree.mc_mem_ttbarfl_weight_chi2 += res.chi2;

                 tree.mc_mem_ttbar_weight += res.weight / xsTTbar;
                 if(res.weight>0)
                    tree.mc_mem_ttbar_weight_logmean += log(res.weight / xsTTbar);
                 tree.mc_mem_ttbar_weight_time += res.time;
                 tree.mc_mem_ttbar_weight_err += res.err / xsTTbar;
                 tree.mc_mem_ttbar_weight_chi2 += res.chi2;
               }
               if (shyp[ih]=="TTbarsl"){
                 tree.mc_mem_ttbarsl_weight += res.weight / xsTTbar;
                 if(res.weight>0)
                    tree.mc_mem_ttbarsl_weight_logmean += log(res.weight / xsTTbar);
                 tree.mc_mem_ttbarsl_weight_time += res.time;
                 tree.mc_mem_ttbarsl_weight_err += res.err / xsTTbar;
                 tree.mc_mem_ttbarsl_weight_chi2 += res.chi2;

                 tree.mc_mem_ttbar_weight += res.weight / xsTTbar;
                 if(res.weight>0)
                    tree.mc_mem_ttbar_weight_logmean += log(res.weight / xsTTbar);
                 tree.mc_mem_ttbar_weight_time += res.time;
                 tree.mc_mem_ttbar_weight_err += res.err / xsTTbar;
                 tree.mc_mem_ttbar_weight_chi2 += res.chi2;
               }

               for(int index_hyp=0; index_hyp<7; index_hyp++){
                  for(int index_cat=0; index_cat<12; index_cat++){
                     if(index_cat<=6)
                        sel_tmp = (tree.is_3l_TTH_SR)*(tree.catJets==index_cat)*(shyp[ih]==index_Hyp[index_hyp]);
                     if(index_cat>=7 && index_cat<=11)
                        sel_tmp = (tree.is_2lss_TTH_SR)*(tree.catJets==index_cat)*(shyp[ih]==index_Hyp[index_hyp]);
                     if(sel_tmp){
                        if(res.weight>0){
                           tree.hMEAllWeights[index_hyp][index_cat]->Fill(res.weight/index_XS[index_hyp],tree.weight);
                           tree.hMEAllWeights_nlog[index_hyp][index_cat]->Fill(-log(res.weight/index_XS[index_hyp]),tree.weight);
                        }
                        else{
                           tree.hMEAllWeights[index_hyp][index_cat]->Fill(1e-300,tree.weight);
                           tree.hMEAllWeights_nlog[index_hyp][index_cat]->Fill(-log(1e-300),tree.weight);
                        }
                     }
                  }
               }

             } //end combcheck
    //cout << "Lepton0Id="<<multiLepton.Leptons.at(0).Id<<" Lepton1Id="<<multiLepton.Leptons.at(1).Id << " Lepton2Id="<<multiLepton.Leptons.at(2).Id<<endl;
    //cout << "Lepton0Pt="<<multiLepton.Leptons.at(0).P4.Pt()<<" Lepton1Pt="<<multiLepton.Leptons.at(1).P4.Pt() << " Lepton2Pt="<<multiLepton.Leptons.at(2).P4.Pt()<<endl;

	    iperm++;

           if (doPermutationLep) {
             if (shyp[ih]!="TTbarsl") permreslep = multiLepton.DoPermutation(&multiLepton.Leptons);
             else permreslep = multiLepton.DoPermutationLinear(&multiLepton.Leptons);
           }
         } while (doPermutationLep && permreslep);
         if (doPermutationJet) {
           if (multiLepton.kCatJets!=kCat_3l_2b_1j && multiLepton.kCatJets!=kCat_3l_1b_1j && multiLepton.kCatJets!=kCat_2lss_2b_3j && multiLepton.kCatJets!=kCat_2lss_1b_3j && multiLepton.kCatJets!=kCat_2lss_2b_2j) permresjet = multiLepton.DoPermutation(&multiLepton.Jets);
	   else permresjet = multiLepton.DoPermutationMissingJet("jet");
         }
       } while (doPermutationJet && permresjet);
       if (doPermutationBjet) {
	 if (multiLepton.kCatJets!=kCat_3l_1b_2j && multiLepton.kCatJets!=kCat_3l_1b_1j && multiLepton.kCatJets!=kCat_4l_1b && multiLepton.kCatJets!=kCat_2lss_1b_4j && multiLepton.kCatJets!=kCat_2lss_1b_3j) permresbjet = multiLepton.DoPermutation(&multiLepton.Bjets);
	 else permresbjet = multiLepton.DoPermutationMissingJet("bjet");
       }
     } while (doPermutationBjet && permresbjet);

   } //end hypothesis for loop

   for (int ih=0; ih<nhyp; ih++) {
     cout<< "hyp "<< shyp[ih]<<" nPermAllowed="<<nHypAllowed[index[ih]]<<" (permutation or kinematics), nPermForbidden="<<nNullResult[index[ih]]<<endl;
     if (nHypAllowed[ih]+nNullResult[ih]>0)
	tree.hMEPhaseSpace_ErrorTot[ih]->Scale(1./(double)(nHypAllowed[ih]+nNullResult[ih]));
     if (nHypAllowed[ih]>0)
	tree.hMEPhaseSpace_ErrorTot_Pass[ih]->Scale(1./(double)(nHypAllowed[ih]));
     if (nNullResult[ih]>0)
        tree.hMEPhaseSpace_ErrorTot_Fail[ih]->Scale(1./(double)(nNullResult[ih]));
   }

   ievent++;
 

   if (nhyp>=0){
     if (tree.mc_mem_ttz_weight>0){
       tree.mc_mem_ttz_weight_avg = tree.mc_mem_ttz_weight / (nHypAllowed[0]+nNullResult[0]);
       tree.mc_mem_ttz_weight_max = weight_max[0] / xsTTLL;
       tree.mc_mem_ttz_weight /= nHypAllowed[0];
       tree.mc_mem_ttz_weight_logmean /= nHypAllowed[0];
       tree.mc_mem_ttz_weight_log = log(tree.mc_mem_ttz_weight);
       tree.mc_mem_ttz_weight_err /= nHypAllowed[0];
       tree.mc_mem_ttz_weight_chi2 /= nHypAllowed[0];
     }
     else {
       tree.mc_mem_ttz_weight = 1e-300;//std::numeric_limits< double >::min();
       tree.mc_mem_ttz_weight_log = log(tree.mc_mem_ttz_weight);
       tree.mc_mem_ttz_weight_logmean = log(1e-300);
       tree.mc_mem_ttz_weight_max = 1e-300;
     }
     if (tree.mc_mem_tthfl_weight>0){       
       tree.mc_mem_tthfl_weight_avg = tree.mc_mem_tthfl_weight / (nHypAllowed[1]+nNullResult[1]);
       tree.mc_mem_tthfl_weight_max = weight_max[1] / xsTTH;
       tree.mc_mem_tthfl_weight /= nHypAllowed[1];
       tree.mc_mem_tthfl_weight_logmean /= nHypAllowed[1];
       tree.mc_mem_tthfl_weight_log = log(tree.mc_mem_tthfl_weight);
       tree.mc_mem_tthfl_weight_err /= nHypAllowed[1];
       tree.mc_mem_tthfl_weight_chi2 /= nHypAllowed[1];
     }
     else {
       tree.mc_mem_tthfl_weight = 1e-300;//std::numeric_limits< double >::min();
       tree.mc_mem_tthfl_weight_log = log(tree.mc_mem_tthfl_weight);
       tree.mc_mem_tthfl_weight_logmean = log(1e-300);
       tree.mc_mem_tthfl_weight_max = 1e-300;
     }
     if (tree.mc_mem_tthsl_weight>0){
       tree.mc_mem_tthsl_weight_avg = tree.mc_mem_tthsl_weight / (nHypAllowed[2]+nNullResult[2]);
       tree.mc_mem_tthsl_weight_max = weight_max[2] / xsTTH;
       tree.mc_mem_tthsl_weight /= nHypAllowed[2];
       tree.mc_mem_tthsl_weight_logmean /= nHypAllowed[2];
       tree.mc_mem_tthsl_weight_log = log(tree.mc_mem_tthsl_weight);
       tree.mc_mem_tthsl_weight_err /= nHypAllowed[2];
       tree.mc_mem_tthsl_weight_chi2 /= nHypAllowed[2];  
     }
     else {
       tree.mc_mem_tthsl_weight = 1e-300;//std::numeric_limits< double >::min();
       tree.mc_mem_tthsl_weight_log = log(tree.mc_mem_tthsl_weight);
       tree.mc_mem_tthsl_weight_logmean = log(1e-300);
       tree.mc_mem_tthsl_weight_max = 1e-300;
     }
     int nHypAllowed_TTH = nHypAllowed[1]+nHypAllowed[2];
     cout << "nHypAllowed_TTH="<<nHypAllowed_TTH<<endl;
     if (nHypAllowed_TTH>0){
       cout << "Greater than 0"<<endl;
       tree.mc_mem_tth_weight_avg = tree.mc_mem_tth_weight / ((nHypAllowed[1]+nNullResult[1])+nHypAllowed[2]+nNullResult[2]);
       tree.mc_mem_tth_weight_max = ((weight_max[1]>weight_max[2])?weight_max[1]:weight_max[2]) / xsTTH;
       tree.mc_mem_tth_weight /= nHypAllowed_TTH;
       tree.mc_mem_tth_weight_logmean /= nHypAllowed_TTH;
       tree.mc_mem_tth_weight_log = log(tree.mc_mem_tth_weight);
       tree.mc_mem_tth_weight_err /= nHypAllowed_TTH;
       tree.mc_mem_tth_weight_chi2 /= nHypAllowed_TTH;
     }
     else {
       cout << "Equal to 0"<<endl;
       tree.mc_mem_tth_weight = 1e-300;//std::numeric_limits< double >::min();
       tree.mc_mem_tth_weight_log = log(tree.mc_mem_tth_weight);
       tree.mc_mem_tth_weight_logmean = log(1e-300);
      tree.mc_mem_tth_weight_max = 1e-300;
     }
     if (nHypAllowed[3]>0){
       tree.mc_mem_ttw_weight_avg = tree.mc_mem_ttw_weight / (nHypAllowed[3]+nNullResult[3]);
       tree.mc_mem_ttw_weight_max = weight_max[3] / xsTTW;
       tree.mc_mem_ttw_weight /= nHypAllowed[3];
       tree.mc_mem_ttw_weight_logmean /= nHypAllowed[3];
       tree.mc_mem_ttw_weight_log = log(tree.mc_mem_ttw_weight);
       tree.mc_mem_ttw_weight_err /= nHypAllowed[3];
       tree.mc_mem_ttw_weight_chi2 /= nHypAllowed[3];
     }
     else {
       tree.mc_mem_ttw_weight = 1e-300;//std::numeric_limits< double >::min();
       tree.mc_mem_ttw_weight_log = log(tree.mc_mem_ttw_weight);
       tree.mc_mem_ttw_weight_logmean = log(1e-300);
       tree.mc_mem_ttw_weight_max = 1e-300;
     }
     if (nHypAllowed[4]>0){
       tree.mc_mem_ttwjj_weight_avg = tree.mc_mem_ttwjj_weight / (nHypAllowed[4]+nNullResult[4]);
       tree.mc_mem_ttwjj_weight_max = weight_max[4] / xsTTW;
       tree.mc_mem_ttwjj_weight /= nHypAllowed[4];
       tree.mc_mem_ttwjj_weight_logmean /= nHypAllowed[4];
       tree.mc_mem_ttwjj_weight_log = log(tree.mc_mem_ttwjj_weight);
       tree.mc_mem_ttwjj_weight_err /= nHypAllowed[4];
       tree.mc_mem_ttwjj_weight_chi2 /= nHypAllowed[4];
     }
     else {
       tree.mc_mem_ttwjj_weight = 1e-300;//std::numeric_limits< double >::min();
       tree.mc_mem_ttwjj_weight_log = log(tree.mc_mem_ttwjj_weight);
       tree.mc_mem_ttwjj_weight_logmean = log(1e-300);
       tree.mc_mem_ttwjj_weight_max = 1e-300;
     }
     if (tree.mc_mem_ttbarfl_weight>0){
       tree.mc_mem_ttbarfl_weight_avg = tree.mc_mem_ttbarfl_weight / (nHypAllowed[1]+nNullResult[1]);
       tree.mc_mem_ttbarfl_weight_max = weight_max[1] / xsTTbar;
       tree.mc_mem_ttbarfl_weight /= nHypAllowed[1];
       tree.mc_mem_ttbarfl_weight_logmean /= nHypAllowed[1];
       tree.mc_mem_ttbarfl_weight_log = log(tree.mc_mem_ttbarfl_weight);
       tree.mc_mem_ttbarfl_weight_err /= nHypAllowed[1];
       tree.mc_mem_ttbarfl_weight_chi2 /= nHypAllowed[1];
     }
     else {
       tree.mc_mem_ttbarfl_weight = 1e-300;//std::numeric_limits< double >::min();
       tree.mc_mem_ttbarfl_weight_log = log(tree.mc_mem_ttbarfl_weight);
       tree.mc_mem_ttbarfl_weight_logmean = log(1e-300);
       tree.mc_mem_ttbarfl_weight_max = 1e-300;
     }
     if (tree.mc_mem_ttbarsl_weight>0){
       tree.mc_mem_ttbarsl_weight_avg = tree.mc_mem_ttbarsl_weight / (nHypAllowed[2]+nNullResult[2]);
       tree.mc_mem_ttbarsl_weight_max = weight_max[2] / xsTTbar;
       tree.mc_mem_ttbarsl_weight /= nHypAllowed[2];
       tree.mc_mem_ttbarsl_weight_logmean /= nHypAllowed[2];
       tree.mc_mem_ttbarsl_weight_log = log(tree.mc_mem_ttbarsl_weight);
       tree.mc_mem_ttbarsl_weight_err /= nHypAllowed[2];
       tree.mc_mem_ttbarsl_weight_chi2 /= nHypAllowed[2];
     }
     else {
       tree.mc_mem_ttbarsl_weight = 1e-300;//std::numeric_limits< double >::min();
       tree.mc_mem_ttbarsl_weight_log = log(tree.mc_mem_ttbarsl_weight);
       tree.mc_mem_ttbarsl_weight_logmean = log(1e-300);
       tree.mc_mem_ttbarsl_weight_max = 1e-300;
     }
     int nHypAllowed_TTbar = nHypAllowed[5]+nHypAllowed[6];
     cout << "nHypAllowed_TTbar="<<nHypAllowed_TTbar<<endl;
     if (nHypAllowed_TTbar>0){
       cout << "Greater than 0"<<endl;
       tree.mc_mem_ttbar_weight_avg = tree.mc_mem_ttbar_weight / ((nHypAllowed[5]+nNullResult[5])+nHypAllowed[6]+nNullResult[6]);
       tree.mc_mem_ttbar_weight_max = ((weight_max[5]>weight_max[6])?weight_max[5]:weight_max[6]) / xsTTbar;
       tree.mc_mem_ttbar_weight /= nHypAllowed_TTbar;
       tree.mc_mem_ttbar_weight_logmean /= nHypAllowed_TTbar;
       tree.mc_mem_ttbar_weight_log = log(tree.mc_mem_ttbar_weight);
       tree.mc_mem_ttbar_weight_err /= nHypAllowed_TTbar;
       tree.mc_mem_ttbar_weight_chi2 /= nHypAllowed_TTbar;
     }
     else {
       cout << "Equal to 0"<<endl;
       tree.mc_mem_ttbar_weight = 1e-300;//std::numeric_limits< double >::min();
       tree.mc_mem_ttbar_weight_log = log(tree.mc_mem_ttbar_weight);
       tree.mc_mem_ttbar_weight_logmean = log(1e-300);
      tree.mc_mem_ttbar_weight_max = 1e-300;
     }

     //tree.mc_mem_ttz_tthfl_likelihood = xsTTLL*tree.mc_mem_ttz_weight / (xsTTLL*tree.mc_mem_ttz_weight + xsTTH*tree.mc_mem_tthfl_weight);
     //tree.mc_mem_ttz_tthsl_likelihood = xsTTLL*tree.mc_mem_ttz_weight / (xsTTLL*tree.mc_mem_ttz_weight + xsTTH*tree.mc_mem_tthsl_weight);
     //tree.mc_mem_ttw_tthfl_likelihood = xsTTW*tree.mc_mem_ttw_weight / (xsTTW*tree.mc_mem_ttw_weight + xsTTH*tree.mc_mem_tthfl_weight);
     //tree.mc_mem_ttw_tthsl_likelihood = xsTTW*tree.mc_mem_ttw_weight / (xsTTW*tree.mc_mem_ttw_weight + xsTTH*tree.mc_mem_tthsl_weight);

     if (tree.mc_mem_ttz_weight>0 && tree.mc_mem_tth_weight>0){
       tree.mc_mem_ttz_tth_likelihood = xsTTLL*tree.mc_mem_ttz_weight / (xsTTLL*tree.mc_mem_ttz_weight + xsTTH*tree.mc_mem_tth_weight);
       tree.mc_mem_ttz_tth_likelihood_nlog = -log(xsTTLL*tree.mc_mem_ttz_weight / (xsTTLL*tree.mc_mem_ttz_weight + xsTTH*tree.mc_mem_tth_weight));
       tree.mc_mem_ttz_tth_likelihood_max = xsTTLL*tree.mc_mem_ttz_weight_max / (xsTTLL*tree.mc_mem_ttz_weight_max + xsTTH*tree.mc_mem_tth_weight_max);
       tree.mc_mem_ttz_tth_likelihood_avg = xsTTLL*tree.mc_mem_ttz_weight_avg / (xsTTLL*tree.mc_mem_ttz_weight_avg + xsTTH*tree.mc_mem_tth_weight_avg);
     }
     if (tree.mc_mem_ttw_weight>0 && tree.mc_mem_tth_weight>0){
       tree.mc_mem_ttw_tth_likelihood = xsTTW*tree.mc_mem_ttw_weight / (xsTTW*tree.mc_mem_ttw_weight + xsTTH*tree.mc_mem_tth_weight);
       tree.mc_mem_ttw_tth_likelihood_nlog = -log(xsTTW*tree.mc_mem_ttw_weight / (xsTTW*tree.mc_mem_ttw_weight + xsTTH*tree.mc_mem_tth_weight));
       tree.mc_mem_ttw_tth_likelihood_max = xsTTW*tree.mc_mem_ttw_weight_max / (xsTTW*tree.mc_mem_ttw_weight_max + xsTTH*tree.mc_mem_tth_weight_max);
       tree.mc_mem_ttw_tth_likelihood_avg = xsTTW*tree.mc_mem_ttw_weight_avg / (xsTTW*tree.mc_mem_ttw_weight_avg + xsTTH*tree.mc_mem_tth_weight_avg);
     }
     if (tree.mc_mem_ttwjj_weight>0 && tree.mc_mem_ttwjj_weight>0){
       tree.mc_mem_ttwjj_tth_likelihood = xsTTW*tree.mc_mem_ttwjj_weight / (xsTTW*tree.mc_mem_ttwjj_weight + xsTTH*tree.mc_mem_tth_weight);
       tree.mc_mem_ttwjj_tth_likelihood_nlog = -log(xsTTW*tree.mc_mem_ttwjj_weight / (xsTTW*tree.mc_mem_ttwjj_weight + xsTTH*tree.mc_mem_tth_weight));
       tree.mc_mem_ttwjj_tth_likelihood_max = xsTTW*tree.mc_mem_ttwjj_weight_max / (xsTTW*tree.mc_mem_ttwjj_weight_max + xsTTH*tree.mc_mem_tth_weight_max);
       tree.mc_mem_ttwjj_tth_likelihood_avg = xsTTW*tree.mc_mem_ttwjj_weight_avg / (xsTTW*tree.mc_mem_ttwjj_weight_avg + xsTTH*tree.mc_mem_tth_weight_avg);
     }
     if (tree.mc_mem_ttbar_weight>0 && tree.mc_mem_tth_weight>0){
       tree.mc_mem_ttbar_tth_likelihood = xsTTbar*tree.mc_mem_ttbar_weight / (xsTTbar*tree.mc_mem_ttbar_weight + xsTTH*tree.mc_mem_tth_weight);
       tree.mc_mem_ttbar_tth_likelihood_nlog = -log(xsTTbar*tree.mc_mem_ttbar_weight / (xsTTbar*tree.mc_mem_ttbar_weight + xsTTH*tree.mc_mem_tth_weight));
       tree.mc_mem_ttbar_tth_likelihood_max = xsTTbar*tree.mc_mem_ttbar_weight_max / (xsTTbar*tree.mc_mem_ttbar_weight_max + xsTTH*tree.mc_mem_tth_weight_max);
       tree.mc_mem_ttbar_tth_likelihood_avg = xsTTbar*tree.mc_mem_ttbar_weight_avg / (xsTTbar*tree.mc_mem_ttbar_weight_avg + xsTTH*tree.mc_mem_tth_weight_avg);
     }

     if (tree.mc_mem_ttz_weight>0 && tree.mc_mem_ttw_weight>0 && tree.mc_mem_tth_weight>0){
       tree.mc_mem_ttv_tth_likelihood = (xsTTLL*tree.mc_mem_ttz_weight + xsTTW*tree.mc_mem_ttw_weight) / (xsTTLL*tree.mc_mem_ttz_weight + xsTTW*tree.mc_mem_ttw_weight + xsTTH*tree.mc_mem_tth_weight);
       tree.mc_mem_ttv_tth_likelihood_nlog = -log((xsTTLL*tree.mc_mem_ttz_weight + xsTTW*tree.mc_mem_ttw_weight) / (xsTTLL*tree.mc_mem_ttz_weight + xsTTW*tree.mc_mem_ttw_weight + xsTTH*tree.mc_mem_tth_weight));
       tree.mc_mem_ttv_tth_likelihood_max = (xsTTLL*tree.mc_mem_ttz_weight_max + xsTTW*tree.mc_mem_ttw_weight_max) / (xsTTLL*tree.mc_mem_ttz_weight_max + xsTTW*tree.mc_mem_ttw_weight_max + xsTTH*tree.mc_mem_tth_weight_max);    
       tree.mc_mem_ttv_tth_likelihood_avg = (xsTTLL*tree.mc_mem_ttz_weight_avg + xsTTW*tree.mc_mem_ttw_weight_avg) / (xsTTLL*tree.mc_mem_ttz_weight_avg + xsTTW*tree.mc_mem_ttw_weight_avg + xsTTH*tree.mc_mem_tth_weight_avg);
     }
     if (tree.mc_mem_ttz_weight>0 && tree.mc_mem_ttwjj_weight>0 && tree.mc_mem_tth_weight>0){
       tree.mc_mem_ttvjj_tth_likelihood = (xsTTLL*tree.mc_mem_ttz_weight + xsTTW*tree.mc_mem_ttwjj_weight) / (xsTTLL*tree.mc_mem_ttz_weight + xsTTW*tree.mc_mem_ttwjj_weight + xsTTH*tree.mc_mem_tth_weight);
       tree.mc_mem_ttvjj_tth_likelihood_nlog = -log((xsTTLL*tree.mc_mem_ttz_weight + xsTTW*tree.mc_mem_ttwjj_weight) / (xsTTLL*tree.mc_mem_ttz_weight + xsTTW*tree.mc_mem_ttwjj_weight + xsTTH*tree.mc_mem_tth_weight));
       tree.mc_mem_ttvjj_tth_likelihood_max = (xsTTLL*tree.mc_mem_ttz_weight_max + xsTTW*tree.mc_mem_ttwjj_weight_max) / (xsTTLL*tree.mc_mem_ttz_weight_max + xsTTW*tree.mc_mem_ttwjj_weight_max + xsTTH*tree.mc_mem_tth_weight_max);
       tree.mc_mem_ttvjj_tth_likelihood_avg = (xsTTLL*tree.mc_mem_ttz_weight_avg + xsTTW*tree.mc_mem_ttwjj_weight_avg) / (xsTTLL*tree.mc_mem_ttz_weight_avg + xsTTW*tree.mc_mem_ttwjj_weight_avg + xsTTH*tree.mc_mem_tth_weight_avg);
     }

     cout << "PRINT DISCRIMINATION"<<endl;
     if (tree.mc_mem_tth_weight>0){
       cout << "TTH nHypAllowed="<<nHypAllowed_TTH<<endl;
       cout << "TTH weight="<<tree.mc_mem_tth_weight<<"; log(weight)="<<tree.mc_mem_tth_weight_log<<endl;
     }
     if (tree.mc_mem_ttz_weight>0){
       cout << "TTZ nHypAllowed="<<nHypAllowed[0]<<endl;
       cout << "TTZ weight="<<tree.mc_mem_ttz_weight<<" ; log(weight)="<<tree.mc_mem_ttz_weight_log<<endl;
     }
     if (tree.mc_mem_ttw_weight>0){
       cout << "TTW nHypAllowed="<<nHypAllowed[3]<<endl;
       cout << "TTW weight="<<tree.mc_mem_ttw_weight<<" log(weight)="<<tree.mc_mem_ttw_weight_log<<endl;
     }
     if (tree.mc_mem_ttwjj_weight>0){
       cout << "TTWJJ nHypAllowed="<<nHypAllowed[4]<<endl;
       cout << "TTWJJ weight="<<tree.mc_mem_ttwjj_weight<<" log(weight)="<<tree.mc_mem_ttwjj_weight_log<<endl;
     }
     if (tree.mc_mem_ttbar_weight>0){
       cout << "TTbar nHypAllowed="<<nHypAllowed_TTbar<<endl;
       cout << "TTbar weight="<<tree.mc_mem_ttbar_weight<<" ; log(weight)="<<tree.mc_mem_ttbar_weight_log<<endl;
     }
     if (tree.mc_mem_tth_weight>0 && tree.mc_mem_ttz_weight)
       cout << "TTHvsTTZ discrim="<<tree.mc_mem_ttz_tth_likelihood<<"; -log(discrim)="<< tree.mc_mem_ttz_tth_likelihood_nlog<<endl;
     if (tree.mc_mem_ttw_weight>0 && tree.mc_mem_tth_weight>0)
       cout << "TTHvsTTW discrim="<<tree.mc_mem_ttw_tth_likelihood<<"; -log(discrim)="<< tree.mc_mem_ttw_tth_likelihood_nlog<<endl;
     if (tree.mc_mem_ttwjj_tth_likelihood>0 && tree.mc_mem_tth_weight>0)
       cout << "TTHvsTTWJJ discrim="<<tree.mc_mem_ttwjj_tth_likelihood<<"; -log(discrim)="<< tree.mc_mem_ttwjj_tth_likelihood_nlog<<endl;
     if (tree.mc_mem_tth_weight>0 && tree.mc_mem_ttz_weight>0 && tree.mc_mem_ttw_weight>0)
       cout << "TTHvsTTV discrim="<<tree.mc_mem_ttv_tth_likelihood<<"; -log(discrim)="<< tree.mc_mem_ttv_tth_likelihood_nlog<<endl;
     if (tree.mc_mem_tth_weight>0 && tree.mc_mem_ttz_weight>0 && tree.mc_mem_ttwjj_weight>0)
       cout << "TTHvsTTVjj discrim="<<tree.mc_mem_ttvjj_tth_likelihood<<"; -log(discrim)="<< tree.mc_mem_ttvjj_tth_likelihood_nlog<<endl;
     
     cout << "End of event ---"<<endl;
     cout << endl;

 
   }

   tree.tOutput->Fill();

  }

  tree.fOutput->cd();

  for (int ihist=0; ihist<nErrHist; ihist++) tree.hMEPhaseSpace_Error[ihist]->Write();
  for (int ihyp=0; ihyp<nhyp; ihyp++) {
    tree.hMEPhaseSpace_ErrorTot[ihyp]->Write();    
    tree.hMEPhaseSpace_ErrorTot_Pass[ihyp]->Write();
    tree.hMEPhaseSpace_ErrorTot_Fail[ihyp]->Write();
  }
  for(int index_hyp=0; index_hyp<7; index_hyp++){
     for(int index_cat=0; index_cat<12; index_cat++){
        tree.hMEAllWeights[index_hyp][index_cat]->Write();
        tree.hMEAllWeights_nlog[index_hyp][index_cat]->Write();
     }
  }

  tree.tOutput->Write();
  tree.fOutput->Close();

 }

}
