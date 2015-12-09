
#ifndef MULTILEPTON_H
#define MULTILEPTON_H

#include <algorithm>

#define kJetPair_HighestPt 0
#define kJetPair_MwClosest 1
#define kJetPair_MjjLowest 2

struct Particle {
  int Id;
  TLorentzVector P4;
};

class MultiLepton
{
  public:
  MultiLepton();
  ~MultiLepton();

  vector<Particle> Bjets;
  vector<Particle> Leptons;
  vector<Particle> Jets;
  vector<Particle> AllJets;
  vector<Particle> JetsHighestPt;
  vector<Particle> JetsClosestMw;
  vector<Particle> JetsLowestMjj;

  TLorentzVector Ptot;
  TLorentzVector mET;

  double *xL;
  double *xU;
  int nParam;

  double mB;
  double JetTFfracmin;
  double JetTFfracmax;
  double NeutMaxE;

  void FillParticle(string, int, TLorentzVector);
  void FillParticlesHypothesis(int, MEPhaseSpace**);   

  void DoSort(vector<Particle>*);
  int DoPermutation(vector<Particle>*);
  int CheckPermutationHyp(int);

  void SwitchJetsFromAllJets(int);

  void FillTTHFullyLepHyp(MEPhaseSpace**);
  void FillTTLLHyp(MEPhaseSpace**);
  void FillTTHSemiLepHyp(MEPhaseSpace**);
  void FillTTWHyp(MEPhaseSpace**);

  void ReadIntegrationBoundaries(int);

  struct ComparePt{
    bool operator() (const Particle& PA, const Particle& PB) const {
      return PA.P4.Pt()>PB.P4.Pt();
    }
  };

  private:
};

MultiLepton::MultiLepton(){

  nParam = 11;
  xL = new double[nParam];
  xU = new double[nParam];
  mB = 4.7;
  JetTFfracmin = 0.65;//0.65; //next try 0.5 (5 sigma at 100 GeV for a 0.2*Egen gaussian width)
  JetTFfracmax = 2.0;//2.0; //next try 5.0 (4.5 sigma at 100 GeV)
  NeutMaxE = 300.; //try 500 or 1000 ?
}

MultiLepton::~MultiLepton(){

}

void MultiLepton::FillParticle(string Type, int id, TLorentzVector p4){

  Particle p;
  p.Id = id;
  p.P4 = p4;
  if (Type=="lepton") Leptons.push_back(p);
  if (Type=="jet") Jets.push_back(p);
  if (Type=="bjet") Bjets.push_back(p);
  if (Type=="alljet") AllJets.push_back(p);
  if (Type=="jetHighestPt") JetsHighestPt.push_back(p);
  if (Type=="jetClosestMw") JetsClosestMw.push_back(p);
  if (Type=="jetLowestMjj") JetsLowestMjj.push_back(p);

  return;
}

void MultiLepton::DoSort(vector<Particle>* particles)
{
  std::sort(((*particles).begin()), ((*particles).end()), MultiLepton::ComparePt());
  return;
}

int MultiLepton::DoPermutation(vector<Particle>* particles)
{
  int res = std::next_permutation(((*particles).begin()), ((*particles).end()), MultiLepton::ComparePt());
  return res;
}


void MultiLepton::FillParticlesHypothesis(int kMode, MEPhaseSpace** meIntegrator)
{

  if (kMode==kMEM_TTLL_TopAntitopDecay) FillTTLLHyp(meIntegrator);
  if (kMode==kMEM_TTH_TopAntitopHiggsDecay) FillTTHFullyLepHyp(meIntegrator);
  if (kMode==kMEM_TTH_TopAntitopHiggsSemiLepDecay) FillTTHSemiLepHyp(meIntegrator);
  if (kMode==kMEM_TTW_TopAntitopDecay || kMode==kMEM_TTWJJ_TopAntitopDecay) FillTTWHyp(meIntegrator);

  ReadIntegrationBoundaries(kMode);

  return;
}

void MultiLepton::ReadIntegrationBoundaries(int kMode){

  int nparam=5;
  if (kMode==kMEM_TTLL_TopAntitopDecay) nparam = 5;
  if (kMode==kMEM_TTH_TopAntitopHiggsDecay || kMode==kMEM_TTH_TopAntitopHiggsSemiLepDecay) nparam = 10;
  if (kMode==kMEM_TTH_TopAntitopHiggsSemiLepDecay) nparam = 10;
  if (kMode==kMEM_TTW_TopAntitopDecay) nparam = 9;
  if (kMode==kMEM_TTWJJ_TopAntitopDecay) nparam = 11;

  for (int i=0; i<nparam; i++) cout << "Var "<<i<<" xL="<<xL[i]<<" xU="<<xU[i]<<endl;

  return;
}

int MultiLepton::CheckPermutationHyp(int kMode){

  int check = 1;

  if (kMode==kMEM_TTLL_TopAntitopDecay){
    if (Leptons[0].Id != -Leptons[1].Id) check=0; // opposite sign same flavour to make a Z
    if (!(Leptons[0].Id>0 && Leptons[1].Id<0)) check=0; // ME needs l+ l- ordering
  }

  if (kMode==kMEM_TTH_TopAntitopHiggsDecay){
    if (!(Leptons[0].Id>0 && Leptons[1].Id<0)) check=0; //opposite sign, and ME needs l+ l- ordering
  }

  if (kMode==kMEM_TTH_TopAntitopHiggsSemiLepDecay){
    if (!(Leptons[1].Id>0 && Leptons[2].Id<0)) check=0; //opposite sign leptons from tops
  }

  if (kMode==kMEM_TTW_TopAntitopDecay || kMode==kMEM_TTWJJ_TopAntitopDecay){
    if (!(Leptons[1].Id>0 && Leptons[2].Id<0)) check=0; //opposite sign leptons from tops
  }

  return check;
}

void MultiLepton::SwitchJetsFromAllJets(int kMode){

  //kMode==0: highest scalar pt pair
  //kMode==1: pair with mass closest to mW

  Jets.clear();

  if (kMode==0){
    Jets.push_back(JetsHighestPt[0]);
    Jets.push_back(JetsHighestPt[1]);
  }
  if (kMode==1){
    Jets.push_back(JetsClosestMw[0]);
    Jets.push_back(JetsClosestMw[1]);
  }
  if (kMode==2){
    Jets.push_back(JetsLowestMjj[0]);
    Jets.push_back(JetsLowestMjj[1]);
  }

  cout << "Using Jets mode "<<kMode<<endl;

  return;
}

void MultiLepton::FillTTHFullyLepHyp(MEPhaseSpace** meIntegrator)
{
 
  (*meIntegrator)->MEMFix_HiggsFullLep.Lep1_E = Leptons[0].P4.E();
  (*meIntegrator)->MEMFix_HiggsFullLep.Lep1_Theta = Leptons[0].P4.Theta();
  (*meIntegrator)->MEMFix_HiggsFullLep.Lep1_Phi = Leptons[0].P4.Phi();
  (*meIntegrator)->MEMFix_HiggsFullLep.Lep2_E = Leptons[1].P4.E();
  (*meIntegrator)->MEMFix_HiggsFullLep.Lep2_Theta = Leptons[1].P4.Theta();
  (*meIntegrator)->MEMFix_HiggsFullLep.Lep2_Phi = Leptons[1].P4.Phi();
    
  (*meIntegrator)->MEMFix_TopHad.Bjet_Theta = Bjets[0].P4.Theta();
  (*meIntegrator)->MEMFix_TopHad.Bjet_Phi = Bjets[0].P4.Phi();
  (*meIntegrator)->MEMFix_TopHad.Jet1_Theta = Jets[0].P4.Theta();
  (*meIntegrator)->MEMFix_TopHad.Jet1_Phi = Jets[0].P4.Phi();
  (*meIntegrator)->MEMFix_TopHad.Jet2_Theta = Jets[1].P4.Theta();
  (*meIntegrator)->MEMFix_TopHad.Jet2_Phi = Jets[1].P4.Phi();
  (*meIntegrator)->MEMFix_TopHad.TopSign = (Leptons[2].Id>0)?-1:1;

  (*meIntegrator)->MEMFix_TopLep.Bjet_Theta = Bjets[0].P4.Theta();
  (*meIntegrator)->MEMFix_TopLep.Bjet_Phi = Bjets[0].P4.Phi();
  (*meIntegrator)->MEMFix_TopLep.Lep_E = Leptons[2].P4.E();
  (*meIntegrator)->MEMFix_TopLep.Lep_Theta = Leptons[2].P4.Theta();
  (*meIntegrator)->MEMFix_TopLep.Lep_Phi = Leptons[2].P4.Phi();
  (*meIntegrator)->MEMFix_TopLep.TopSign = (Leptons[2].Id>0)?1:-1;  

  (*meIntegrator)->MeasuredVarForTF.Jet1_E = Jets[0].P4.E();
  (*meIntegrator)->MeasuredVarForTF.Jet1_Eta = Jets[0].P4.Eta();
  (*meIntegrator)->MeasuredVarForTF.Jet2_E = Jets[1].P4.E();
  (*meIntegrator)->MeasuredVarForTF.Jet2_Eta = Jets[1].P4.Eta();
  (*meIntegrator)->MeasuredVarForTF.Bjet1_E =  Bjets[0].P4.E();
  (*meIntegrator)->MeasuredVarForTF.Bjet1_Eta =  Bjets[0].P4.Eta();
  (*meIntegrator)->MeasuredVarForTF.Bjet2_E =  Bjets[1].P4.E();
  (*meIntegrator)->MeasuredVarForTF.Bjet2_Eta =  Bjets[1].P4.Eta();
  (*meIntegrator)->MeasuredVarForTF.Recoil_Px = -Ptot.Px();
  (*meIntegrator)->MeasuredVarForTF.Recoil_Py = -Ptot.Py();
  (*meIntegrator)->MeasuredVarForTF.mET_Px = mET.Px();
  (*meIntegrator)->MeasuredVarForTF.mET_Py = mET.Py();

    xL[0] = (Bjets[0].P4.E()*JetTFfracmin<mB)?mB:Bjets[0].P4.E()*JetTFfracmin; //TopHad, Bjet_E
    xL[1] = JetTFfracmin*Jets[0].P4.E(); //TopHad, Jet1_E
    xL[2] = (Bjets[1].P4.E()*JetTFfracmin<mB)?mB:Bjets[1].P4.E()*JetTFfracmin; //TopLep, Bjet_E
    xL[3] = 0; //TopLep, Neut_Theta
    xL[4] = 0; //TopLep, Neut_Phi
    xL[5] = 0; //HiggsFullLep, Neut1_E
    xL[6] = 0; //HiggsFullLep, Neut1_Theta
    xL[7] = 0; //HiggsFullLep, Neut1_Phi
    xL[8] = 0; //HiggsFullLep, Neut2_Theta
    xL[9] = 0; //HiggsFullLep, Neut2_Phi

    xU[0] = JetTFfracmax*Bjets[0].P4.E();//meIntegrator->comEnergy; //TopHad, Bjet_E
    xU[1] = JetTFfracmax*Jets[0].P4.E();//meIntegrator->comEnergy; //TopHad, Jet1_E
    xU[2] = JetTFfracmax*Bjets[1].P4.E();//meIntegrator->comEnergy; //TopLep, Bjet_E
    xU[3] = TMath::Pi(); //TopLep, Neut_Theta
    xU[4] = 2*TMath::Pi(); //TopLep, Neut_Phi
    xU[5] = NeutMaxE;//meIntegrator->comEnergy; //HiggsFullLep, Neut1_E
    xU[6] = TMath::Pi(); //HiggsFullLep, Neut1_Theta
    xU[7] = 2*TMath::Pi(); //HiggsFullLep, Neut1_Phi
    xU[8] = TMath::Pi(); //HiggsFullLep, Neut2_Theta
    xU[9] = 2*TMath::Pi(); //HiggsFullLep, Neut2_Phi


  return;
}

void MultiLepton::FillTTLLHyp(MEPhaseSpace** meIntegrator){
  FillTTHFullyLepHyp(meIntegrator);


    xL[0] = (Bjets[0].P4.E()*JetTFfracmin<mB)?mB:Bjets[0].P4.E()*JetTFfracmin; //TopHad, Bjet_E
    xL[1] = Jets[0].P4.E()*JetTFfracmin; //TopHad, Jet1_E
    xL[2] = (Bjets[1].P4.E()*JetTFfracmin<mB)?mB:Bjets[1].P4.E()*JetTFfracmin; //TopLep, Bjet_E
    xL[3] = 0; //TopLep, Neut_Theta
    xL[4] = 0; //TopLep, Neut_Phi

    xU[0] = JetTFfracmax*Bjets[0].P4.E();//meIntegrator->comEnergy; //TopHad, Bjet_E
    xU[1] = JetTFfracmax*Jets[0].P4.E();//meIntegrator->comEnergy; //TopHad, Jet1_E
    xU[2] = JetTFfracmax*Bjets[1].P4.E();//meIntegrator->comEnergy; //TopLep, Bjet_E
    xU[3] = TMath::Pi(); //TopLep, Neut_Theta
    xU[4] = 2*TMath::Pi(); //TopLep, Neut_Phi

}

void MultiLepton::FillTTHSemiLepHyp(MEPhaseSpace** meIntegrator)
{

  (*meIntegrator)->MEMFix_HiggsSemiLep.Lep1_E = Leptons[0].P4.E();
  (*meIntegrator)->MEMFix_HiggsSemiLep.Lep1_Theta = Leptons[0].P4.Theta();
  (*meIntegrator)->MEMFix_HiggsSemiLep.Lep1_Phi = Leptons[0].P4.Phi();
  (*meIntegrator)->MEMFix_HiggsSemiLep.Jet1_Theta = Jets[0].P4.Theta();
  (*meIntegrator)->MEMFix_HiggsSemiLep.Jet1_Phi = Jets[0].P4.Phi();
  (*meIntegrator)->MEMFix_HiggsSemiLep.Jet2_Theta = Jets[1].P4.Theta();
  (*meIntegrator)->MEMFix_HiggsSemiLep.Jet2_Phi = Jets[1].P4.Phi();
  (*meIntegrator)->MEMFix_HiggsSemiLep.LepSign = (Leptons[0].Id>0)?1:-1;

  (*meIntegrator)->MEMFix_TopLep.Bjet_Theta = Bjets[0].P4.Theta();
  (*meIntegrator)->MEMFix_TopLep.Bjet_Phi = Bjets[0].P4.Phi();
  (*meIntegrator)->MEMFix_TopLep.Lep_E = Leptons[1].P4.E();
  (*meIntegrator)->MEMFix_TopLep.Lep_Theta = Leptons[1].P4.Theta();
  (*meIntegrator)->MEMFix_TopLep.Lep_Phi = Leptons[1].P4.Phi();
  (*meIntegrator)->MEMFix_TopLep.TopSign = (Leptons[1].Id>0)?1:-1;

  (*meIntegrator)->MEMFix_TopLep2.Bjet_Theta = Bjets[1].P4.Theta();
  (*meIntegrator)->MEMFix_TopLep2.Bjet_Phi = Bjets[1].P4.Phi();
  (*meIntegrator)->MEMFix_TopLep2.Lep_E = Leptons[2].P4.E();
  (*meIntegrator)->MEMFix_TopLep2.Lep_Theta = Leptons[2].P4.Theta();
  (*meIntegrator)->MEMFix_TopLep2.Lep_Phi = Leptons[2].P4.Phi();
  (*meIntegrator)->MEMFix_TopLep2.TopSign = (Leptons[2].Id>0)?1:-1;

  (*meIntegrator)->MeasuredVarForTF.Jet1_E = Jets[0].P4.E();
  (*meIntegrator)->MeasuredVarForTF.Jet1_Eta = Jets[0].P4.Eta();
  (*meIntegrator)->MeasuredVarForTF.Jet2_E = Jets[1].P4.E();
  (*meIntegrator)->MeasuredVarForTF.Jet2_Eta = Jets[1].P4.Eta();
  (*meIntegrator)->MeasuredVarForTF.Bjet1_E =  Bjets[0].P4.E();
  (*meIntegrator)->MeasuredVarForTF.Bjet1_Eta =  Bjets[0].P4.Eta();
  (*meIntegrator)->MeasuredVarForTF.Bjet2_E =  Bjets[1].P4.E();
  (*meIntegrator)->MeasuredVarForTF.Bjet2_Eta =  Bjets[1].P4.Eta();
  (*meIntegrator)->MeasuredVarForTF.Recoil_Px = -Ptot.Px();
  (*meIntegrator)->MeasuredVarForTF.Recoil_Py = -Ptot.Py();
  (*meIntegrator)->MeasuredVarForTF.mET_Px = mET.Px();
  (*meIntegrator)->MeasuredVarForTF.mET_Py = mET.Py();
  
    xL[0] = (Bjets[0].P4.E()*JetTFfracmin<mB)?mB:Bjets[0].P4.E()*JetTFfracmin; //TopLep, Bjet_E
    xL[1] = 0; //TopLep, Neut_Theta
    xL[2] = 0; //TopLep, Neut_Phi
    xL[3] = (Bjets[1].P4.E()*JetTFfracmin<mB)?mB:Bjets[1].P4.E()*JetTFfracmin; //TopLep2, Bjet_E
    xL[4] = 0; //TopLep2, Neut_Theta
    xL[5] = 0; //TopLep2, Neut_Phi
    xL[6] = 0; //HiggsSemiLep, Neut_Theta
    xL[7] = 0; //HiggsSemiLep, Neut_Phi
    xL[8] = Jets[0].P4.E()*JetTFfracmin; //HiggsSemiLep, Jet1_E
    xL[9] = Jets[1].P4.E()*JetTFfracmin; //HiggsSemiLep, Jet2_E

    xU[0] = JetTFfracmax*Bjets[0].P4.E();//meIntegrator->comEnergy; //TopLep, Bjet_E
    xU[1] = TMath::Pi(); //TopLep, Neut_Theta
    xU[2] = 2.*TMath::Pi(); //TopLep, Neut_Phi
    xU[3] = JetTFfracmax*Bjets[1].P4.E();//meIntegrator->comEnergy; //TopLep2, Bjet_E
    xU[4] = TMath::Pi(); //TopLep2, Neut_Theta
    xU[5] = 2*TMath::Pi(); //TopLep2, Neut_Phi
    xU[6] = TMath::Pi(); //HiggsSemiLep, Neut_Theta
    xU[7] = 2*TMath::Pi(); //HiggsSemiLep, Neut_Phi
    xU[8] = JetTFfracmax*Jets[0].P4.E();//meIntegrator->comEnergy; //HiggsSemiLep, Jet1_E
    xU[9] = JetTFfracmax*Jets[1].P4.E();//meIntegrator->comEnergy //HiggsSemiLep, Jet2_E

}

void MultiLepton::FillTTWHyp(MEPhaseSpace** meIntegrator)
{

  (*meIntegrator)->MEMFix_HiggsSemiLep.Lep1_E = Leptons[0].P4.E();
  (*meIntegrator)->MEMFix_HiggsSemiLep.Lep1_Theta = Leptons[0].P4.Theta();
  (*meIntegrator)->MEMFix_HiggsSemiLep.Lep1_Phi = Leptons[0].P4.Phi();
  (*meIntegrator)->MEMFix_HiggsSemiLep.Jet1_Theta = Jets[0].P4.Theta();
  (*meIntegrator)->MEMFix_HiggsSemiLep.Jet1_Phi = Jets[0].P4.Phi();
  (*meIntegrator)->MEMFix_HiggsSemiLep.Jet2_Theta = Jets[1].P4.Theta();
  (*meIntegrator)->MEMFix_HiggsSemiLep.Jet2_Phi = Jets[1].P4.Phi();
  (*meIntegrator)->MEMFix_HiggsSemiLep.LepSign = (Leptons[0].Id>0)?1:-1;

  (*meIntegrator)->MEMFix_TopLep.Bjet_Theta = Bjets[0].P4.Theta();
  (*meIntegrator)->MEMFix_TopLep.Bjet_Phi = Bjets[0].P4.Phi();
  (*meIntegrator)->MEMFix_TopLep.Lep_E = Leptons[1].P4.E();
  (*meIntegrator)->MEMFix_TopLep.Lep_Theta = Leptons[1].P4.Theta();
  (*meIntegrator)->MEMFix_TopLep.Lep_Phi = Leptons[1].P4.Phi();
  (*meIntegrator)->MEMFix_TopLep.TopSign = (Leptons[1].Id>0)?1:-1;

  (*meIntegrator)->MEMFix_TopLep2.Bjet_Theta = Bjets[1].P4.Theta();
  (*meIntegrator)->MEMFix_TopLep2.Bjet_Phi = Bjets[1].P4.Phi();
  (*meIntegrator)->MEMFix_TopLep2.Lep_E = Leptons[2].P4.E();
  (*meIntegrator)->MEMFix_TopLep2.Lep_Theta = Leptons[2].P4.Theta();
  (*meIntegrator)->MEMFix_TopLep2.Lep_Phi = Leptons[2].P4.Phi();
  (*meIntegrator)->MEMFix_TopLep2.TopSign = (Leptons[2].Id>0)?1:-1;

  (*meIntegrator)->MeasuredVarForTF.Jet1_E = Jets[0].P4.E();
  (*meIntegrator)->MeasuredVarForTF.Jet1_Eta = Jets[0].P4.Eta();
  (*meIntegrator)->MeasuredVarForTF.Jet2_E = Jets[1].P4.E();
  (*meIntegrator)->MeasuredVarForTF.Jet2_Eta = Jets[1].P4.Eta();
  (*meIntegrator)->MeasuredVarForTF.Bjet1_E =  Bjets[0].P4.E();
  (*meIntegrator)->MeasuredVarForTF.Bjet1_Eta =  Bjets[0].P4.Eta();
  (*meIntegrator)->MeasuredVarForTF.Bjet2_E =  Bjets[1].P4.E();
  (*meIntegrator)->MeasuredVarForTF.Bjet2_Eta =  Bjets[1].P4.Eta();
  (*meIntegrator)->MeasuredVarForTF.Recoil_Px = -Ptot.Px();
  (*meIntegrator)->MeasuredVarForTF.Recoil_Py = -Ptot.Py();
  (*meIntegrator)->MeasuredVarForTF.mET_Px = mET.Px();
  (*meIntegrator)->MeasuredVarForTF.mET_Py = mET.Py();


    xL[0] = (Bjets[0].P4.E()*JetTFfracmin<mB)?mB:Bjets[0].P4.E()*JetTFfracmin; //TopLep, Bjet_E
    xL[1] = 0; //TopLep, Neut_Theta
    xL[2] = 0; //TopLep, Neut_Phi
    xL[3] = (Bjets[1].P4.E()*JetTFfracmin<mB)?mB:Bjets[1].P4.E()*JetTFfracmin; //TopLep2, Bjet_E
    xL[4] = 0; //TopLep2, Neut_Theta
    xL[5] = 0; //TopLep2, Neut_Phi
    xL[6] = 0;//(74.-Leptons[0].P4.E()>0)?(74.-Leptons[0].P4.E()):0; //HiggsSemiLep, Neut_E
    xL[7] = 0; //HiggsSemiLep, Neut_Theta
    xL[8] = 0; //HiggsSemiLep, Neut_Phi
    xL[9] = Jets[0].P4.E()*JetTFfracmin; //Jet1_E
    xL[10] = Jets[1].P4.E()*JetTFfracmin; //Jet2_E

    xU[0] = JetTFfracmax*Bjets[0].P4.E();//meIntegrator->comEnergy; //TopLep, Bjet_E
    xU[1] = TMath::Pi(); //TopLep, Neut_Theta
    xU[2] = 2.*TMath::Pi(); //TopLep, Neut_Phi
    xU[3] = JetTFfracmax*Bjets[1].P4.E();//meIntegrator->comEnergy; //TopLep2, Bjet_E
    xU[4] = TMath::Pi(); //TopLep2, Neut_Theta
    xU[5] = 2*TMath::Pi(); //TopLep2, Neut_Phi
    xU[6] = NeutMaxE; //HiggsSemiLep, Neut_E
    xU[7] = TMath::Pi(); //HiggsSemiLep, Neut_Theta
    xU[8] = 2*TMath::Pi();//HiggsSemiLep, Neut_Phi
    xU[9] = JetTFfracmax*Jets[0].P4.E(); //Jet1_E
    xU[10] = JetTFfracmax*Jets[1].P4.E(); //Jet2_E

  return;
}

#endif
