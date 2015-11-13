#include "include/NtupleProducer.h"

ClassImp(Event)
    
Event::Event()
{
}

Event::~Event()
{
}

void Event::read()
{
   _id               = ntP->ev_id;
   _run              = ntP->ev_run;
   _lumi             = ntP->ev_lumi;
   _rho              = ntP->ev_rho;

   _pv_n             = ntP->nvertex;
   _pv_z             = ntP->pv_z;
   _pv_zError        = ntP->pv_zError;
   
   _metpt            = ntP->met_pt;
   _metphi           = ntP->met_phi;
   _metsumet         = ntP->met_sumet;
   
   _mc_weight        = ntP->mc_weight;
   _mc_ptHat         = ntP->mc_ptHat;
   _mc_pu_trueNumInt = ntP->mc_pu_trueNumInt;

   int trigger_comb = 0;

   if( ntP->trigger_pass.at(94)  == 1) {trigger_comb = trigger_comb + 1   ;} // [94]  HLT_Mu17_Mu8_DZ_v1
   if( ntP->trigger_pass.at(216) == 1) {trigger_comb = trigger_comb + 10  ;} // [216] HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1
   if( ntP->trigger_pass.at(219) == 1) {trigger_comb = trigger_comb + 100 ;} // [219] HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1
   if( ntP->trigger_pass.at(221) == 1) {trigger_comb = trigger_comb + 1000;} // [221] HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1
   // for more triggers, see /opt/sbg/data/data2/cms/xcoubez/PhD2ndYear/Production/Production_7_4_12_patch4/dev/JetInformation/CMSSW_7_4_12_patch4/src/IPHCFlatTree/FlatTreeProducer/test/TriggerList.out...

   _trigger_pass = trigger_comb;
}

void Event::init()
{
   _id          = -666;
   _run         = -666;
   _lumi        = -666;
   _rho         = -666;

   _pv_n      = -1;
   _pv_z      = -666;
   _pv_zError = -666;

   _metpt       = -666;
   _metphi      = -666;
   _metsumet    = -666;
   
   _mc_weight   = -666;
   _mc_ptHat         = -666;
   _mc_pu_trueNumInt = -666;

    _trigger_pass = -666;

   _tth_channel = -666;
}
