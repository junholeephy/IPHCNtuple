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

   //std::cout << "Taille menu: " << ntP->trigger_pass->size() << std::endl;
   
   //for( int i = 0; i < ntP->trigger->size(); i++)
   //{
   //    std::cout << "Trigger [" << i << "]: " << ntP->trigger_name->at(i) << std::endl;
   //}
       
   if( ntP->trigger_pass->at(19) == 1) {trigger_comb = trigger_comb + 1   ;} // [19]  HLT_Mu17_Mu8_DZ_v1
   if( ntP->trigger_pass->at(60) == 1) {trigger_comb = trigger_comb + 10  ;} // [60] HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1
   if( ntP->trigger_pass->at(63) == 1) {trigger_comb = trigger_comb + 100 ;} // [63] HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1
   if( ntP->trigger_pass->at(65) == 1) {trigger_comb = trigger_comb + 1000;} // [65] HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1
    //for more triggers, see /opt/sbg/data/data2/cms/xcoubez/PhD2ndYear/Production/Production_7_4_12_patch4/dev/JetInformation/CMSSW_7_4_12_patch4/src/IPHCFlatTree/FlatTreeProducer/test/TriggerList.out...

   //std::cout << "Trigger combination = " << trigger_comb << std::endl;
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
