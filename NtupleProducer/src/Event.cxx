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
    // Event
    _id               = ntP->ev_id;
    _run              = ntP->ev_run;
    _lumi             = ntP->ev_lumi;
    _rho              = ntP->ev_rho;

    // pv
    _pv_n             = ntP->nvertex;
    _pv_z             = ntP->pv_z;
    _pv_zError        = ntP->pv_zError;

    // MET
    _metpt            = ntP->met_pt;
    _metphi           = ntP->met_phi;
    _metsumet         = ntP->met_sumet;

    _metNoHF_pt        = ntP->metNoHF_pt;
    _metNoHF_phi       = ntP->metNoHF_phi;
    _metNoHF_sumet     = ntP->metNoHF_sumet;

    _mc_weight        = ntP->mc_weight;
    _mc_ptHat         = ntP->mc_ptHat;
    _mc_pu_trueNumInt = ntP->mc_pu_trueNumInt;

    // trigger

    //std::cout << "Taille menu: " << ntP->trigger_pass->size() << std::endl;

    //for( int i = 0; i < ntP->trigger->size(); i++)
    //{
    //    std::cout << "Trigger [" << i << "]: " << ntP->trigger_name->at(i) << std::endl;
    //}

    int trigger_comb        = 0;

    if( ntP->trigger_pass->at(19) == 1) {trigger_comb = trigger_comb + 1   ;} // [19]  HLT_Mu17_Mu8_DZ_v1
    if( ntP->trigger_pass->at(60) == 1) {trigger_comb = trigger_comb + 10  ;} // [60] HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1
    if( ntP->trigger_pass->at(63) == 1) {trigger_comb = trigger_comb + 100 ;} // [63] HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1
    if( ntP->trigger_pass->at(65) == 1) {trigger_comb = trigger_comb + 1000;} // [65] HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1
    //for more triggers, see /opt/sbg/data/data2/cms/xcoubez/PhD2ndYear/Production/Production_7_4_12_patch4/dev/JetInformation/CMSSW_7_4_12_patch4/src/IPHCFlatTree/FlatTreeProducer/test/TriggerList.out...

    int trigger_pass_byname = 0;

    for( int i = 0; i < ntP->trigger->size(); i++)
    {
        if ( ntP->trigger_name->at(i) == "HLT_Mu17_Mu8_DZ_v1"                                ) {trigger_pass_byname = trigger_pass_byname + 1     ;}
        if ( ntP->trigger_name->at(i) == "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1"      ) {trigger_pass_byname = trigger_pass_byname + 10    ;}
        if ( ntP->trigger_name->at(i) == "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1" ) {trigger_pass_byname = trigger_pass_byname + 100   ;}
        if ( ntP->trigger_name->at(i) == "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1") {trigger_pass_byname = trigger_pass_byname + 1000  ;}
    }

    //std::cout << "Trigger combination = " << trigger_comb << std::endl;
    _trigger_pass        = trigger_comb;
    _trigger_pass_byname = trigger_pass_byname;

    // discriminant vs tt
    _disc_TT = 0;
}

void Event::init()
{
    _id               = -666;
    _run              = -666;
    _lumi             = -666;
    _rho              = -666;

    _pv_n             = -1;
    _pv_z             = -666;
    _pv_zError        = -666;

    _metpt            = -666;
    _metphi           = -666;
    _metsumet         = -666;

    _metNoHF_pt        = -888;
    _metNoHF_phi       = -888;
    _metNoHF_sumet     = -888; 

    _mc_weight        = -888;
    _mc_ptHat         = -888;
    _mc_pu_trueNumInt = -888;

    _trigger_pass        = -888;
    _trigger_pass_byname = -888;

    _tth_channel      = -888;

    _disc_TT          = -888;
}
