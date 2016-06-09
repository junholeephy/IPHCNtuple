#include "include/NtupleProducer.h"
#include <string.h>

ClassImp(Event)

Event::Event()
{
}

Event::~Event()
{
}

void Event::read(bool isdata)
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
    _metcov00	      = ntP->met_cov00;
    _metcov01         = ntP->met_cov01;
    _metcov10         = ntP->met_cov10;
    _metcov11         = ntP->met_cov11;

    _metNoHF_pt        = ntP->metNoHF_pt;
    _metNoHF_phi       = ntP->metNoHF_phi;
    _metNoHF_sumet     = ntP->metNoHF_sumet;

    if (!isdata)
    {

        _weight_scale_muF0p5 = ntP->weight_scale_muF0p5;
        _weight_scale_muF2   = ntP->weight_scale_muF2;
        _weight_scale_muR0p5 = ntP->weight_scale_muR0p5;
        _weight_scale_muR2   = ntP->weight_scale_muR2;

        _mc_weight           = ntP->mc_weight;
        _mc_ptHat            = ntP->mc_ptHat;
        _mc_pu_trueNumInt    = ntP->mc_pu_trueNumInt;
	
        _pdf_weights = *ntP->mc_pdfweights;
	//_pdf_ids = *ntP->mc_pdfids;
	
	}

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

        trigger_comb        = 0;
        if( ntP->trigger_pass->at(61) == 1) {trigger_comb = trigger_comb + 1     ;} // [61]  HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v1
        if( ntP->trigger_pass->at(65) == 1) {trigger_comb = trigger_comb + 2     ;} // [65]  HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1
        if( ntP->trigger_pass->at(63) == 1) {trigger_comb = trigger_comb + 5     ;} // [63]  HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1
        if( ntP->trigger_pass->at(60) == 1) {trigger_comb = trigger_comb + 10    ;} // [60]  HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1
        if( ntP->trigger_pass->at(22) == 1) {trigger_comb = trigger_comb + 20    ;} // [22]  HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1
        if( ntP->trigger_pass->at(24) == 1) {trigger_comb = trigger_comb + 50    ;} // [24]  HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1
        if( ntP->trigger_pass->at(10) == 1) {trigger_comb = trigger_comb + 100   ;} // [10]  HLT_IsoMu20_v1
        if( ntP->trigger_pass->at(15) == 1) {trigger_comb = trigger_comb + 200   ;} // [15]  HLT_IsoTkMu20_v1
        if( ntP->trigger_pass->at(69) == 1) {trigger_comb = trigger_comb + 500   ;} // [69]  HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v1
        // for data HLT_Ele23_WPLoose_Gsf_v1

        int trigger_pass_byname = 0;
        for( int i = 0; i < ntP->trigger->size(); i++)
        {
            if ( ntP->trigger_name->at(i) == "HLT_Mu17_Mu8_DZ_v1"                                ) {trigger_pass_byname = trigger_pass_byname + 1     ;}
            if ( ntP->trigger_name->at(i) == "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1"      ) {trigger_pass_byname = trigger_pass_byname + 10    ;}
            if ( ntP->trigger_name->at(i) == "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1" ) {trigger_pass_byname = trigger_pass_byname + 100   ;}
            if ( ntP->trigger_name->at(i) == "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1") {trigger_pass_byname = trigger_pass_byname + 1000  ;}
        }

        trigger_pass_byname     = 0;
        for( int i = 0; i < ntP->trigger->size(); i++)
        {
            std::string currentpath ("Nopathsofar");
            if( ntP->trigger_pass->at(i) == 1) { currentpath = ntP->trigger_name->at(i); }

            std::string eee   ("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v1"            );

            std::string me    ("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1"  );
            std::string em    ("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1"   );
            std::string ee    ("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1"        );
            std::string mm    ("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1"              );
            std::string mmTk  ("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1"            );
            std::string m     ("HLT_IsoMu20_v1"                                      );
            std::string mTk   ("HLT_IsoTkMu20_v1"                                    );
            std::string e     ("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v1"                 );
            std::string eData ("HLT_Ele23_WPLoose_Gsf_v1"                            );

            std::string mme   ("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v1"                  );
            std::string eem   ("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v1"                 );
            std::string mmm   ("HLT_TripleMu_12_10_5_v1"                             );

            if(currentpath.compare(eee)  == 0) {trigger_pass_byname = trigger_pass_byname + 1    ;}
            if(currentpath.compare(me)   == 0) {trigger_pass_byname = trigger_pass_byname + 2    ;}
            if(currentpath.compare(em)   == 0) {trigger_pass_byname = trigger_pass_byname + 5    ;}
            if(currentpath.compare(ee)   == 0) {trigger_pass_byname = trigger_pass_byname + 10   ;}
            if(currentpath.compare(mm)   == 0) {trigger_pass_byname = trigger_pass_byname + 20   ;}
            if(currentpath.compare(mmTk) == 0) {trigger_pass_byname = trigger_pass_byname + 50   ;}
            if(currentpath.compare(m)    == 0) {trigger_pass_byname = trigger_pass_byname + 100  ;}
            if(currentpath.compare(mTk)  == 0) {trigger_pass_byname = trigger_pass_byname + 200  ;}
            if(currentpath.compare(e)    == 0) {trigger_pass_byname = trigger_pass_byname + 500  ;}

            if(currentpath.compare(eData)  == 0) {trigger_pass_byname = trigger_pass_byname + 1000   ;}
            if(currentpath.compare(mme)    == 0) {trigger_pass_byname = trigger_pass_byname + 2000   ;}
            if(currentpath.compare(eem)    == 0) {trigger_pass_byname = trigger_pass_byname + 5000   ;}
            if(currentpath.compare(mmm)    == 0) {trigger_pass_byname = trigger_pass_byname + 10000  ;}
        }

        int trigger_pass_byname_1      = 0; 
        int trigger_pass_byname_1_noDz = 0;

        for( int i = 0; i < ntP->trigger->size(); i++)
        {
            std::string currentpath ("Nopathsofar");
            if( ntP->trigger_pass->at(i) == 1) { currentpath = ntP->trigger_name->at(i); }

            std::size_t eee   = currentpath.find("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v"            );
            std::size_t me    = currentpath.find("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"  );
            std::size_t em    = currentpath.find("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v"   );
            std::size_t ee    = currentpath.find("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"        );
            std::size_t mm    = currentpath.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"              );
            std::size_t mmTk  = currentpath.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"            );
            std::size_t m     = currentpath.find("HLT_IsoMu20_v"                                      );
            std::size_t mTk   = currentpath.find("HLT_IsoTkMu20_v"                                    );
            std::size_t e     = currentpath.find("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v"                 );
            std::size_t eData = currentpath.find("HLT_Ele23_WPLoose_Gsf_v"                            );
            std::size_t mme   = currentpath.find("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v"                  );
            std::size_t eem   = currentpath.find("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v"                 );
            std::size_t mmm   = currentpath.find("HLT_TripleMu_12_10_5_v"                             );

            std::size_t ee_noDz   = currentpath.find("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v"	  );
            std::size_t mm_noDz   = currentpath.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v"		  );
            std::size_t mmTk_noDz = currentpath.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v" 	  );

            if(eee  != std::string::npos) {trigger_pass_byname_1 = trigger_pass_byname_1 + 1     ;}
            if(me   != std::string::npos) {trigger_pass_byname_1 = trigger_pass_byname_1 + 2     ;}
            if(em   != std::string::npos) {trigger_pass_byname_1 = trigger_pass_byname_1 + 5     ;}
            if(ee   != std::string::npos) {trigger_pass_byname_1 = trigger_pass_byname_1 + 10    ;}
            if(mm   != std::string::npos) {trigger_pass_byname_1 = trigger_pass_byname_1 + 20    ;}
            if(mmTk != std::string::npos) {trigger_pass_byname_1 = trigger_pass_byname_1 + 50    ;}
            if(m    != std::string::npos) {trigger_pass_byname_1 = trigger_pass_byname_1 + 100   ;}
            if(mTk  != std::string::npos) {trigger_pass_byname_1 = trigger_pass_byname_1 + 200   ;}
            if(e    != std::string::npos) {trigger_pass_byname_1 = trigger_pass_byname_1 + 500   ;}
            if(eData !=std::string::npos) {trigger_pass_byname_1 = trigger_pass_byname_1 + 1000  ;}
            if(mme   !=std::string::npos) {trigger_pass_byname_1 = trigger_pass_byname_1 + 2000  ;}
            if(eem   !=std::string::npos) {trigger_pass_byname_1 = trigger_pass_byname_1 + 5000  ;}
            if(mmm   !=std::string::npos) {trigger_pass_byname_1 = trigger_pass_byname_1 + 10000 ;}

            if(ee_noDz   != std::string::npos) {trigger_pass_byname_1_noDz = trigger_pass_byname_1_noDz + 10 ;}
            if(mm_noDz   != std::string::npos) {trigger_pass_byname_1_noDz = trigger_pass_byname_1_noDz + 20 ;}
            if(mmTk_noDz != std::string::npos) {trigger_pass_byname_1_noDz = trigger_pass_byname_1_noDz + 50 ;}


        }

        //std::cout << "Trigger combination = " << trigger_comb << std::endl;
        _trigger_pass               = trigger_comb;
        _trigger_pass_byname        = trigger_pass_byname;
        _trigger_pass_byname_1      = trigger_pass_byname_1;
        _trigger_pass_byname_1_noDz = trigger_pass_byname_1_noDz;

        // discriminant vs tt
        _disc_TT = 0;
}

void Event::init()
{
    _id                    = -888;
    _run                   = -888;
    _lumi                  = -888;
    _rho                   = -888;

    _pv_n                  = -1;
    _pv_z                  = -888;
    _pv_zError             = -888;

    _metpt                 = -888;
    _metphi                = -888;
    _metsumet              = -888;
    _metcov00              = -888;
    _metcov01              = -888;
    _metcov10              = -888;
    _metcov11              = -888;

    _metNoHF_pt            = -888;
    _metNoHF_phi           = -888;
    _metNoHF_sumet         = -888; 

    _weight_scale_muF0p5   = -888;
    _weight_scale_muF2     = -888;
    _weight_scale_muR0p5   = -888;
    _weight_scale_muR2     = -888;

    _mc_weight             = -888;
    _mc_ptHat              = -888;
    _mc_pu_trueNumInt      = -888;

    _trigger_pass          = -888;
    _trigger_pass_byname   = -888;
    _trigger_pass_byname_1 = -888;

    _tth_channel           = -888;

    _disc_TT               = -888;
    
    _pdf_weights.clear();
    //_pdf_ids.clear();
}
