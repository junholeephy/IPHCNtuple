#include "include/NtupleProducer.h"
#include <iostream>
#include <iomanip>

ClassImp(Muon)
    
Muon::Muon()
{
}

Muon::~Muon()
{
}

void Muon::read()
{
    _ID = idx;

    // general informations
    if( CHECK(ntP->mu_E)                 ) _E        = ntP->mu_E->at(idx);
    if( CHECK(ntP->mu_pt)                ) _pt       = ntP->mu_pt->at(idx);
    if( CHECK(ntP->mu_eta)               ) _eta      = ntP->mu_eta->at(idx);   
    if( CHECK(ntP->mu_phi)               ) _phi      = ntP->mu_phi->at(idx);   
    if( CHECK(ntP->mu_m)                 ) _m        = ntP->mu_m->at(idx);

    // preselection variables
    if( CHECK(ntP->mu_innerTrack_PV_dxy) ) _dxy      = ntP->mu_innerTrack_PV_dxy->at(idx);
    if( CHECK(ntP->mu_innerTrack_PV_dz)  ) _dz       = ntP->mu_innerTrack_PV_dz->at(idx);
    // placeholder miniIso
    // placeholder SIP
    if( CHECK(ntP->mu_isPFMuon)          ) _isPFMuon = ntP->mu_isPFMuon->at(idx);

    // cut-based selection additionnal variables
    if( CHECK(ntP->mu_lepMVA_sip3d)      ) _sip3d    = ntP->mu_lepMVA_sip3d->at(idx);
    // placholder Medium ID
    // placeholder track.pt/errortrack.pt

    // mva-based selection additionnal variables
    // placeholder personal MVA

    // jet related variables
    //_jetPtRel   = ntP->mu_lepMVA_jetPtRelv2->at(idx);
    //_jetCSV     = ntP->mu_lepMVA_jetBTagCSV->at(idx);
    //_jetPtRatio = ntP->mu_lepMVA_jetPtRatio->at(idx);

    // more variables
    if( CHECK(ntP->mu_id)                ) _id       = ntP->mu_id->at(idx);

}

void Muon::init()
{
    // general informations
    _E                 = -666;
    _pt                = -666;
    _eta               = -666;
    _phi               = -666;
    _m                 = -666;
    _charge            = 0;
    _id                = 0;    

    // Id
    _isLoose           = 0;
    _isTight           = 0;
    _isLooseMVA        = 0;
    _isTightMVA        = 0;

    // variables for Id
    _dxy               = -666;
    _dz                = -666;
    _iso               = -666;

    // more variables
    _sip3d             = -666;
}

bool Muon::sel()
{

    float miniIso            = ntP->mu_miniIsoTTH->at(idx);
    float sip3d              = ntP->mu_ip3d->at(idx) / ntP->mu_ip3dErr->at(idx);
    bool  isLoose            = ntP->mu_isLooseMuon->at(idx);

    // preselection
    bool pass_pt      = (_pt         >  5    );
    bool pass_eta     = (fabs(_eta)  <  2.4  );
    bool pass_dxy     = (fabs(_dxy)  <  0.05 );
    bool pass_dz      = (fabs(_dz)   <  0.1  );
    bool pass_miniIso = (miniIso     <  0.4  );
    bool pass_SIP     = (fabs(sip3d) <  8    );
    bool pass_isLoose = (isLoose             );

    bool isPreselectionMuon = ( pass_pt      &&
                                pass_eta     &&
                                pass_dxy     &&
                                pass_dz      &&
                                pass_miniIso &&
                                pass_SIP     &&
                                pass_isLoose );

    // mva-based selection
    //bool pass_mva      = (_mva       > 0.65);
    //bool pass_isMedium = (_isMedium        );

    //bool isSelectionElectron        = ( isPreselectionElectron &&
                                        //pass_mva               &&
                                        //pass_isMedium          &&
    //                                  )

    cout<<std::setiosflags(ios::fixed)<<setprecision(5);
    // synchronization printout
    /*if( isPreselectionMuon ) std::cout                                              << std::setw(10)
                                       << nt->NtEvent->at(0).id()                   << std::setw(10)
                                       << _pt                                       << std::setw(10)
                                       << _eta                                      << std::setw(10)
                                       << _phi                                      << std::setw(10)
                                       << _E                                        << std::setw(5)
                                       << ntP->mu_id->at(idx)                       << std::setw(5)
                                       << ntP->mu_charge->at(idx)                   << std::setw(15)
                                       << miniIso                                   << std::setw(15)
                                       << ntP->mu_lepMVA_miniRelIsoCharged->at(idx) << std::setw(15)
                                       << ntP->mu_lepMVA_miniRelIsoNeutral->at(idx) << std::setw(10)
                                       << ntP->mu_lepMVA_jetPtRelv2->at(idx)        << std::setw(10)
                                       //<< 0.0                                       << std::setw(10)
                                       << ntP->mu_lepMVA_jetBTagCSV->at(idx)        << std::setw(10)
                                       << ntP->mu_lepMVA_jetPtRatio->at(idx)        << std::setw(10)
                                       << fabs(sip3d)                               << std::setw(10)
                                       << fabs(_dxy)                                << std::setw(10)
                                       << _dz                                       << std::setw(21)
                                       << ntP->mu_segmentCompatibility->at(idx)     << std::endl;
    */
    return isPreselectionMuon;
}

float Muon::effectiveArea(int dr,float eta)
{
    float effArea = 0.;

    if( fabs(eta) >= 0 && fabs(eta) < 0.8 )
    {
        effArea = 0.0735;
    }   
    else if( fabs(eta) >= 0.8 && fabs(eta) < 1.3 )
    {
        effArea = 0.0619;
    }   
    else if( fabs(eta) >= 1.3 && fabs(eta) < 2.0 )
    {
        effArea = 0.0465;
    }   
    else if( fabs(eta) >= 2.0 && fabs(eta) < 2.2 )
    {
        effArea = 0.0433;
    }   
    else if( fabs(eta) >= 2.2 && fabs(eta) < 2.5 )
    {
        effArea = 0.0577;
    }   

    return effArea;
}
