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
    _E        = ntP->mu_E->at(idx);
    _pt       = ntP->mu_pt->at(idx);
    _eta      = ntP->mu_eta->at(idx);	
    _phi      = ntP->mu_phi->at(idx);	
    _m        = ntP->mu_m->at(idx);
    _charge   = ntP->mu_charge->at(idx);
    _id       = ntP->mu_id->at(idx);

    // 
    _isLoose           = ntP->mu_isLooseMuon->at(idx);
    _isMedium	       = ntP->mu_isMediumMuon->at(idx);
    _isTight	       = ntP->mu_isTightMuon->at(idx);
    _isPFMuon          = ntP->mu_isPFMuon->at(idx);

    // preselection variables
    _dxy               = ntP->mu_innerTrack_PV_dxy->at(idx);
    _dz                = ntP->mu_innerTrack_PV_dz->at(idx);
    _iso               = ntP->mu_miniIsoTTH->at(idx);
    _sip3d             = ntP->mu_ip3d->at(idx) / ntP->mu_ip3dErr->at(idx);
    _bestTrack_pt      = ntP->mu_bestTrack_pt->at(idx);
    _bestTrack_ptError = ntP->mu_bestTrack_ptError->at(idx);
    _dB3D	       = ntP->mu_dB3D->at(idx);
    _edB3D	       = ntP->mu_edB3D->at(idx);
   
    // mva-based selection additionnal variables
    _lepMVA	      = ntP->mu_lepMVA->at(idx);
    _lepMVA_Moriond16 = ntP->mu_lepMVA_Moriond16->at(idx);

    _lepMVA_miniRelIsoCharged    = ntP->mu_lepMVA_miniRelIsoCharged->at(idx);
    _lepMVA_miniRelIsoNeutral    = ntP->mu_lepMVA_miniRelIsoNeutral->at(idx);
    _lepMVA_jetPtRelv2           = ntP->mu_lepMVA_jetPtRelv2->at(idx);  
    _lepMVA_jetDR                = ntP->mu_lepMVA_jetDR->at(idx);
    _lepMVA_jetPtRatio           = ntP->mu_lepMVA_jetPtRatio->at(idx);
    _lepMVA_jetBTagCSV           = ntP->mu_lepMVA_jetBTagCSV->at(idx);
    _lepMVA_sip3d                = ntP->mu_lepMVA_sip3d->at(idx);
    _lepMVA_dxy                  = ntP->mu_lepMVA_dxy->at(idx);
    _lepMVA_dz                   = ntP->mu_lepMVA_dz->at(idx);
    _lepMVA_mvaId                = ntP->mu_lepMVA_mvaId->at(idx);
    _lepMVA_eta                  = ntP->mu_lepMVA_eta->at(idx);
    _lepMVA_jetNDauChargedMVASel = ntP->mu_lepMVA_jetNDauChargedMVASel->at(idx);
 
    // more variables

}

void Muon::init()
{
    // general informations
    
    _fakeType          = -666.;
    
    _E                 = -666.;
    _pt                = -666.;
    _eta               = -666.;
    _phi               = -666.;
    _m                 = -666.;
    _charge            = 0;
    _id                = 0;    

    // Id
    _isLoose           = 0;
    _isMedium          = 0;
    _isTight           = 0;
    _isPFMuon          = 0;
  
    _isLooseTTH        = 0;
    _isFakeableTTH     = 0;
    _isTightTTH        = 0;
    
    // variables for Id
    
    _dxy                = -666.;
    _dz                 = -666.;
    _iso                = -666.; 
    _sip3d              = -666.;   
    _bestTrack_pt       = -666.;
    _bestTrack_ptError  = -666.;
    _dB3D               = -666.;
    _edB3D              = -666.;
    
    // more variables
   
    _lepMVA            = -666.; 
    _lepMVA_Moriond16  = -666.; 
    
      
    _lepMVA_miniRelIsoCharged    = -666.;
    _lepMVA_miniRelIsoNeutral    = -666.;
    _lepMVA_jetPtRelv2           = -666.;
    _lepMVA_jetDR                = -666.;
    _lepMVA_jetPtRatio           = -666.;
    _lepMVA_jetBTagCSV           = -666.;
    _lepMVA_sip3d                = -666.;
    _lepMVA_dxy                  = -666.;
    _lepMVA_dz                   = -666.;
    _lepMVA_mvaId                = -666.;
    _lepMVA_eta                  = -666.;
    _lepMVA_jetNDauChargedMVASel = -666.;
           
}

bool Muon::sel()
{

    // Loose
    bool pass_pt      = (_pt          >  5    );
    bool pass_eta     = (fabs(_eta)   <  2.4  );
    bool pass_dxy     = (fabs(_dxy)   <  0.05 );
    bool pass_dz      = (fabs(_dz)    <  0.1  );
    bool pass_miniIso = (_iso         <  0.4  );
    bool pass_SIP     = (fabs(_sip3d) <  8    );
    bool pass_isLoose = (_isLoose             );

    bool isLooseTTH = ( pass_pt      &&
                        pass_eta     &&
                        pass_dxy     &&
                        pass_dz      &&
                        pass_miniIso &&
                        pass_SIP     &&
                        pass_isLoose );

    _isLooseTTH = isLooseTTH;

    
    // Fakeable
    
  
    bool pass_lepMVA_Moriond16  = _lepMVA_Moriond16 > 0.75 ;
    bool pass_lepMVA_jetBTagCSV089 = _lepMVA_jetBTagCSV < 0.89;
    
    bool pass_lepMVA_jetBtagCSVPtRatio = false;
    
    if (!pass_lepMVA_Moriond16 && _lepMVA_jetPtRatio > 0.3 && _lepMVA_jetBTagCSV < 0.605) pass_lepMVA_jetBtagCSVPtRatio = true;
    if (pass_lepMVA_Moriond16 && pass_lepMVA_jetBTagCSV089) pass_lepMVA_jetBtagCSVPtRatio = true;
    
    _isFakeableTTH = isLooseTTH && pass_lepMVA_jetBtagCSVPtRatio;

    // Tight
  
    _isTightTTH = isLooseTTH && pass_lepMVA_Moriond16 && _isMedium && pass_lepMVA_jetBTagCSV089;

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
    return isLooseTTH;
}

