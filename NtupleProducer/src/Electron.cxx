#include "include/NtupleProducer.h"
#include <iostream>
#include <iomanip>

ClassImp(Electron)
    
Electron::Electron()
{
}

Electron::~Electron()
{
}

void Electron::read()
{
    _ID = idx;

    // general informations
    if( CHECK(ntP->el_E)                             ) _E         = ntP->el_E->at(idx);
    if( CHECK(ntP->el_pt)                            ) _pt        = ntP->el_pt->at(idx);
    if( CHECK(ntP->el_eta)                           ) _eta       = ntP->el_eta->at(idx);
    if( CHECK(ntP->el_phi)                           ) _phi       = ntP->el_phi->at(idx);
    if( CHECK(ntP->el_m)                             ) _m         = ntP->el_m->at(idx);

    // preselection variables
    if( CHECK(ntP->el_gsfTrack_PV_dxy)               ) _dxy       = ntP->el_gsfTrack_PV_dxy->at(idx);
    if( CHECK(ntP->el_gsfTrack_PV_dz)                ) _dz        = ntP->el_gsfTrack_PV_dz->at(idx);
    // placeholder miniIso
    // placeholder SIP
    if( CHECK(ntP->el_looseCBId)                     ) _isLoose   = ntP->el_looseCBId->at(idx);
    if( CHECK(ntP->el_numberOfLostHits)              ) _nlosthits = ntP->el_numberOfLostHits->at(idx);

    // cut-based selection additionnal variables
    if( CHECK(ntP->el_lepMVA_sip3d)                  ) _sip3d     = ntP->el_lepMVA_sip3d->at(idx);
    if( CHECK(ntP->el_mediumCBId)                    ) _isMedium  = ntP->el_mediumCBId->at(idx);
    if( CHECK(ntP->el_passConversionVeto)            ) _passCV    = ntP->el_passConversionVeto->at(idx);
    if( CHECK(ntP->el_isGsfCtfScPixChargeConsistent) ) _isPCC     = ntP->el_isGsfCtfScPixChargeConsistent->at(idx);

    // mva-based selection additionnal variables
    // placeholder personnal MVA

    // more variables
    if( CHECK(ntP->el_ip3d)                          ) _ip3d      = ntP->el_ip3d->at(idx);
    if( CHECK(ntP->el_ip3dErr)                       ) _ip3dErr   = ntP->el_ip3dErr->at(idx);
}

void Electron::init()
{

    // general informations
    _E            = -666;
    _pt           = -666;
    _eta          = -666;
    _phi          = -666;
    _m            = -666;
    _charge       = 0;
    _id           = 0;

    // Id
    _isLoose      = -666;
    _isMedium     = -666;

    // variables for Id
    _dxy          = -666;
    _dz           = -666;

    // more variables

}

bool Electron::sel()
{

    float miniIso            = ntP->el_miniIsoTTH->at(idx);
    float SIP                = fabs(_ip3d/_ip3dErr);
    bool  isLoose            = false;

    //std::cout << "Ele Non Trig MVA discriminant = " << ntP->el_mvaNonTrigV0->at(idx) << std::endl;

    /*if (_pt <= 10)
    {
        if      (fabs(_eta) < 0.8  )  { isLoose = ( ntP->el_mvaNonTrigV0->at(idx) > -0.11 ); }
        else if (fabs(_eta) < 1.479)  { isLoose = ( ntP->el_mvaNonTrigV0->at(idx) > -0.55 ); }
        else                          { isLoose = ( ntP->el_mvaNonTrigV0->at(idx) > -0.60 ); }
    }
    else
    {
         if      (fabs(_eta) < 0.8  ) { isLoose = ( ntP->el_mvaNonTrigV0->at(idx) > -0.16 ); }
         else if (fabs(_eta) < 1.479) { isLoose = ( ntP->el_mvaNonTrigV0->at(idx) > -0.65 ); }
         else                         { isLoose = ( ntP->el_mvaNonTrigV0->at(idx) > -0.74 ); }    
    }*/

    if      (fabs(_eta) < 0.8  ) { isLoose = ( ntP->el_mvaNonTrigV0->at(idx) > -0.70 ); }
    else if (fabs(_eta) < 1.479) { isLoose = ( ntP->el_mvaNonTrigV0->at(idx) > -0.83 ); }
    else                         { isLoose = ( ntP->el_mvaNonTrigV0->at(idx) > -0.92 ); }  

    // preselection
    bool pass_pt       = (_pt        > 7   );
    bool pass_eta      = (fabs(_eta) < 2.5 );
    bool pass_dxy      = (fabs(_dxy) < 0.05);
    bool pass_dz       = (fabs(_dz)  < 0.1 );
    bool pass_miniIso  = (miniIso    < 0.4 );
    bool pass_SIP      = (SIP        < 8   );
    bool pass_isLoose  = (isLoose          );
    bool pass_losthits = (_nlosthits < 2   );
    bool pass_CV       = (_passCV          );

    bool pass_muOverlap = 1;
    int nMuon = nt->NtMuon->size();
    for(int im=0;im<nMuon;im++)
    {
        float dr = GetDeltaR(_eta,_phi,nt->NtMuon->at(im).eta(),nt->NtMuon->at(im).phi());
        if( dr < 0.05 ) pass_muOverlap = 0; //&& nt->NtMuon->at(im).isLoose() ) pass_muOverlap = 0;
    }

    bool isPreselectionElectron     = ( pass_pt        &&
                                        pass_eta       &&
                                        pass_dxy       &&
                                        pass_dz        &&
                                        pass_miniIso   &&
                                        pass_SIP       &&
                                        pass_isLoose   &&
                                        pass_losthits  &&
                                        pass_CV        &&
                                        pass_muOverlap );


    // mva-based selection
    //bool pass_mva      = (_mva       > 0.65);
    //bool pass_isMedium = (_isMedium        );

    //bool isSelectionElectron        = ( isPreselectionElectron &&
                                        //pass_mva               &&
                                        //pass_isMedium          &&
    //                                  )

    cout<<std::setiosflags(ios::fixed)<<setprecision(5);
    // synchronization printout
    //if( isPreselectionElectron ) 
    /*std::cout                                                     << std::setw(10)
                                           << nt->NtEvent->at(0).id()                          << std::setw(10)
                                           << _pt                                              << std::setw(10)
                                           << _eta                                             << std::setw(10)
                                           << _phi                                             << std::setw(10)
                                           << _E                                               << std::setw(5)
                                           << ntP->el_id->at(idx)                              << std::setw(5)
                                           << ntP->el_charge->at(idx)                          << std::setw(15)
                                           << ntP->el_miniIsoTTH->at(idx)                      << std::setw(15)
                                           << ntP->el_lepMVA_miniRelIsoCharged->at(idx)        << std::setw(15)
                                           << ntP->el_lepMVA_miniRelIsoNeutral->at(idx)        << std::setw(10)
                                           << ntP->el_lepMVA_jetPtRelv2->at(idx)               << std::setw(10)
                                           //<< 0.0                                              << std::setw(10)
                                           << ntP->el_lepMVA_jetBTagCSV->at(idx)               << std::setw(10)
                                           << ntP->el_lepMVA_jetPtRatio->at(idx)               << std::setw(10)
                                           << SIP                                              << std::setw(10)
                                           << fabs(_dxy)                                       << std::setw(10)
                                           << fabs(_dz)                                        << std::setw(21)
                                           << fabs( ntP->el_mvaNonTrigV0->at(idx) )            << std::endl;
     */
    
     /*std::cout << "Event:           " << nt->NtEvent->at(0).id()
               << " pass_pt:        " << pass_pt
               << " pass_eta:       " << pass_eta
               << " pass_dxy:       " << pass_dxy
               << " pass_dz:        " << pass_dz
               << " pass_miniIso:   " << pass_miniIso
               << " pass_SIP:       " << pass_SIP
               << " pass_isLoose:   " << pass_isLoose
               << " pass_losthits:  " << pass_losthits
               << " pass_CV:        " << pass_CV
               << " pass_muOverlap: " << pass_muOverlap << std::endl;*/

    return isPreselectionElectron;
}

float Electron::effectiveArea(int dr,float eta)
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
