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
    if( CHECK(ntP->el_charge)                        ) _charge    = ntP->el_charge->at(idx);
    if( CHECK(ntP->el_id)                            ) _id        = ntP->el_id->at(idx);

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
    if( CHECK(ntP->el_lepMVA)                        ) _lepMVA           = ntP->el_lepMVA->at(idx);
    if( CHECK(ntP->el_lepMVA_Moriond16)              ) _lepMVA_Moriond16 = ntP->el_lepMVA_Moriond16->at(idx);
   
    if( CHECK(ntP->el_lepMVA_neuRelIso)	             ) _lepMVA_neuRelIso = ntP->el_lepMVA_neuRelIso->at(idx);
    if( CHECK(ntP->el_lepMVA_chRelIso) 	             ) _lepMVA_chRelIso = ntP->el_lepMVA_chRelIso->at(idx);
    if( CHECK(ntP->el_lepMVA_jetDR)		             ) _lepMVA_jetDR = ntP->el_lepMVA_jetDR->at(idx);
    if( CHECK(ntP->el_lepMVA_jetPtRatio) 	         ) _lepMVA_jetPtRatio = ntP->el_lepMVA_jetPtRatio->at(idx);
    if( CHECK(ntP->el_lepMVA_jetBTagCSV)	         ) _lepMVA_jetBTagCSV = ntP->el_lepMVA_jetBTagCSV->at(idx);
    if( CHECK(ntP->el_lepMVA_sip3d) 	             ) _lepMVA_sip3d = ntP->el_lepMVA_sip3d->at(idx);
    if( CHECK(ntP->el_lepMVA_dxy)		             ) _lepMVA_dxy = ntP->el_lepMVA_dxy->at(idx);
    if( CHECK(ntP->el_lepMVA_dz) 	                 ) _lepMVA_dz = ntP->el_lepMVA_dz->at(idx);
    if( CHECK(ntP->el_lepMVA_mvaId)		             ) _lepMVA_mvaId = ntP->el_lepMVA_mvaId->at(idx);
    if( CHECK(ntP->el_lepMVA_eta) 	                 ) _lepMVA_eta = ntP->el_lepMVA_eta->at(idx);
    if( CHECK(ntP->el_lepMVA_jetNDauChargedMVASel)   ) _lepMVA_jetNDauChargedMVASel = ntP->el_lepMVA_jetNDauChargedMVASel->at(idx);
 
    // more variables
    if( CHECK(ntP->el_ip3d)                          ) _ip3d      = ntP->el_ip3d->at(idx);
    if( CHECK(ntP->el_ip3dErr)                       ) _ip3dErr   = ntP->el_ip3dErr->at(idx);

    // fake rate related closure test variables
    //if( CHECK(                                          _full5x5_sigmaIetaIeta
    //if( CHECK(                                          _superCluster_eta
    //if( CHECK(                                          _hadronicOverEm
    //if( CHECK(                                          _deltaEtaSuperClusterTrackAtVtx
    //if( CHECK(                                          _deltaPhiSuperClusterTrackAtVtx
    //if( CHECK(                                          _ecalEnergy
    //if( CHECK(                                          _eSuperClusterOverP

}

void Electron::init()
{

    // general informations
    _E                              = -888;
    _pt                             = -888;
    _eta                            = -888;
    _phi                            = -888;
    _m                              = -888;
    _charge                         = 0;
    _id                             = 0;

    // Id
    _isLoose                        = -888;
    _isMedium                       = -888;
    _isLooseTTH                     = -888;
    _isFakeableTTH                  = -888;
    _isTightTTH                     = -888;


    // variables for Id
    _dxy                            = -888;
    _dz                             = -888;
    _passCV                         = -888;
    _isPCC                          = -888;

    // more variables
    _lepMVA                         = -888; 
    _lepMVA_Moriond16               = -888;  
    
    _lepMVA_neuRelIso               = -888;
    _lepMVA_chRelIso                = -888;
    _lepMVA_jetDR                   = -888;
    _lepMVA_jetPtRatio              = -888;
    _lepMVA_jetBTagCSV              = -888;
    _lepMVA_sip3d                   = -888;
    _lepMVA_dxy                     = -888;
    _lepMVA_dz                      = -888;
    _lepMVA_mvaId                   = -888;
    _lepMVA_eta                     = -888;
    _lepMVA_jetNDauChargedMVASel    = -888;

    //_full5x5_sigmaIetaIeta          = -1000;
    //_superCluster_eta               = -1000;
    //_hadronicOverEm                 = -1000;
    //_deltaEtaSuperClusterTrackAtVtx = -1000;
    //_deltaPhiSuperClusterTrackAtVtx = -1000;
    //_ecalEnergy                     = -1000;
    //_eSuperClusterOverP             = -1000;
}

bool Electron::sel()
{

    float miniIso            = ntP->el_miniIsoTTH->at(idx);
    float SIP                = fabs(_ip3d/_ip3dErr);
    bool  isLoose            = false;

    if      (fabs(_eta) < 0.8  ) { isLoose = ( ntP->el_mvaNonTrigV0->at(idx) > -0.70 ); }
    else if (fabs(_eta) < 1.479) { isLoose = ( ntP->el_mvaNonTrigV0->at(idx) > -0.83 ); }
    else                         { isLoose = ( ntP->el_mvaNonTrigV0->at(idx) > -0.92 ); }  

    // Loose
    bool pass_pt       = (_pt        > 7   );
    bool pass_eta      = (fabs(_eta) < 2.5 );
    bool pass_dxy      = (fabs(_dxy) < 0.05);
    bool pass_dz       = (fabs(_dz)  < 0.1 );
    bool pass_miniIso  = (miniIso    < 0.4 );
    bool pass_SIP      = (SIP        < 8   );
    bool pass_isLoose  = (isLoose          );
    bool pass_losthits = (_nlosthits < 2   );

    bool pass_muOverlap = 1;
    int nMuon = nt->NtMuon->size();
    for(int im=0;im<nMuon;im++)
    {
        float dr = GetDeltaR(_eta,_phi,nt->NtMuon->at(im).eta(),nt->NtMuon->at(im).phi());
        if( dr < 0.05 ) pass_muOverlap = 0; //&& nt->NtMuon->at(im).isLoose() ) pass_muOverlap = 0;
    }

    bool isLooseTTH     = ( pass_pt          &&
                            pass_eta         &&
                            pass_dxy         &&
                            pass_dz          &&
                            pass_miniIso     &&
                            pass_SIP         &&
                            pass_isLoose     &&
                            pass_losthits    &&
                            pass_muOverlap   );

    _isLooseTTH = isLooseTTH; // OK

    // preselection electron for pt > 30 (aiming at closure for the fake rate method whatever)
    // from https://github.com/peruzzim/cmg-cmssw/blob/works_260116/CMGTools/TTHAnalysis/python/tools/emulateElectronTriggerCuts.py#L1-L8
    // def _susy2lss_idEmu_cuts(lep):
    //if (abs(lep.pdgId())!=11): return True OK
    //if (lep.full5x5_sigmaIetaIeta()>=(0.011 if abs(lep.superCluster().eta())<1.479 else 0.031)): return False
    //if (lep.hadronicOverEm()>=0.08): return False
    //if (abs(lep.deltaEtaSuperClusterTrackAtVtx())>=0.01): return False
    //if (abs(lep.deltaPhiSuperClusterTrackAtVtx())>=(0.04 if abs(lep.superCluster().eta())<1.479 else 0.08)): return False
    //if (abs((1.0/lep.ecalEnergy() - lep.eSuperClusterOverP()/lep.ecalEnergy()) if lep.ecalEnergy()>0. else 9e9)>=0.01): return False
    //return True

    bool cond_closuretest = true;
    if ( _pt > 30 )
    {
        if (  ( (ntP->el_sigmaIetaIeta->at(idx) >= 0.011) && (abs(ntP->el_superCluster_eta->at(idx)) <  1.479) )
                || ( (ntP->el_sigmaIetaIeta->at(idx) >= 0.031) && (abs(ntP->el_superCluster_eta->at(idx)) >= 1.479) )                                       ) cond_closuretest = false;
        if (  ntP->el_hadronicOverEm->at(idx)                                                                                                 >= 0.08  ) cond_closuretest = false;
        if (  abs( ntP->el_deltaEtaSuperClusterTrackAtVtx->at(idx) )                                                                          >= 0.01  ) cond_closuretest = false;
        if (  ( (abs( ntP->el_deltaPhiSuperClusterTrackAtVtx->at(idx)) >= 0.04) && (abs(ntP->el_superCluster_eta->at(idx) <  1.479) ) )
                || ( (abs( ntP->el_deltaPhiSuperClusterTrackAtVtx->at(idx)) >= 0.08) && (abs(ntP->el_superCluster_eta->at(idx) >= 1.479) ) )                ) cond_closuretest = false;
        if (ntP->el_correctedEcalEnergy->at(idx) == 0.)
        {
            cond_closuretest = false;
        }
        else
        {
            if( abs( 1. / ntP->el_correctedEcalEnergy->at(idx) - ntP->el_eSuperClusterOverP->at(idx) / ntP->el_correctedEcalEnergy->at(idx) ) >= 0.01)
            {cond_closuretest = false;}
        }
    }

    // Fakeable

    pass_losthits = (_nlosthits == 0 );
    pass_pt       = (_pt        >  10);

    bool pass_lepMVA_Moriond16  = _lepMVA_Moriond16 > 0.75 ;
    
    bool pass_lepMVA_jetBTagCSV089 = _lepMVA_jetBTagCSV < 0.89;
    
    bool pass_lepMVA_jetBtagCSVPtRatio = false;
    
    if (!pass_lepMVA_Moriond16 && _lepMVA_jetPtRatio > 0.3 && _lepMVA_jetBTagCSV < 0.605) pass_lepMVA_jetBtagCSVPtRatio = true;
    if (pass_lepMVA_Moriond16 && pass_lepMVA_jetBTagCSV089) pass_lepMVA_jetBtagCSVPtRatio = true;
    
    bool isFakeableTTH     = ( pass_pt          &&
                               pass_eta         &&
                               pass_dxy         &&
                               pass_dz          &&
                               pass_miniIso     &&
                               pass_SIP         &&
                               pass_isLoose     &&
                               pass_losthits    &&
                               cond_closuretest &&
                               pass_muOverlap   && 
			       pass_lepMVA_jetBtagCSVPtRatio
			       );
    
    _isFakeableTTH = isFakeableTTH;

    // Tight

    bool pass_CV       = (_passCV          );

    bool isTightTTH     = ( pass_pt               &&
                            pass_eta              &&
                            pass_dxy              &&
                            pass_dz               &&
                            pass_miniIso          &&
                            pass_SIP              &&
                            pass_isLoose          &&
                            pass_losthits         &&
                            pass_CV               &&
                            cond_closuretest      &&
                            pass_muOverlap        && 
			    pass_lepMVA_Moriond16 &&
			    pass_lepMVA_jetBTagCSV089);

    _isTightTTH = isTightTTH;

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

    return isLooseTTH;
}
