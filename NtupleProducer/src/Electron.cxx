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
     _ID                                = idx;
    
     _E 	                            = ntP->el_E->at(idx);
     _pt	                            = ntP->el_pt->at(idx);
     _eta	                            = ntP->el_eta->at(idx);
     _phi	                            = ntP->el_phi->at(idx);
     _m 	                            = ntP->el_m->at(idx);
     _charge	                        = ntP->el_charge->at(idx);
     _id	                            = ntP->el_id->at(idx);
    
     _dxy	                            = ntP->el_gsfTrack_PV_dxy->at(idx);
     _dz	                            = ntP->el_gsfTrack_PV_dz->at(idx);
     _nlosthits                         = ntP->el_numberOfLostHits->at(idx);
     _ip3d	                            = ntP->el_ip3d->at(idx);
     _ip3dErr	                        = ntP->el_ip3dErr->at(idx);

     _isLoose	                        = ntP->el_looseCBId->at(idx);
     _isMedium                          = ntP->el_mediumCBId->at(idx);
     _passCV	                        = ntP->el_passConversionVeto->at(idx);
     _isPCC	                            = ntP->el_isGsfCtfScPixChargeConsistent->at(idx);

     _lepMVA	                        = ntP->el_lepMVA->at(idx);
     _lepMVA_TTH                        = ntP->el_lepMVA_Moriond16->at(idx);
     _mvaNonTrigV0                      = ntP->el_mvaNonTrigV0->at(idx);
 
     _lepMVA_miniRelIsoCharged          = ntP->el_lepMVA_miniRelIsoCharged->at(idx);
     _lepMVA_miniRelIsoNeutral          = ntP->el_lepMVA_miniRelIsoNeutral->at(idx);
     _lepMVA_jetPtRelv2                 = ntP->el_lepMVA_jetPtRelv2->at(idx);  
     //_lepMVA_jetDR                    = ntP->el_lepMVA_jetDR->at(idx);
     _lepMVA_jetPtRatio                 = ntP->el_lepMVA_jetPtRatio->at(idx);
     _lepMVA_jetBTagCSV                 = ntP->el_lepMVA_jetBTagCSV->at(idx);
     _lepMVA_sip3d                      = ntP->el_lepMVA_sip3d->at(idx);
     _lepMVA_dxy                        = ntP->el_lepMVA_dxy->at(idx);
     _lepMVA_dz                         = ntP->el_lepMVA_dz->at(idx);
     _lepMVA_mvaId                      = ntP->el_lepMVA_mvaId->at(idx);
     _lepMVA_eta                        = ntP->el_lepMVA_eta->at(idx);
     _lepMVA_jetNDauChargedMVASel       = ntP->el_lepMVA_jetNDauChargedMVASel->at(idx);

     _miniIso			                = ntP->el_miniIsoTTH->at(idx);
  
     _sigmaIetaIeta		                = ntP->el_sigmaIetaIeta->at(idx);    
     _superCluster_eta  	            = ntP->el_superCluster_eta->at(idx);     
     _hadronicOverEm		            = ntP->el_hadronicOverEm->at(idx);    
     _deltaEtaSuperClusterTrackAtVtx    = ntP->el_deltaEtaSuperClusterTrackAtVtx->at(idx);
     _deltaPhiSuperClusterTrackAtVtx    = ntP->el_deltaPhiSuperClusterTrackAtVtx->at(idx);
     _eSuperClusterOverP	            = ntP->el_eSuperClusterOverP->at(idx);	
     _correctedEcalEnergy	            = ntP->el_correctedEcalEnergy->at(idx);
     _ecalEnergy	                    = ntP->el_ecalEnergy->at(idx);
    
     _trackMomentumError	            = ntP->el_trackMomentumError->at(idx);
     _tightCharge                       = ntP->el_isGsfCtfScPixChargeConsistent->at(idx) + ntP->el_isGsfScPixChargeConsistent->at(idx);
}

void Electron::init()
{

    // general informations
    _E                              = 0.;
    _pt                             = 0.;
    _eta                            = -888;
    _phi                            = -888;
    _m                              = 0.;
    _charge                         = 0;
    _id                             = 0;

    // Id
    _isLooseCBId           	        = 0;
    _isMediumCBId          	        = 0;
    _isLoose               	        = 0;
    _isMedium                  	    = 0;
    _isTight                   	    = 0;
    _isLooseMVA                	    = 0;
    _isTightMVA             	    = 0;
    _isLooseTTH            	        = 0;
    _isFakeableTTH             	    = 0;
    _isTightTTH                	    = 0;


    // variables for Id
    _dxy                            = -888;
    _dz                             = -888; 
    _miniIso                        = -888;
    _nlosthits                      = -888;
   
    _passCV                         = 0;
    _isPCC                          = 0;
    _passPtEta                      = 0;
    _ip3d                           = -888.;
    _ip3dErr                        = -888.;

    // more variables
    _lepMVA                         = -888; 
    _lepMVA_TTH                     = -888;  
    
    
    _lepMVA_miniRelIsoCharged       = -888.;
    _lepMVA_miniRelIsoNeutral       = -888.;
    _lepMVA_jetPtRelv2              = -888.;  
    //_lepMVA_jetDR                   = -888.;
    _lepMVA_jetPtRatio              = -888.;
    _lepMVA_jetBTagCSV              = -888.;
    _lepMVA_sip3d                   = -888.;
    _lepMVA_dxy                     = -888.;
    _lepMVA_dz                      = -888.;
    _lepMVA_mvaId                   = -888.;
    
    _lepMVA_eta                     = -888.;
    _lepMVA_jetNDauChargedMVASel    = -888.;
    
    _passChargeFlip                 = 0;
    _hasMatchedConversion           = 0;
    _isGsfCtfScPixChargeConsistent  = 0;
   
    _sigmaIetaIeta                  = -888.;
    _hadronicOverEm                 = -888.;
    _correctedEcalEnergy            = -888.; 
    _ecalEnergy                     = -888.;
    _eSuperClusterOverP             = -888.;
    _deltaEtaSuperClusterTrackAtVtx = -888.;
    _deltaPhiSuperClusterTrackAtVtx = -888.;
    _see                            = -888.;
    _superCluster_eta               = -888.;
	
    _trackMomentumError             = -888;      
    _tightCharge                    = 0;
    _mvaNonTrigV0                   = -888;
}

bool Electron::sel()
{

    float SIP     = fabs(_ip3d/_ip3dErr);
    
    //MVA ID
    bool  isLoose = false;
    
    if      (fabs(_eta) < 0.8  ) { isLoose = ( _mvaNonTrigV0 > -0.70 ); }
    else if (fabs(_eta) < 1.479) { isLoose = ( _mvaNonTrigV0 > -0.83 ); }
    else                         { isLoose = ( _mvaNonTrigV0 > -0.92 ); }  


    // Loose
    bool pass_pt       = (_pt        > 7   );
    bool pass_eta      = (fabs(_eta) < 2.5 );
    bool pass_dxy      = (fabs(_dxy) < 0.05);
    bool pass_dz       = (fabs(_dz)  < 0.1 );
    bool pass_miniIso  = (_miniIso    < 0.4 );
    bool pass_SIP      = (SIP        < 8   );
    bool pass_isLoose  = (isLoose          );
    bool pass_losthits = (_nlosthits < 2   );
    
    //
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
        if (  ( _sigmaIetaIeta >= 0.011 && fabs(_superCluster_eta) <  1.479 )
           || ( _sigmaIetaIeta >= 0.030 && fabs(_superCluster_eta) >= 1.479 )  ) cond_closuretest = false;

        if (  ( ntP->el_hadronicOverEm->at(idx) >= 0.10 && fabs(_superCluster_eta) <  1.479 )
           || ( ntP->el_hadronicOverEm->at(idx) >= 0.07 && fabs(_superCluster_eta) >= 1.479 )  ) cond_closuretest = false;

        if (  ( fabs(_deltaEtaSuperClusterTrackAtVtx) >= 0.010 && fabs(_superCluster_eta) <  1.479 ) 
           || ( fabs(_deltaEtaSuperClusterTrackAtVtx) >= 0.008 && fabs(_superCluster_eta) >= 1.479 )  ) cond_closuretest = false;
        if (  ( fabs(_deltaPhiSuperClusterTrackAtVtx) >= 0.04  && fabs(_superCluster_eta) <  1.479 ) 
           || ( fabs(_deltaPhiSuperClusterTrackAtVtx) >= 0.07  && fabs(_superCluster_eta) >= 1.479 )  ) cond_closuretest = false;

        if (_ecalEnergy == 0.)
        {
            cond_closuretest = false;
        }
        else
        {
            if ( ( 1. / _ecalEnergy - _eSuperClusterOverP / _ecalEnergy ) < -0.005 )
            {cond_closuretest = false;}
            if ( ( 1. / _ecalEnergy - _eSuperClusterOverP / _ecalEnergy ) >= 0.010 && fabs(_superCluster_eta) <  1.479)
            {cond_closuretest = false;}
            if ( ( 1. / _ecalEnergy - _eSuperClusterOverP / _ecalEnergy ) >= 0.005 && fabs(_superCluster_eta) >= 1.479)
            {cond_closuretest = false;}
        }
    }

    // Fakeable

    pass_losthits = (_nlosthits == 0 );
    pass_pt       = (_pt        >  10); // should be 0.85 * pt(jet) for fakeable object cf v4 of note

    bool pass_lepMVA_TTH  = _lepMVA_TTH > 0.75 ;
    
    bool pass_lepMVA_jetBTagCSV089 = _lepMVA_jetBTagCSV < 0.89;
    
    bool pass_lepMVA_jetBtagCSVPtRatio = false;
    
    if (!pass_lepMVA_TTH && _lepMVA_jetPtRatio > 0.3 && _lepMVA_jetBTagCSV < 0.605) pass_lepMVA_jetBtagCSVPtRatio = true;
    if ( pass_lepMVA_TTH && pass_lepMVA_jetBTagCSV089) pass_lepMVA_jetBtagCSVPtRatio = true;
    
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

    bool pass_CV            = (_passCV         );
    bool pass_tightCharge   = (_tightCharge  >1);
  
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
			                pass_lepMVA_TTH       &&
			                pass_lepMVA_jetBTagCSV089 && 
			                pass_tightCharge);

    _isTightTTH = isTightTTH;

    cout<<std::setiosflags(ios::fixed)<<setprecision(5);
    
    // synchronization printout
    if( false ) //isLooseTTH )
    {
        std::cout   << nt->NtEvent->at(0).id()                          << " "
                    << _pt                                              << " "
                    << _eta                                             << " "
                    << _phi                                             << " "
                    << _E                                               << " "
                    << ntP->el_id->at(idx)                              << " "
                    << ntP->el_charge->at(idx)                          << " "
                    << ntP->el_miniIsoTTH->at(idx)                      << " "
                    << ntP->el_lepMVA_miniRelIsoCharged->at(idx)        << " "
                    << ntP->el_lepMVA_miniRelIsoNeutral->at(idx)        << " "
                    << ntP->el_lepMVA_jetPtRelv2->at(idx)               << " "
                    //<< 0.0                                              << " "
                    << ntP->el_lepMVA_jetBTagCSV->at(idx)               << " "
                    << ntP->el_lepMVA_jetPtRatio->at(idx)               << " "
                    << SIP                                              << " "
                    << fabs(_dxy)                                       << " "
                    << fabs(_dz)                                        << " "
                    << fabs( ntP->el_mvaNonTrigV0->at(idx) )            << " "
                    << _lepMVA_TTH                                      << " "
                    << std::endl;
    
        /*std::cout   << "Event:           "              << nt->NtEvent->at(0).id()
                    << " pass_pt:        "              << pass_pt
                    << " pass_eta:       "              << pass_eta
                    << " pass_dxy:       "              << pass_dxy
                    << " pass_dz:        "              << pass_dz
                    << " pass_miniIso:   "              << pass_miniIso
                    << " pass_SIP:       "              << pass_SIP
                    << " pass_isLoose:   "              << pass_isLoose
                    << " pass_losthits:  "              << pass_losthits
                    << " pass_CV:        "              << pass_CV
                    << " pass_muOverlap: "              << pass_muOverlap 
                    << " pass_lepMVA_TTH: "             << pass_lepMVA_TTH
                    << " pass_lepMVA_jetBTagCSV089: "   << pass_lepMVA_jetBTagCSV089
                    << " pass_tightCharge: "            << pass_tightCharge
                    << " is Tight: "                    << isTightTTH
                    << std::endl;*/
    }

    return isLooseTTH;
}
