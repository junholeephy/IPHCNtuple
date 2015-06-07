#include "include/NtupleProducer.h"

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
   
   if( CHECK(ntP->el_E) )                _E   = ntP->el_E->at(idx);
   if( CHECK(ntP->el_pt) )               _pt  = ntP->el_pt->at(idx);
   if( CHECK(ntP->el_eta) )              _eta = ntP->el_eta->at(idx);   
   if( CHECK(ntP->el_phi) )              _phi = ntP->el_phi->at(idx);   
   if( CHECK(ntP->el_m) )                _m   = ntP->el_m->at(idx);
   if( CHECK(ntP->el_dxy) )                _dxy   = ntP->el_dxy->at(idx);
   if( CHECK(ntP->el_dz) )                _dz   = ntP->el_dz->at(idx);
   
   if( CHECK(ntP->el_charge) )           _charge   = ntP->el_charge->at(idx);
   if( CHECK(ntP->el_id) )           _id   = ntP->el_id->at(idx);
   
   if( CHECK(ntP->el_lepMVA) )           _lepMVA   = ntP->el_lepMVA->at(idx);
   
   _passChargeFlip = ntP->el_isGsfCtfScPixChargeConsistent->at(idx);
   
   if( CHECK(ntP->el_lepMVA_neuRelIso) ) _lepMVA_neuRelIso = ntP->el_lepMVA_neuRelIso->at(idx);
   if( CHECK(ntP->el_lepMVA_chRelIso) ) _lepMVA_chRelIso = ntP->el_lepMVA_chRelIso->at(idx);
   if( CHECK(ntP->el_lepMVA_jetDR) ) _lepMVA_jetDR = ntP->el_lepMVA_jetDR->at(idx);
   if( CHECK(ntP->el_lepMVA_jetPtRatio) ) _lepMVA_jetPtRatio = ntP->el_lepMVA_jetPtRatio->at(idx);
   if( CHECK(ntP->el_lepMVA_jetBTagCSV) ) _lepMVA_jetBTagCSV = ntP->el_lepMVA_jetBTagCSV->at(idx);
   if( CHECK(ntP->el_lepMVA_sip3d) ) _lepMVA_sip3d = ntP->el_lepMVA_sip3d->at(idx);
   if( CHECK(ntP->el_lepMVA_dxy) ) _lepMVA_dxy = ntP->el_lepMVA_dxy->at(idx);
   if( CHECK(ntP->el_lepMVA_dz) ) _lepMVA_dz = ntP->el_lepMVA_dz->at(idx);
   if( CHECK(ntP->el_lepMVA_mvaId) ) _lepMVA_mvaId = ntP->el_lepMVA_mvaId->at(idx);

   _deltaEtaSuperClusterTrackAtVtx = ntP->el_deltaEtaSuperClusterTrackAtVtx->at(idx);
   _deltaPhiSuperClusterTrackAtVtx = ntP->el_deltaPhiSuperClusterTrackAtVtx->at(idx);
   _see = ntP->el_see->at(idx);
   _hadronicOverEm = ntP->el_hadronicOverEm->at(idx);
   _scleta = ntP->el_scleta->at(idx);
   
   _dB3D = ntP->el_dB3D->at(idx);
   _edB3D = ntP->el_edB3D->at(idx);
   
   _hasMatchedConversion = ntP->el_hasMatchedConversion->at(idx);
}

void Electron::init()
{
   _E        = -666;
   _pt       = -666;
   _eta      = -666;
   _phi      = -666;
   _m        = -666;
   _dxy        = -666;
   _dz        = -666;
   _iso        = -666;
   _isLoose   = 0;
   _isTight   = 0;
   _isLooseMVA   = 0;
   _isTightMVA   = 0;
   
   _charge   = 0;
   _id   = 0;
   
   _fakeType = -1;
//   _chargeFlip = 0;
   
   _lepMVA   = -666;
   _passChargeFlip = 0;
   _hasMatchedConversion = 0;
   
   _lepMVA_neuRelIso = -666;
   _lepMVA_chRelIso = -666;
   _lepMVA_jetDR = -666;
   _lepMVA_jetPtRatio = -666;
   _lepMVA_jetBTagCSV = -666;
   _lepMVA_sip3d = -666;
   _lepMVA_dxy = -666;
   _lepMVA_dz = -666;
   _lepMVA_mvaId = -666;

   _deltaEtaSuperClusterTrackAtVtx = -666;
   _deltaPhiSuperClusterTrackAtVtx = -666;
   _see = -666;
   _hadronicOverEm = -666;
   _scleta = -666;
   
   _dB3D = -666;
   _edB3D = -666;
}

void Electron::sel()
{   
   bool pass_mvaId_loose = 0;
   bool pass_mvaId_tight = 0;
   
   if( _pt > 10. )
     {	
	if( fabs(_eta) < 0.8 ) pass_mvaId_loose = (ntP->el_mvaNonTrigV0->at(idx) > 0.35);
	else if( fabs(_eta) >= 0.8 && fabs(_eta) < 1.479 ) pass_mvaId_loose = (ntP->el_mvaNonTrigV0->at(idx) > 0.20);
	else pass_mvaId_loose = (ntP->el_mvaNonTrigV0->at(idx) > -0.52);

	if( fabs(_eta) < 0.8 ) pass_mvaId_tight = (ntP->el_mvaNonTrigV0->at(idx) > 0.73);
	else if( fabs(_eta) >= 0.8 && fabs(_eta) < 1.479 ) pass_mvaId_tight = (ntP->el_mvaNonTrigV0->at(idx) > 0.57);
	else pass_mvaId_tight = (ntP->el_mvaNonTrigV0->at(idx) > 0.05);
     }
	
   bool pass_pt = (_pt > 7.);
   bool pass_eta = (fabs(_eta) < 2.5);
   bool pass_dxy = (fabs(_dxy) < 0.05);
   bool pass_dz = (fabs(_dz) < 0.1);
   bool pass_numberOfHits = (ntP->el_numberOfHits->at(idx) <= 1);
   bool pass_numberOfHits_tight = (ntP->el_numberOfHits->at(idx) == 0);

   _pass_isGsfCtfScPixChargeConsistent = (ntP->el_isGsfCtfScPixChargeConsistent->at(idx));
   
   float sumChargedHadronPt = ntP->el_pfIso_sumChargedHadronPt->at(idx);
   float sumNeutralHadronEt = ntP->el_pfIso_sumNeutralHadronEt->at(idx);
   float sumPhotonEt = ntP->el_pfIso_sumPhotonEt->at(idx);
   float sumPUPt = ntP->el_pfIso_sumPUPt->at(idx);

   float rho = nt->NtEvent->at(0).rho();
   float effArea = effectiveArea(30,_eta);
   
   _iso = (sumChargedHadronPt + std::max(sumNeutralHadronEt+sumPhotonEt-rho*effArea,float(0.)))/_pt;
   bool pass_iso_loose = (_iso < 0.5);
   bool pass_iso_tight = (_iso < 0.1);
   
   bool pass_muOverlap = 1;
   int nMuon = nt->NtMuon->size();
   for(int im=0;im<nMuon;im++)
     {
	float dr = GetDeltaR(_eta,_phi,nt->NtMuon->at(im).eta(),nt->NtMuon->at(im).phi());
	if( dr < 0.05 && nt->NtMuon->at(im).isTight() ) pass_muOverlap = 0;
     }  

   bool pass_sip3d = ( _edB3D > 0. ) ? (fabs(_dB3D/_edB3D) < 4.) : 0;

   _passPtEta = (pass_pt && pass_eta);
   
   _isLoose = (
	       pass_dxy &&
	       pass_dz &&
	       pass_numberOfHits &&
	       pass_iso_loose &&
	       pass_muOverlap &&
	       pass_mvaId_loose
	      );

   _isTight = (
	       _isLoose &&
	       pass_sip3d &&
	       pass_mvaId_tight &&
	       pass_numberOfHits_tight &&
	       !_hasMatchedConversion &&
	       pass_iso_tight
	      );
   
   _isLooseMVA = (_lepMVA > 0.5);
   _isTightMVA = (_lepMVA > 0.8);

   for(int id=0;id<evdebug->size();id++)
     {	
	if( nt->NtEvent->at(0).id() == evdebug->at(id) )
	  {
	     std::cout << "Electron #" << _ID << std::endl;
	     std::cout << "   pt=" << _pt << " eta=" << _eta << " phi=" << _phi << std::endl;
	     std::cout << "   passPtEta=" << _passPtEta << std::endl;
	     std::cout << "   isLoose=" << _isLoose << std::endl;
	     std::cout << "   isTight=" << _isTight << std::endl;
	     std::cout << "   isLooseMVA=" << _isLooseMVA << std::endl;
	     std::cout << "   isTightMVA=" << _isTightMVA << std::endl;
	     std::cout << "   mva=" << ntP->el_mvaNonTrigV0->at(idx) << std::endl;
	     std::cout << "   lepMVA=" << _lepMVA << std::endl;
	     std::cout << "   iso=" << _iso << std::endl;
	     std::cout << "   pass_mvaId_loose=" << pass_mvaId_loose << std::endl;
	     std::cout << "   pass_mvaId_tight=" << pass_mvaId_tight << std::endl;
	     std::cout << "   pass_iso_tight=" << pass_iso_tight << std::endl;
	     std::cout << "   pass_muOverlap=" << pass_muOverlap << std::endl;
	     std::cout << "   pass_sip3d=" << pass_sip3d << std::endl;
	     std::cout << "   pass_numberOfHits_tight=" << pass_numberOfHits_tight << std::endl;
	     std::cout << "   _hasMatchedConversion=" << !_hasMatchedConversion << std::endl;
	     std::cout << "   lepMVA_neuRelIso=" << _lepMVA_neuRelIso << std::endl;
	     std::cout << "   lepMVA_chRelIso=" << _lepMVA_chRelIso << std::endl;
	     std::cout << "   lepMVA_jetDR=" << _lepMVA_jetDR << std::endl;
	     std::cout << "   lepMVA_jetPtRatio=" << _lepMVA_jetPtRatio << std::endl;
	     std::cout << "   lepMVA_jetBTagCSV=" << _lepMVA_jetBTagCSV << std::endl;
	     std::cout << "   lepMVA_sip3d=" << _lepMVA_sip3d << std::endl;
	     std::cout << "   lepMVA_dxy=" << _lepMVA_dxy << std::endl;
	     std::cout << "   lepMVA_dz=" << _lepMVA_dz << std::endl;
	     std::cout << "   lepMVA_mvaId=" << _lepMVA_mvaId << std::endl;
	  }
     }
}

float Electron::effectiveArea(int dr,float eta)
{
   float effArea = 0.;
   
   if( fabs(eta) >= 0 && fabs(eta) < 0.8 )
     {
	if( dr == 30 ) effArea = 0.1013;
	else if( dr == 40 ) effArea = 0.1830;
     }   
   else if( fabs(eta) >= 0.8 && fabs(eta) < 1.3 )
     {
	if( dr == 30 ) effArea = 0.0988;
	else if( dr == 40 ) effArea = 0.1734;
     }   
   else if( fabs(eta) >= 1.3 && fabs(eta) < 2.0 )
     {
	if( dr == 30 ) effArea = 0.0572;
	else if( dr == 40 ) effArea = 0.1077;
     }   
   else if( fabs(eta) >= 2.0 && fabs(eta) < 2.2 )
     {
	if( dr == 30 ) effArea = 0.0842;
	else if( dr == 40 ) effArea = 0.1565;
     }   
   else if( fabs(eta) >= 2.2 && fabs(eta) < 2.5 )
     {
	if( dr == 30 ) effArea = 0.1530;
	else if( dr == 40 ) effArea = 0.2680;
     }   
   
   return effArea;
}
