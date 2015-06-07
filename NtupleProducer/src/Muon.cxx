#include "include/NtupleProducer.h"

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
   
   if( CHECK(ntP->mu_E) )                _E   = ntP->mu_E->at(idx);
   if( CHECK(ntP->mu_pt) )               _pt  = ntP->mu_pt->at(idx);
   if( CHECK(ntP->mu_eta) )              _eta = ntP->mu_eta->at(idx);   
   if( CHECK(ntP->mu_phi) )              _phi = ntP->mu_phi->at(idx);   
   if( CHECK(ntP->mu_m) )                _m   = ntP->mu_m->at(idx);
   if( CHECK(ntP->mu_innerTrack_dxy) )                _dxy   = ntP->mu_innerTrack_dxy->at(idx);
   if( CHECK(ntP->mu_innerTrack_dz) )                _dz   = ntP->mu_innerTrack_dz->at(idx);
   
   if( CHECK(ntP->mu_charge) )                _charge   = ntP->mu_charge->at(idx);
   if( CHECK(ntP->mu_id) )                _id   = ntP->mu_id->at(idx);
   
   if( CHECK(ntP->mu_lepMVA) )                _lepMVA   = ntP->mu_lepMVA->at(idx);  

   _hasInnerTrack = ntP->mu_hasInnerTrack->at(idx);
   
   _isTightMuonOld = ntP->mu_isTightMuon->at(idx);
   bool goodGlb = (ntP->mu_isGlobalMuon->at(idx) &&
		   ntP->mu_globalTrack_normalizedChi2->at(idx) < 3 &&
		   ntP->mu_combinedQuality_chi2LocalPosition->at(idx) < 12 &&
		   ntP->mu_combinedQuality_trkKink->at(idx) < 20);
   _isTightMuon = ntP->mu_innerTrack_validFraction->at(idx) >= 0.8 &&
     ntP->mu_segmentCompatibility->at(idx) >= (goodGlb ? 0.303 : 0.451);
   
   _dB3D = ntP->mu_dB3D->at(idx);
   _edB3D = ntP->mu_edB3D->at(idx);
   
   if( ntP->mu_hasInnerTrack->at(idx) )
     {
	double sig = (ntP->mu_innerTrack_pt->at(idx) != 0.) ? ntP->mu_innerTrack_ptError->at(idx)/ntP->mu_innerTrack_pt->at(idx) : 666.;
	_passChargeFlip = (sig < 0.2);
     }
   
   if( CHECK(ntP->mu_lepMVA_neuRelIso) )                _lepMVA_neuRelIso   = ntP->mu_lepMVA_neuRelIso->at(idx);
   if( CHECK(ntP->mu_lepMVA_chRelIso) )                _lepMVA_chRelIso   = ntP->mu_lepMVA_chRelIso->at(idx);
   if( CHECK(ntP->mu_lepMVA_jetDR) )                _lepMVA_jetDR   = ntP->mu_lepMVA_jetDR->at(idx);
   if( CHECK(ntP->mu_lepMVA_jetPtRatio) )                _lepMVA_jetPtRatio   = ntP->mu_lepMVA_jetPtRatio->at(idx);
   if( CHECK(ntP->mu_lepMVA_jetBTagCSV) )                _lepMVA_jetBTagCSV   = ntP->mu_lepMVA_jetBTagCSV->at(idx);
   if( CHECK(ntP->mu_lepMVA_sip3d) )                _lepMVA_sip3d   = ntP->mu_lepMVA_sip3d->at(idx);
   if( CHECK(ntP->mu_lepMVA_dxy) )                _lepMVA_dxy   = ntP->mu_lepMVA_dxy->at(idx);
   if( CHECK(ntP->mu_lepMVA_dz) )                _lepMVA_dz   = ntP->mu_lepMVA_dz->at(idx);
   if( CHECK(ntP->mu_lepMVA_mvaId) )                _lepMVA_mvaId   = ntP->mu_lepMVA_mvaId->at(idx);
}

void Muon::init()
{
   _E        = -666;
   _pt       = -666;
   _eta      = -666;
   _phi      = -666;
   _m        = -666;
   _dxy        = -666;
   _dz        = -666;
   _iso        = -666;
   _isLoose        = 0;
   _isTight        = 0;
   _isLooseMVA        = 0;
   _isTightMVA        = 0;
   
   _isTightMuonOld        = 0;
   _isTightMuon        = 0;
   
   _charge        = 0;
   _id        = 0;
   
   _fakeType = -1;
//   _chargeFlip = 0;
   
   _lepMVA = -666;
   _passChargeFlip = 0;

   _hasInnerTrack = 0;
   
   _lepMVA_neuRelIso = -666;
   _lepMVA_chRelIso = -666;
   _lepMVA_jetDR = -666;
   _lepMVA_jetPtRatio = -666;
   _lepMVA_jetBTagCSV = -666;
   _lepMVA_sip3d = -666;
   _lepMVA_dxy = -666;
   _lepMVA_dz = -666;
   _lepMVA_mvaId = -666;
   
   _dB3D = -666;
   _edB3D = -666;
}

void Muon::sel()
{
   bool pass_loose = (ntP->mu_isPFMuon->at(idx) && 
		      (ntP->mu_isGlobalMuon->at(idx) ||
			  ntP->mu_isTrackerMuon->at(idx)));
   bool pass_pt = (_pt > 5.);
   bool pass_eta = (fabs(_eta) < 2.4);
   bool pass_dxy = (fabs(_dxy) < 0.05);
   bool pass_dz = (fabs(_dz) < 0.1);

   float sumChargedHadronPt = ntP->mu_pfIso03_sumChargedHadronPt->at(idx);
   float sumNeutralHadronEt = ntP->mu_pfIso03_sumNeutralHadronEt->at(idx);
   float sumPhotonEt = ntP->mu_pfIso03_sumPhotonEt->at(idx);
   float sumPUPt = ntP->mu_pfIso03_sumPUPt->at(idx);

   float rho = nt->NtEvent->at(0).rho();
   float effArea = effectiveArea(30,_eta);
   
   _iso = (sumChargedHadronPt + std::max(sumNeutralHadronEt+sumPhotonEt-rho*effArea,float(0.)))/_pt;
   bool pass_iso_loose = (_iso < 0.5);
   bool pass_iso_tight = (_iso < 0.1);

   bool pass_sip3d = ( _edB3D > 0. ) ? (fabs(_dB3D/_edB3D) < 4.) : 0;

   _passPtEta = (pass_pt && pass_eta);
   
   _isLoose = (
	       pass_loose &&
	       pass_dxy &&
	       pass_dz &&
	       pass_iso_loose
	      );
   
   _isTight = (
	       _isLoose &&
	       _isTightMuon &&
	       pass_sip3d &&
	       pass_iso_tight
	      );

   _isLooseMVA = (_lepMVA > 0.5);
   _isTightMVA = (_lepMVA > 0.8);
   
   for(int id=0;id<evdebug->size();id++)
     {	
	if( nt->NtEvent->at(0).id() == evdebug->at(id) )
	  {
	     std::cout << "Muon #" << _ID << std::endl;
	     std::cout << "   pt=" << _pt << " eta=" << _eta << " phi=" << _phi << std::endl;
	     std::cout << "   passPtEta=" << _passPtEta << std::endl;
	     std::cout << "   isLoose=" << _isLoose << std::endl;
	     std::cout << "   isTight=" << _isTight << std::endl;
	     std::cout << "   iso=" << _iso << std::endl;
	     std::cout << "   rho=" << rho << std::endl;
	     std::cout << "   effArea=" << effArea << std::endl;
	     std::cout << "   pass_loose=" << pass_loose << std::endl;
	     std::cout << "   pass_iso_tight=" << pass_iso_tight << std::endl;
	     std::cout << "   pass_dxy=" << pass_dxy << std::endl;
	     std::cout << "   pass_dz=" << pass_dz << std::endl;
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

float Muon::effectiveArea(int dr,float eta)
{
   float effArea = 0.;
   
   if( fabs(eta) >= 0 && fabs(eta) < 0.8 )
     {
	if( dr == 30 ) effArea = 0.0913;
	else if( dr == 40 ) effArea = 0.1564;
     }   
   else if( fabs(eta) >= 0.8 && fabs(eta) < 1.3 )
     {
	if( dr == 30 ) effArea = 0.0765;
	else if( dr == 40 ) effArea = 0.1325;
     }   
   else if( fabs(eta) >= 1.3 && fabs(eta) < 2.0 )
     {
	if( dr == 30 ) effArea = 0.0546;
	else if( dr == 40 ) effArea = 0.0913;
     }   
   else if( fabs(eta) >= 2.0 && fabs(eta) < 2.2 )
     {
	if( dr == 30 ) effArea = 0.0728;
	else if( dr == 40 ) effArea = 0.1212;
     }   
   else if( fabs(eta) >= 2.2 && fabs(eta) < 2.5 )
     {
	if( dr == 30 ) effArea = 0.1177;
	else if( dr == 40 ) effArea = 0.2085;
     }   
   
   return effArea;
}
