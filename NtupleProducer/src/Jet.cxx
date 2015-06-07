#include "include/NtupleProducer.h"

ClassImp(Jet)
    
Jet::Jet()
{
}

Jet::~Jet()
{
}

void Jet::read()
{
   _ID = idx;
   
   if( CHECK(ntP->jet_E) )                _E   = ntP->jet_E->at(idx);
   if( CHECK(ntP->jet_pt) )               _pt  = ntP->jet_pt->at(idx);
   if( CHECK(ntP->jet_eta) )              _eta = ntP->jet_eta->at(idx);   
   if( CHECK(ntP->jet_phi) )              _phi = ntP->jet_phi->at(idx);   
   if( CHECK(ntP->jet_m) )                _m   = ntP->jet_m->at(idx);
   
   if( CHECK(ntP->jet_ntrk) )             _ntrk   = ntP->jet_ntrk->at(idx);
   
   if( CHECK(ntP->jet_CSV) )                _CSV   = ntP->jet_CSV->at(idx);
   if( CHECK(ntP->jet_CSVv2) )                _CSVv2   = ntP->jet_CSVv2->at(idx);
}

void Jet::init()
{
   _E        = -666;
   _pt       = -666;
   _eta      = -666;
   _phi      = -666;
   _m        = -666;
   
   _ntrk        = -666;
   
   _CSV        = -666.;
   _CSVv2        = -666.;

   _isTight   = 0;
}

void Jet::sel()
{   
   bool pass_pt = (_pt > 25.);
   bool pass_eta = (fabs(_eta) < 2.4);
   bool pass_jetId = 0;
   
   float pileupJetId = ntP->jet_pileupJetId->at(idx);
   
   if( fabs(_eta) >= 0 && fabs(_eta) < 2.5 )
     if( pileupJetId > -0.63 ) pass_jetId = 1;
   else if( fabs(_eta) >= 2.5 && fabs(_eta) < 2.75 )
     if( pileupJetId > -0.60 ) pass_jetId = 1;
   else if( fabs(_eta) >= 2.75 && fabs(_eta) < 3.0 )
     if( pileupJetId > -0.55 ) pass_jetId = 1;
   else if( fabs(_eta) >= 3.0 && fabs(_eta) < 5.2 )
     if( pileupJetId > -0.45 ) pass_jetId = 1;

   bool pass_muOverlap = 1;
   int nMuon = nt->NtMuon->size();
   for(int im=0;im<nMuon;im++)
     {
	float dr = GetDeltaR(_eta,_phi,nt->NtMuon->at(im).eta(),nt->NtMuon->at(im).phi());
	if( dr < 0.4 && nt->NtMuon->at(im).pt() > 10. && nt->NtMuon->at(im).isTight() ) pass_muOverlap = 0;
     }  

   bool pass_elOverlap = 1;
   int nElectron = nt->NtElectron->size();
   for(int ie=0;ie<nElectron;ie++)
     {
	float dr = GetDeltaR(_eta,_phi,nt->NtElectron->at(ie).eta(),nt->NtElectron->at(ie).phi());
	if( dr < 0.4 && nt->NtElectron->at(ie).pt() > 10. && nt->NtElectron->at(ie).isTight() ) pass_elOverlap = 0;
/*	for(int id=0;id<evdebug->size();id++)
	  {	
	     if( nt->NtEvent->at(0).id() == evdebug->at(id) )
	       {
		  std::cout << "el(pt=" << nt->NtElectron->at(ie).pt() << ") dr=" << dr << " overlap=" << pass_elOverlap << std::endl;
	       }
	  }*/
     }  
   
   _isTight = (
//	       pass_pt &&
//	       pass_eta &&
	       pass_jetId &&
	       pass_muOverlap &&
	       pass_elOverlap
	      );

   for(int id=0;id<evdebug->size();id++)
     {	
	if( nt->NtEvent->at(0).id() == evdebug->at(id) )
	  {
	     if( _isTight && _pt > 25. &&
		 fabs(_eta) < 2.4 )
	       {		  
		  std::cout << "Jet #" << _ID << std::endl;
		  std::cout << "   pt=" << _pt << " eta=" << _eta << " phi=" << _phi << std::endl;
		  std::cout << "   pileupJetId=" << pileupJetId << std::endl;
		  std::cout << "   isTight=" << _isTight << std::endl;
		  std::cout << "   pass_jetId=" << pass_jetId << std::endl;
		  std::cout << "   pass_elOverlap=" << pass_elOverlap << std::endl;
		  std::cout << "   pass_muOverlap=" << pass_muOverlap << std::endl;
	       }	     
	  }	
     }
}
