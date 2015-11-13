#ifndef LEPTON_H
#define LEPTON_H

#include "Electron.h"
#include "Muon.h"

#include "TLorentzVector.h"

class Lepton
{
   
 public:
   
   Lepton();
   virtual ~Lepton();

   float pt() {return _pt;};
   float eta() {return _eta;};
   float phi() {return _phi;};
   float E() {return _E;};
   
   TLorentzVector p4() {return _p4;};
   
   int id() {return _id;};
   
   int idx() {return _idx;};
   bool isElectron() {return _isElectron;};
   
   int charge() {return _charge;};
   
   template <class T> void setLepton(T *lep, int idx, bool isE)
     {
	_pt = lep->pt();
	_eta = lep->eta();
	_phi = lep->phi();
	_E = lep->E();
	
	_p4.SetPtEtaPhiE(_pt,_eta,_phi,_E);
	
	_id = lep->id();
	
	_idx = idx;
	_isElectron = isE;
	_charge = lep->charge();
     }
   
 protected:

   float _pt;
   float _eta;
   float _phi;
   float _E;
   
   int _id;
   
   TLorentzVector _p4;
   
   int _idx;
   bool _isElectron;
   
   int _charge;
};

#endif
