#include "../include/TransferFunc.h"
//#include "../include/RangesTF.h"

#include <boost/bind.hpp>                                                                                                   
#include <boost/function.hpp>                                                                                               
#include <boost/type_traits.hpp>

TransferFunc::TransferFunc()
{
   _v_ElectronTight = new std::vector<Electron>;
   _v_MuonTight = new std::vector<Muon>;
   _v_JetTight = new std::vector<Jet>;
}

TransferFunc::~TransferFunc()
{
   delete _v_ElectronTight;
   delete _v_MuonTight;
   delete _v_JetTight;
}

void TransferFunc::init()
{
   rnd = new TRandom3(666);
   
   _fout = new TFile("hist/output.root","RECREATE");

   _tr = new TTree("tr","tree");

   _tr->Branch("elTruth_n",&m_elTruth_n,"elTruth_n/I");
   _tr->Branch("elTruth_pt",&m_elTruth_pt,"elTruth_pt[elTruth_n]/F");
   _tr->Branch("elTruth_E",&m_elTruth_E,"elTruth_E[elTruth_n]/F");
   _tr->Branch("elTruth_eta",&m_elTruth_eta,"elTruth_eta[elTruth_n]/F");
   _tr->Branch("elRec_pt",&m_elRec_pt,"elRec_pt[elTruth_n]/F");
   _tr->Branch("elRec_E",&m_elRec_E,"elRec_E[elTruth_n]/F");
   _tr->Branch("elRec_eta",&m_elRec_eta,"elRec_eta[elTruth_n]/F");
   _tr->Branch("elTruth_label",&m_elTruth_label,"elTruth_label[elTruth_n]/I");

   _tr->Branch("muTruth_n",&m_muTruth_n,"muTruth_n/I");
   _tr->Branch("muTruth_pt",&m_muTruth_pt,"muTruth_pt[muTruth_n]/F");
   _tr->Branch("muTruth_E",&m_muTruth_E,"muTruth_E[muTruth_n]/F");
   _tr->Branch("muTruth_eta",&m_muTruth_eta,"muTruth_eta[muTruth_n]/F");
   _tr->Branch("muRec_pt",&m_muRec_pt,"muRec_pt[muTruth_n]/F");
   _tr->Branch("muRec_E",&m_muRec_E,"muRec_E[muTruth_n]/F");
   _tr->Branch("muRec_eta",&m_muRec_eta,"muRec_eta[muTruth_n]/F");
   _tr->Branch("muTruth_label",&m_muTruth_label,"muTruth_label[muTruth_n]/I");

   _tr->Branch("qTruth_n",&m_qTruth_n,"qTruth_n/I");
   _tr->Branch("qTruth_pt",&m_qTruth_pt,"qTruth_pt[qTruth_n]/F");
   _tr->Branch("qTruth_E",&m_qTruth_E,"qTruth_E[qTruth_n]/F");
   _tr->Branch("qTruth_eta",&m_qTruth_eta,"qTruth_eta[qTruth_n]/F");
   _tr->Branch("qRec_pt",&m_qRec_pt,"qRec_pt[qTruth_n]/F");
   _tr->Branch("qRec_E",&m_qRec_E,"qRec_E[qTruth_n]/F");
   _tr->Branch("qRec_eta",&m_qRec_eta,"qRec_eta[qTruth_n]/F");
   _tr->Branch("qTruth_label",&m_qTruth_label,"qTruth_label[qTruth_n]/I");

   _tr->Branch("bTruth_n",&m_bTruth_n,"bTruth_n/I");
   _tr->Branch("bTruth_pt",&m_bTruth_pt,"bTruth_pt[bTruth_n]/F");
   _tr->Branch("bTruth_E",&m_bTruth_E,"bTruth_E[bTruth_n]/F");
   _tr->Branch("bTruth_eta",&m_bTruth_eta,"bTruth_eta[bTruth_n]/F");
   _tr->Branch("bRec_pt",&m_bRec_pt,"bRec_pt[bTruth_n]/F");
   _tr->Branch("bRec_E",&m_bRec_E,"bRec_E[bTruth_n]/F");
   _tr->Branch("bRec_eta",&m_bRec_eta,"bRec_eta[bTruth_n]/F");
   _tr->Branch("bTruth_label",&m_bTruth_label,"bTruth_label[bTruth_n]/I");
   
   std::cout << "Initialisation done" << std::endl;
}

void TransferFunc::run()
{
   int nPart = _v_Truth->at(0).mc_truth_n();
   std::vector<int> labTruth = _v_Truth->at(0).mc_truth_label();
   
   m_elTruth_n = 0;
   m_muTruth_n = 0;
   m_qTruth_n = 0;
   m_bTruth_n = 0;
   
   for(int i=0;i<nPart;i++)
     {
	float l_elRec_pt = -666;
	float l_elRec_E = -666;
	float l_elRec_eta = -666;

	float l_muRec_pt = -666;
	float l_muRec_E = -666;
	float l_muRec_eta = -666;

	float l_qRec_pt = -666;
	float l_qRec_E = -666;
	float l_qRec_eta = -666;

	float l_bRec_pt = -666;
	float l_bRec_E = -666;
	float l_bRec_eta = -666;
	
	// hWl,hZl,tWl
	if( labTruth[i] == 120 || labTruth[i] == 130 ||
	    labTruth[i] == 140 || labTruth[i] == 141 ||
	    labTruth[i] == 150 || labTruth[i] == 151 ||
	    labTruth[i] == 220 || labTruth[i] == 330 )
	  {
	     TLorentzVector *vlPart = &(_v_Truth->at(0).mc_truth_p4()[i]);
	     int idPart = _v_Truth->at(0).mc_truth_id()[i];
	     float etaPart = vlPart->PseudoRapidity();
	     float phiPart = vlPart->Phi();
	     
	     if( abs(idPart) == 11 )
	       {
		  m_elTruth_pt[m_elTruth_n] = vlPart->Pt();
		  m_elTruth_E[m_elTruth_n] = vlPart->E();
		  m_elTruth_eta[m_elTruth_n] = vlPart->PseudoRapidity();
		  m_elTruth_label[m_elTruth_n] = idPart;
		  
		  float drMin = 10E+10;
		  int idMin = -666;
		  int nElec = _v_Electron->size();
		  for(int ie=0;ie<nElec;ie++)
		    {
		       if( !_v_Electron->at(ie).isLoose() ) continue;
		       
		       float ptEl = _v_Electron->at(ie).pt();
		       float etaEl = _v_Electron->at(ie).eta();
		       float phiEl = _v_Electron->at(ie).phi();
		       
		       float dr = GetDeltaR(etaPart,phiPart,etaEl,phiEl);
		       float ptDiff = (m_elTruth_pt[m_elTruth_n] > 0.) ? fabs(ptEl-m_elTruth_pt[m_elTruth_n])/m_elTruth_pt[m_elTruth_n] : 10E+10;
//		       if( dr < drMin && ptDiff < 0.3 )
		       if( dr < drMin )
			 {
			    drMin = dr;
			    idMin = ie;
			 }		       
		    }
		  if( drMin < 0.5 )
		    {
		       l_elRec_pt = _v_Electron->at(idMin).pt();
		       l_elRec_E = _v_Electron->at(idMin).E();
		       l_elRec_eta = _v_Electron->at(idMin).eta();
		    }
		  
		  m_elRec_pt[m_elTruth_n] = l_elRec_pt;
		  m_elRec_E[m_elTruth_n] = l_elRec_E;
		  m_elRec_eta[m_elTruth_n] = l_elRec_eta;
		  
		  m_elTruth_n++;
	       } // electrons

	     if( abs(idPart) == 13 )
	       {
		  m_muTruth_pt[m_muTruth_n] = vlPart->Pt();
		  m_muTruth_E[m_muTruth_n] = vlPart->E();
		  m_muTruth_eta[m_muTruth_n] = vlPart->PseudoRapidity();
		  m_muTruth_label[m_muTruth_n] = idPart;
		  
		  float drMin = 10E+10;
		  int idMin = -666;
		  int nMuon = _v_Muon->size();
		  for(int im=0;im<nMuon;im++)
		    {
		       if( !_v_Muon->at(im).isLoose() ) continue;
		       
		       float ptMu = _v_Muon->at(im).pt();
		       float etaMu = _v_Muon->at(im).eta();
		       float phiMu = _v_Muon->at(im).phi();
		       
		       float dr = GetDeltaR(etaPart,phiPart,etaMu,phiMu);
		       float ptDiff = (m_muTruth_pt[m_muTruth_n] > 0.) ? fabs(ptMu-m_muTruth_pt[m_muTruth_n])/m_muTruth_pt[m_muTruth_n] : 10E+10;
//		       if( dr < drMin && ptDiff < 0.3 )
		       if( dr < drMin )
			 {
			    drMin = dr;
			    idMin = im;
			 }		       
		    }
		  if( drMin < 0.5 )
		    {
		       l_muRec_pt = _v_Muon->at(idMin).pt();
		       l_muRec_E = _v_Muon->at(idMin).E();
		       l_muRec_eta = _v_Muon->at(idMin).eta();
		    }
		  
		  m_muRec_pt[m_muTruth_n] = l_muRec_pt;
		  m_muRec_E[m_muTruth_n] = l_muRec_E;
		  m_muRec_eta[m_muTruth_n] = l_muRec_eta;
		  
		  m_muTruth_n++;
	       } // muons
	  } // leptons
	
	// hWq,hZq (non-b)
	if( labTruth[i] == 122 || labTruth[i] == 123 ||
	    labTruth[i] == 132 || labTruth[i] == 133 ||
	    labTruth[i] == 142 || labTruth[i] == 143 ||
	    labTruth[i] == 152 || labTruth[i] == 153 ||
	    labTruth[i] == 222 || labTruth[i] == 223 ||
	    labTruth[i] == 332 || labTruth[i] == 333
	  )
	  {
	     TLorentzVector *vlPart = &(_v_Truth->at(0).mc_truth_p4()[i]);
	     int idPart = _v_Truth->at(0).mc_truth_id()[i];
	     float etaPart = vlPart->PseudoRapidity();
	     float phiPart = vlPart->Phi();

	     if( abs(idPart) < 5 )
	       {
		  m_qTruth_pt[m_qTruth_n] = vlPart->Pt();
		  m_qTruth_E[m_qTruth_n] = vlPart->E();
		  m_qTruth_eta[m_qTruth_n] = vlPart->PseudoRapidity();
		  m_qTruth_label[m_qTruth_n] = idPart;
		  
		  float drMin = 10E+10;
		  int idMin = -666;
		  int nJet = _v_Jet->size();
		  for(int ij=0;ij<nJet;ij++)
		    {
		       if( !_v_Jet->at(ij).isTight() ) continue;
		       
		       float ptJet = _v_Jet->at(ij).pt();
		       float etaJet = _v_Jet->at(ij).eta();
		       float phiJet = _v_Jet->at(ij).phi();

/*		       float drMin2 = 10E+10;
		       for(int ijj=0;ijj<nJet;ijj++)
			 {
			    if( ij == ijj ) continue; 
			    float etaJett = _v_Jet->at(ijj).eta();
			    float phiJett = _v_Jet->at(ijj).phi();
			    float dr = GetDeltaR(etaJett,phiJett,etaJet,phiJet);
			    if( dr < drMin2 ) drMin2 = dr;
			 }
		       if( drMin2 < 0.5 ) continue;*/
		       
		       float dr = GetDeltaR(etaPart,phiPart,etaJet,phiJet);
		       float ptDiff = (m_qTruth_pt[m_qTruth_n] > 0.) ? fabs(ptJet-m_qTruth_pt[m_qTruth_n])/m_qTruth_pt[m_qTruth_n] : 10E+10;
//		       if( dr < drMin && ptDiff < 0.3 )
		       if( dr < drMin )
			 {
			    drMin = dr;
			    idMin = ij;
			 }		       
		    }
		  if( drMin < 0.5 )
		    {
		       l_qRec_pt = _v_Jet->at(idMin).pt();
		       l_qRec_E = _v_Jet->at(idMin).E();
		       l_qRec_eta = _v_Jet->at(idMin).eta();
		    }
		  
		  m_qRec_pt[m_qTruth_n] = l_qRec_pt;
		  m_qRec_E[m_qTruth_n] = l_qRec_E;
		  m_qRec_eta[m_qTruth_n] = l_qRec_eta;
		  
		  m_qTruth_n++;
	       }	     
	  }

	// tb
	if( labTruth[i] == 21 || labTruth[i] == 31
	  )
	  {
	     TLorentzVector *vlPart = &(_v_Truth->at(0).mc_truth_p4()[i]);
	     int idPart = _v_Truth->at(0).mc_truth_id()[i];
	     float etaPart = vlPart->PseudoRapidity();
	     float phiPart = vlPart->Phi();

	     if( abs(idPart) == 5 )
	       {
		  m_bTruth_pt[m_bTruth_n] = vlPart->Pt();
		  m_bTruth_E[m_bTruth_n] = vlPart->E();
		  m_bTruth_eta[m_bTruth_n] = vlPart->PseudoRapidity();
		  m_bTruth_label[m_bTruth_n] = idPart;
		  
		  float drMin = 10E+10;
		  int idMin = -666;
		  int nJet = _v_Jet->size();
		  for(int ij=0;ij<nJet;ij++)
		    {
		       if( !_v_Jet->at(ij).isTight() ) continue;		       
		       if( _v_Jet->at(ij).CSV() < 0.244 ) continue;
		       
		       float ptJet = _v_Jet->at(ij).pt();
		       float etaJet = _v_Jet->at(ij).eta();
		       float phiJet = _v_Jet->at(ij).phi();
		       
		       float dr = GetDeltaR(etaPart,phiPart,etaJet,phiJet);
		       float ptDiff = (m_bTruth_pt[m_bTruth_n] > 0.) ? fabs(ptJet-m_bTruth_pt[m_bTruth_n])/m_bTruth_pt[m_bTruth_n] : 10E+10;
//		       if( dr < drMin && ptDiff < 0.3 )
		       if( dr < drMin )
			 {
			    drMin = dr;
			    idMin = ij;
			 }		       
		    }
		  if( drMin < 0.5 )
		    {
		       l_bRec_pt = _v_Jet->at(idMin).pt();
		       l_bRec_E = _v_Jet->at(idMin).E();
		       l_bRec_eta = _v_Jet->at(idMin).eta();
		    }
		  
		  m_bRec_pt[m_bTruth_n] = l_bRec_pt;
		  m_bRec_E[m_bTruth_n] = l_bRec_E;
		  m_bRec_eta[m_bTruth_n] = l_bRec_eta;
		  
		  m_bTruth_n++;
	       }
	  }	
     }
   
   _tr->Fill();
}

void TransferFunc::close()
{
   _fout->Write();
   _fout->Close();
}

float TransferFunc::GetDeltaR(float eta1,float phi1,float eta2,float phi2)
{
   float DeltaPhi = TMath::Abs(phi2 - phi1);
      if (DeltaPhi > 3.141593 ) DeltaPhi -= 2.*3.141593;
   return TMath::Sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}
