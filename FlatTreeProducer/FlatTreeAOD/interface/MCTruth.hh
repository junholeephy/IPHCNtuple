#ifndef MCTRUTH_H
#define MCTRUTH_H

#include <string>
#include <iostream>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FlatTreeProducer/FlatTreeAOD/interface/FlatTree.hh"

#define DEFVAL -666

class MCTruth
{
 public:
   
   MCTruth() {};
   
   void Init(FlatTree &tree);
   
   reco::GenParticle* getUnique(const reco::GenParticle* p,
				bool verbose);
   
   void fillTTHSignalGenParticles(const edm::Event& iEvent,
				  const edm::EventSetup& iSetup,
				  FlatTree& tree,
				  const edm::Handle<std::vector<reco::GenParticle> >& GenParticles);

   void fillTZQSignalGenParticles(const edm::Event& iEvent,
				  const edm::EventSetup& iSetup,
				  FlatTree& tree,
				  const edm::Handle<std::vector<reco::GenParticle> >& GenParticles);
   
   void p4toTLV(reco::Particle::LorentzVector vp4,
		TLorentzVector& tlv);
};

#endif
