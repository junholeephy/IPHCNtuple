#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "FlatTreeProducer/FlatTreeAOD/interface/tinyxml2.h"

#include "FlatTreeProducer/FlatTreeAOD/interface/FlatTree.hh"
#include "FlatTreeProducer/FlatTreeAOD/interface/MCTruth.hh"

//#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimator.h"

using namespace tinyxml2;

class FlatTreeAOD : public edm::EDAnalyzer 
{
 public:
   explicit FlatTreeAOD(const edm::ParameterSet&);
   ~FlatTreeAOD();
   
   static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
   
 private:
   virtual void beginJob() override;
   virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
   virtual void endJob() override;
   
   FlatTree* ftree;
   const edm::Service<TFileService> fs;
   
   TH1D* hcount;

//   EGammaMvaEleEstimator* elecMVA;
//   std::vector<std::string> elecMVACatWeights;
   
   XMLDocument xmlconf;
   
   edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
   edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
   
   edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
   edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
   edm::EDGetTokenT<pat::MuonCollection> muonToken_;
   edm::EDGetTokenT<pat::JetCollection> jetToken_;
};

FlatTreeAOD::FlatTreeAOD(const edm::ParameterSet& iConfig)
{
//   elecMVA = new EGammaMvaEleEstimator();
   
//   elecMVACatWeights.push_back("../../../EgammaAnalysis/ElectronTools/data/Electrons_BDTG_NonTrigV0_Cat1.weights.xml");
//   elecMVACatWeights.push_back("../../../EgammaAnalysis/ElectronTools/data/Electrons_BDTG_NonTrigV0_Cat2.weights.xml");
//   elecMVACatWeights.push_back("../../../EgammaAnalysis/ElectronTools/data/Electrons_BDTG_NonTrigV0_Cat3.weights.xml");
//   elecMVACatWeights.push_back("../../../EgammaAnalysis/ElectronTools/data/Electrons_BDTG_NonTrigV0_Cat4.weights.xml");
//   elecMVACatWeights.push_back("../../../EgammaAnalysis/ElectronTools/data/Electrons_BDTG_NonTrigV0_Cat5.weights.xml");
//   elecMVACatWeights.push_back("../../../EgammaAnalysis/ElectronTools/data/Electrons_BDTG_NonTrigV0_Cat6.weights.xml");
   
//   elecMVA->initialize("BDT",
//		       EGammaMvaEleEstimator::kNonTrig,
//		       true,
//		       elecMVACatWeights);
   
   ftree = new FlatTree(fs->make<TTree>("tree","tree"));
   
   xmlconf.LoadFile("conf.xml");
   XMLElement* tElement = xmlconf.FirstChildElement("var");
   
   for( XMLElement* child=tElement;child!=0;child=child->NextSiblingElement() )
     {	
	std::string vname = child->ToElement()->Attribute("name");
	std::string vsave = child->ToElement()->Attribute("save");
	bool bsave = atoi(vsave.c_str());

	ftree->conf.insert(std::make_pair(vname,bsave));
     }

   ftree->CreateBranches();
   
   hcount = fs->make<TH1D>("hcount","hcount",2,0.,2.);
   
   triggerBits_ = consumes<edm::TriggerResults>(edm::InputTag(std::string("TriggerResults"),std::string(""),std::string("HLT")));
   triggerPrescales_ = consumes<pat::PackedTriggerPrescales>(edm::InputTag(std::string("patTrigger")));
   
   vertexToken_ = consumes<reco::VertexCollection>(edm::InputTag(std::string("offlineSlimmedPrimaryVertices")));
   electronToken_ = consumes<pat::ElectronCollection>(edm::InputTag(std::string("slimmedElectrons")));
   muonToken_ = consumes<pat::MuonCollection>(edm::InputTag(std::string("slimmedMuons")));
   jetToken_ = consumes<pat::JetCollection>(edm::InputTag(std::string("slimmedJets")));
}

FlatTreeAOD::~FlatTreeAOD()
{
}

void FlatTreeAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   
   hcount->SetBinContent(1,hcount->GetBinContent(1)+1);
   
   ftree->Init();
   
   ftree->ev_run = iEvent.id().run();
   ftree->ev_id = iEvent.id().event();
   ftree->ev_lumi = iEvent.id().luminosityBlock();

   // general info
   edm::Handle<GenEventInfoProduct> genEventInfo;
   iEvent.getByLabel("generator",genEventInfo);
   
   ftree->mc_id = genEventInfo->signalProcessID();
   ftree->mc_f1 = genEventInfo->pdf()->id.first;
   ftree->mc_f2 = genEventInfo->pdf()->id.second;
   ftree->mc_x1 = genEventInfo->pdf()->x.first;
   ftree->mc_x2 = genEventInfo->pdf()->x.second;
   ftree->mc_scale = genEventInfo->pdf()->scalePDF;
   if( genEventInfo->binningValues().size() > 0 ) ftree->mc_ptHat = genEventInfo->binningValues()[0];

   // pileup
   edm::Handle<std::vector< PileupSummaryInfo> > pileupInfo;
   iEvent.getByLabel("addPileupInfo",pileupInfo);

   ftree->mc_pu_Npvi = pileupInfo->size();
   
   for(std::vector<PileupSummaryInfo>::const_iterator pvi=pileupInfo->begin();
       pvi!=pileupInfo->end();pvi++)
     {
	signed int n_bc = pvi->getBunchCrossing();
	ftree->mc_pu_BunchCrossing.push_back(n_bc);
	
	if( n_bc == 0 )
	  {
	     ftree->mc_pu_intime_NumInt = pvi->getPU_NumInteractions();
	     ftree->mc_pu_trueNumInt = pvi->getTrueNumInteractions();
	  }	
	else if( n_bc == -1 ) ftree->mc_pu_before_npu = pvi->getPU_NumInteractions();
	else if( n_bc == +1 ) ftree->mc_pu_after_npu  = pvi->getPU_NumInteractions();

	std::vector<float> mc_pu_zpositions;
	std::vector<float> mc_pu_sumpT_lowpT;
	std::vector<float> mc_pu_sumpT_highpT;
	std::vector<int> mc_pu_ntrks_lowpT;
	std::vector<int> mc_pu_ntrks_highpT;
	
	ftree->mc_pu_Nzpositions.push_back(pvi->getPU_zpositions().size());
	
	for( unsigned int ipu=0;ipu<pvi->getPU_zpositions().size();ipu++ )
	  {
	     mc_pu_zpositions.push_back((pvi->getPU_zpositions())[ipu]);
	     mc_pu_sumpT_lowpT.push_back((pvi->getPU_sumpT_lowpT())[ipu]);
	     mc_pu_sumpT_highpT.push_back((pvi->getPU_sumpT_highpT())[ipu]);
	     mc_pu_ntrks_lowpT.push_back((pvi->getPU_ntrks_lowpT())[ipu]);
	     mc_pu_ntrks_highpT.push_back((pvi->getPU_ntrks_highpT())[ipu]);
	  }
	
	ftree->mc_pu_zpositions.push_back(mc_pu_zpositions);
	ftree->mc_pu_sumpT_lowpT.push_back(mc_pu_sumpT_lowpT);
	ftree->mc_pu_sumpT_highpT.push_back(mc_pu_sumpT_highpT);
	ftree->mc_pu_ntrks_lowpT.push_back(mc_pu_ntrks_lowpT);
	ftree->mc_pu_ntrks_highpT.push_back(mc_pu_ntrks_highpT);
     }

   // mc truth
   edm::Handle<reco::GenParticleCollection> genParticlesHandle;                                                          
   iEvent.getByLabel("genParticles",genParticlesHandle);

   bool do_mc_truth_tth = ftree->doWrite("mc_truth_tth");
   bool do_mc_truth_tzq = ftree->doWrite("mc_truth_tzq");
   
   if( do_mc_truth_tth || do_mc_truth_tzq )
     {	
	MCTruth *mc_truth = new MCTruth();
	mc_truth->Init(*ftree);
	if( do_mc_truth_tth ) mc_truth->fillTTHSignalGenParticles(iEvent,iSetup,*ftree,genParticlesHandle);
	if( do_mc_truth_tzq ) mc_truth->fillTZQSignalGenParticles(iEvent,iSetup,*ftree,genParticlesHandle);
	delete mc_truth;
     }
   
   // trigger
   edm::Handle<edm::TriggerResults> triggerBits;
   iEvent.getByToken(triggerBits_,triggerBits);
   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

   edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
   iEvent.getByToken(triggerPrescales_, triggerPrescales);
   
   std::vector<std::string> triggerIdentifiers_;
   triggerIdentifiers_.push_back("HLT_Ele27_WP80_v*");
   
   for (unsigned int j = 0; j < triggerIdentifiers_.size(); ++j) 
     {
	std::string idName = triggerIdentifiers_[j];
	std::string idNameUnstarred = idName;
	bool isStarred = (idName.find("*")!=std::string::npos);
	if( isStarred ) idNameUnstarred.erase( idName.find("*"), 1 );
	
	for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i)
	  {
	     if( (isStarred && names.triggerName(i).find(idNameUnstarred)!=std::string::npos ) ||
		 (!isStarred && names.triggerName(i)==idName) )
		 {
		 }		 
//	     std::cout << "Trigger " << names.triggerName(i) <<
//	       ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
//	       ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)") << std::endl;
	  }
     }      

   // Primary vertex
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vertexToken_,vertices);
   
   reco::Vertex *primVtx = NULL;   
   if( ! vertices->empty() )
     {	
	const reco::Vertex &PV = vertices->front();
	primVtx = (reco::Vertex*)&PV;
     }   
   
   // Electrons
   edm::Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken_,electrons);
   
   int nElec = electrons->size();
   ftree->el_n = nElec;
   
   for(int ie=0;ie<nElec;ie++)
     {
	const pat::Electron& elec = electrons->at(ie);
	
	ftree->el_pt.push_back(elec.pt());
	ftree->el_eta.push_back(elec.eta());
	ftree->el_phi.push_back(elec.phi());
	ftree->el_m.push_back(elec.mass());
	ftree->el_E.push_back(elec.energy());
	ftree->el_id.push_back(elec.pdgId());
	ftree->el_charge.push_back(elec.charge());
	
	ftree->el_scleta.push_back(elec.superCluster()->eta());
	ftree->el_isGsfCtfScPixChargeConsistent.push_back(elec.isGsfCtfScPixChargeConsistent());
	ftree->el_sigmaIetaIeta.push_back(elec.sigmaIetaIeta());
	ftree->el_hadronicOverEm.push_back(elec.hadronicOverEm());
	ftree->el_dr03TkSumPt.push_back(elec.dr03TkSumPt());
	ftree->el_dr03EcalRecHitSumEt.push_back(elec.dr03EcalRecHitSumEt());
	ftree->el_dr03HcalTowerSumEt.push_back(elec.dr03HcalTowerSumEt());
	
	if( elec.gsfTrack().isNonnull() )
	  {
	     ftree->el_numberOfLostHits.push_back(elec.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits());
	  }	
	
/*	bool validKF = false;
	reco::TrackRef myTrackRef = elec.closestCtfTrackRef();
	validKF = (myTrackRef.isAvailable());
	validKF = (myTrackRef.isNonnull());
	
	ftree->el_fbrem.push_back(elec.fbrem());
	if( validKF ) ftree->el_kf_normalizedChi2.push_back(myTrackRef->normalizedChi2());
	if( validKF ) ftree->el_trackerLayersWithMeasurement.push_back(myTrackRef->hitPattern().trackerLayersWithMeasurement());
	if( elec.gsfTrack().isNonnull() ) ftree->el_gsf_normalizedChi2.push_back(elec.gsfTrack()->normalizedChi2());
	ftree->el_deltaEtaSuperClusterTrackAtVtx.push_back(elec.deltaEtaSuperClusterTrackAtVtx());
	ftree->el_deltaPhiSuperClusterTrackAtVtx.push_back(elec.deltaPhiSuperClusterTrackAtVtx());
	ftree->el_deltaEtaSeedClusterTrackAtCalo.push_back(elec.deltaEtaSeedClusterTrackAtCalo());*/
	
	ftree->el_dB3D.push_back(elec.dB());
	ftree->el_edB3D.push_back(elec.edB());

/*	double el_mvaNonTrigV0 = elecMVA->mvaValue( Var_fbrem,
						    Var_kfchi2,
						    Var_kfhits,
						    Var_gsfchi2,
						    Var_deta,
						    Var_dphi,
						    Var_detacalo,    
						    Var_see,
						    Var_spp,
						    Var_etawidth,
						    Var_phiwidth,
						    Var_e1x5e5x5,
						    Var_R9,    
						    Var_HoE,
						    Var_EoP,
						    Var_IoEmIoP,
						    Var_eleEoPout,
						    Var_rho,
						    Var_PreShowerOverRaw,    
						    Var_eta,
						    Var_pt,
						    printDebug);*/
	
	ftree->el_neutralHadronIso.push_back(elec.neutralHadronIso());
	ftree->el_chargedHadronIso.push_back(elec.chargedHadronIso());
	ftree->el_puChargedHadronIso.push_back(elec.puChargedHadronIso());
	ftree->el_ecalIso.push_back(elec.ecalIso());
	ftree->el_hcalIso.push_back(elec.hcalIso());
	ftree->el_particleIso.push_back(elec.particleIso());
	ftree->el_photonIso.push_back(elec.photonIso());
	ftree->el_trackIso.push_back(elec.trackIso());
	
	ftree->el_pfIso_sumChargedHadronPt.push_back(elec.pfIsolationVariables().sumChargedHadronPt);
	ftree->el_pfIso_sumNeutralHadronEt.push_back(elec.pfIsolationVariables().sumNeutralHadronEt);
	ftree->el_pfIso_sumPhotonEt.push_back(elec.pfIsolationVariables().sumPhotonEt);
	ftree->el_pfIso_sumPUPt.push_back(elec.pfIsolationVariables().sumPUPt);
	
	ftree->el_isLoose.push_back(elec.electronID("eidLoose"));
	ftree->el_isTight.push_back(elec.electronID("eidTight"));
	ftree->el_isRobustLoose.push_back(elec.electronID("eidRobustLoose"));
	ftree->el_isRobustTight.push_back(elec.electronID("eidRobustTight"));
	ftree->el_isRobustHighEnergy.push_back(elec.electronID("eidRobustHighEnergy"));
	
	ftree->el_vx.push_back(elec.vx());
	ftree->el_vy.push_back(elec.vy());
	ftree->el_vz.push_back(elec.vz());
	
	ftree->el_isGsf.push_back(elec.gsfTrack().isNonnull());
	
	ftree->el_passConversionVeto.push_back(elec.passConversionVeto());
	
	if( elec.gsfTrack().isNonnull() )
	  {
	     ftree->el_numberOfHits.push_back(elec.gsfTrack()->trackerExpectedHitsInner().numberOfHits());
	     
	     if( primVtx )
	       {		  
		  ftree->el_dxy.push_back(elec.gsfTrack()->dxy(primVtx->position()));
		  ftree->el_dz.push_back(elec.gsfTrack()->dz(primVtx->position()));
		  ftree->el_dxyError.push_back(elec.gsfTrack()->dxyError());
		  ftree->el_dzError.push_back(elec.gsfTrack()->dzError());
	       }	     
	  }
     }   

   // Muons
   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_,muons);
   
   int nMuon = muons->size();
   ftree->mu_n = nMuon;
   
   for(int im=0;im<nMuon;im++)
     {
	const pat::Muon& muon = muons->at(im);
	
	ftree->mu_pt.push_back(muon.pt());
	ftree->mu_eta.push_back(muon.eta());
	ftree->mu_phi.push_back(muon.phi());
	ftree->mu_m.push_back(muon.mass());
	ftree->mu_E.push_back(muon.energy());
	ftree->mu_id.push_back(muon.pdgId());
	ftree->mu_charge.push_back(muon.charge());
	
	ftree->mu_dB3D.push_back(muon.dB());
	ftree->mu_edB3D.push_back(muon.edB());
	
	ftree->mu_neutralHadronIso.push_back(muon.neutralHadronIso());
	ftree->mu_chargedHadronIso.push_back(muon.chargedHadronIso());
	ftree->mu_puChargedHadronIso.push_back(muon.puChargedHadronIso());
	ftree->mu_ecalIso.push_back(muon.ecalIso());
	ftree->mu_hcalIso.push_back(muon.hcalIso());
	ftree->mu_photonIso.push_back(muon.photonIso());
	ftree->mu_trackIso.push_back(muon.trackIso());
	
	ftree->mu_pfIso03_sumChargedHadronPt.push_back(muon.pfIsolationR03().sumChargedHadronPt);
	ftree->mu_pfIso03_sumNeutralHadronEt.push_back(muon.pfIsolationR03().sumNeutralHadronEt);
	ftree->mu_pfIso03_sumPhotonEt.push_back(muon.pfIsolationR03().sumPhotonEt);
	ftree->mu_pfIso03_sumPUPt.push_back(muon.pfIsolationR03().sumPUPt);
	
	ftree->mu_isGlobalMuon.push_back(muon.isGlobalMuon());
	ftree->mu_isTrackerMuon.push_back(muon.isTrackerMuon());
	ftree->mu_isStandAloneMuon.push_back(muon.isStandAloneMuon());
	ftree->mu_isCaloMuon.push_back(muon.isCaloMuon());
	ftree->mu_isPFMuon.push_back(muon.isPFMuon());
	
	ftree->mu_vx.push_back(muon.vx());
	ftree->mu_vy.push_back(muon.vy());
	ftree->mu_vz.push_back(muon.vz());
	
	if( primVtx )
	  {		  
	     if( muon.globalTrack().isNonnull() )
	       {		       
		  ftree->mu_globalTrack_dxy.push_back(muon.globalTrack()->dxy(primVtx->position()));
		  ftree->mu_globalTrack_dz.push_back(muon.globalTrack()->dz(primVtx->position()));
		  ftree->mu_globalTrack_dxyError.push_back(muon.globalTrack()->dxyError());
		  ftree->mu_globalTrack_dzError.push_back(muon.globalTrack()->dzError());
	       }		  
	     if( muon.innerTrack().isNonnull() )
	       {		       
		  ftree->mu_innerTrack_dxy.push_back(muon.innerTrack()->dxy(primVtx->position()));
		  ftree->mu_innerTrack_dz.push_back(muon.innerTrack()->dz(primVtx->position()));
		  ftree->mu_innerTrack_dxyError.push_back(muon.innerTrack()->dxyError());
		  ftree->mu_innerTrack_dzError.push_back(muon.innerTrack()->dzError());
	       }		  
	  }		

	if( muon.innerTrack().isNonnull() )
	  {		       
	     ftree->mu_innerTrack_pt.push_back(muon.innerTrack()->pt());
	     ftree->mu_innerTrack_ptError.push_back(muon.innerTrack()->ptError());
	  }
	
	ftree->mu_numberOfMatches.push_back(muon.numberOfMatches());
	
	if( muon.globalTrack().isNonnull() )
	  {
	     ftree->mu_numberOfValidMuonHits.push_back(muon.globalTrack()->hitPattern().numberOfValidMuonHits());
	  }	
     }   

   // Jets
   edm::Handle<pat::JetCollection> jets;
   iEvent.getByToken(jetToken_,jets);
   
   int nJet = jets->size();
   ftree->jet_n = nJet;
   
   for(int ij=0;ij<nJet;ij++)
     {
	const pat::Jet& jet = jets->at(ij);
	
	ftree->jet_pt.push_back(jet.pt());
	ftree->jet_eta.push_back(jet.eta());
	ftree->jet_phi.push_back(jet.phi());
	ftree->jet_m.push_back(jet.mass());
	ftree->jet_E.push_back(jet.energy());
	
	ftree->jet_CSV.push_back(jet.bDiscriminator("combinedSecondaryVertexBJetTags"));
	ftree->jet_flavour.push_back(jet.partonFlavour());
	
	ftree->jet_neutralHadronEnergy.push_back(jet.neutralHadronEnergy());
	ftree->jet_neutralEmEnergy.push_back(jet.neutralEmEnergy());
	ftree->jet_chargedHadronEnergy.push_back(jet.chargedHadronEnergy());
	ftree->jet_chargedEmEnergy.push_back(jet.chargedEmEnergy());
	ftree->jet_electronEnergy.push_back(jet.electronEnergy());
	ftree->jet_muonEnergy.push_back(jet.muonEnergy());
	ftree->jet_photonEnergy.push_back(jet.photonEnergy());
	
	ftree->jet_pileupJetId.push_back(jet.userFloat("pileupJetId:fullDiscriminant"));
	
	const reco::GenJet* genJet = jet.genJet();
	if( genJet )
	  {
	     ftree->jet_gen_pt.push_back(genJet->pt());
	     ftree->jet_gen_eta.push_back(genJet->eta());
	     ftree->jet_gen_phi.push_back(genJet->phi());
	     ftree->jet_gen_m.push_back(genJet->mass());
	     ftree->jet_gen_E.push_back(genJet->energy());
	     
	     ftree->jet_gen_status.push_back(genJet->status());
	     ftree->jet_gen_id.push_back(genJet->pdgId());
	  }	
     }      
   
   ftree->tree->Fill();
}

void FlatTreeAOD::beginJob()
{
}

void FlatTreeAOD::endJob() 
{
}

void FlatTreeAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{   
   //The following says we do not know what parameters are allowed so do no validation
   // Please change this to state exactly what you do use, even if it is no parameters
   edm::ParameterSetDescription desc;
   desc.setUnknown();
   descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FlatTreeAOD);
