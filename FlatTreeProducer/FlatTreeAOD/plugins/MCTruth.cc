#include "FlatTreeProducer/FlatTreeAOD/interface/MCTruth.hh"

void MCTruth::Init(FlatTree &tree)
{
   tree.mc_truth_tth_channel = DEFVAL;

   // TLV
   
   tree.mc_truth_tth_h0_p4.Clear();

   tree.mc_truth_tth_h0W1_p4.Clear();
   tree.mc_truth_tth_h0W2_p4.Clear();
   tree.mc_truth_tth_h0Wl1_p4.Clear();
   tree.mc_truth_tth_h0Wnu1_p4.Clear();
   tree.mc_truth_tth_h0Wtau1_p4.Clear();
   tree.mc_truth_tth_h0Wnutau1_p4.Clear();
   tree.mc_truth_tth_h0Wtaul1_p4.Clear();
   tree.mc_truth_tth_h0Wtaunu1_p4.Clear();
   tree.mc_truth_tth_h0Wtaunutau1_p4.Clear();
   tree.mc_truth_tth_h0Wl2_p4.Clear();
   tree.mc_truth_tth_h0Wnu2_p4.Clear();
   tree.mc_truth_tth_h0Wtau2_p4.Clear();
   tree.mc_truth_tth_h0Wnutau2_p4.Clear();
   tree.mc_truth_tth_h0Wtaul2_p4.Clear();
   tree.mc_truth_tth_h0Wtaunu2_p4.Clear();
   tree.mc_truth_tth_h0Wtaunutau2_p4.Clear();
   tree.mc_truth_tth_h0Wq11_p4.Clear();
   tree.mc_truth_tth_h0Wq21_p4.Clear();
   tree.mc_truth_tth_h0Wq12_p4.Clear();
   tree.mc_truth_tth_h0Wq22_p4.Clear();
   
   tree.mc_truth_tth_h0Z1_p4.Clear();
   tree.mc_truth_tth_h0Z2_p4.Clear();
   tree.mc_truth_tth_h0Zl11_p4.Clear();
   tree.mc_truth_tth_h0Zl21_p4.Clear();
   tree.mc_truth_tth_h0Ztau11_p4.Clear();
   tree.mc_truth_tth_h0Ztau21_p4.Clear();
   tree.mc_truth_tth_h0Ztaul11_p4.Clear();
   tree.mc_truth_tth_h0Ztaul21_p4.Clear();
   tree.mc_truth_tth_h0Ztaunu11_p4.Clear();
   tree.mc_truth_tth_h0Ztaunu21_p4.Clear();
   tree.mc_truth_tth_h0Ztaunutau11_p4.Clear();
   tree.mc_truth_tth_h0Ztaunutau21_p4.Clear();
   tree.mc_truth_tth_h0Zq11_p4.Clear();
   tree.mc_truth_tth_h0Zq21_p4.Clear();
   tree.mc_truth_tth_h0Zl12_p4.Clear();
   tree.mc_truth_tth_h0Zl22_p4.Clear();
   tree.mc_truth_tth_h0Ztau12_p4.Clear();
   tree.mc_truth_tth_h0Ztau22_p4.Clear();
   tree.mc_truth_tth_h0Ztaul12_p4.Clear();
   tree.mc_truth_tth_h0Ztaul22_p4.Clear();
   tree.mc_truth_tth_h0Ztaunu12_p4.Clear();
   tree.mc_truth_tth_h0Ztaunu22_p4.Clear();
   tree.mc_truth_tth_h0Ztaunutau12_p4.Clear();
   tree.mc_truth_tth_h0Ztaunutau22_p4.Clear();
   tree.mc_truth_tth_h0Zq12_p4.Clear();
   tree.mc_truth_tth_h0Zq22_p4.Clear();
   tree.mc_truth_tth_h0Znu11_p4.Clear();
   tree.mc_truth_tth_h0Znu21_p4.Clear();
   tree.mc_truth_tth_h0Znu12_p4.Clear();
   tree.mc_truth_tth_h0Znu22_p4.Clear();
   
   tree.mc_truth_tth_h0tau1_p4.Clear();
   tree.mc_truth_tth_h0tau2_p4.Clear();
   tree.mc_truth_tth_h0taul1_p4.Clear();
   tree.mc_truth_tth_h0taunutau1_p4.Clear();
   tree.mc_truth_tth_h0taunu1_p4.Clear();
   tree.mc_truth_tth_h0taul2_p4.Clear();
   tree.mc_truth_tth_h0taunutau2_p4.Clear();
   tree.mc_truth_tth_h0taunu2_p4.Clear();
   
   tree.mc_truth_tth_t1_p4.Clear();
   tree.mc_truth_tth_t2_p4.Clear();
   tree.mc_truth_tth_tb1_p4.Clear();
   tree.mc_truth_tth_tb2_p4.Clear();
   
   tree.mc_truth_tth_tW1_p4.Clear();
   tree.mc_truth_tth_tWnu1_p4.Clear();
   tree.mc_truth_tth_tWnutau1_p4.Clear();
   tree.mc_truth_tth_tWl1_p4.Clear();
   tree.mc_truth_tth_tWtau1_p4.Clear();
   tree.mc_truth_tth_tWtaunu1_p4.Clear();
   tree.mc_truth_tth_tWtaunutau1_p4.Clear();
   tree.mc_truth_tth_tWtaul1_p4.Clear();
   tree.mc_truth_tth_tWq11_p4.Clear();
   tree.mc_truth_tth_tWq21_p4.Clear();

   tree.mc_truth_tth_tW2_p4.Clear();
   tree.mc_truth_tth_tWnu2_p4.Clear();
   tree.mc_truth_tth_tWnutau2_p4.Clear();
   tree.mc_truth_tth_tWl2_p4.Clear();
   tree.mc_truth_tth_tWtau2_p4.Clear();
   tree.mc_truth_tth_tWtaunu2_p4.Clear();
   tree.mc_truth_tth_tWtaunutau2_p4.Clear();
   tree.mc_truth_tth_tWtaul2_p4.Clear();
   tree.mc_truth_tth_tWq12_p4.Clear();
   tree.mc_truth_tth_tWq22_p4.Clear();

   tree.mc_truth_tth_j1_p4.Clear();
   tree.mc_truth_tth_j2_p4.Clear();
   tree.mc_truth_tth_j3_p4.Clear();

   // pdgId
   
   tree.mc_truth_tth_h0_id = DEFVAL;

   tree.mc_truth_tth_h0W1_id = DEFVAL;
   tree.mc_truth_tth_h0W2_id = DEFVAL;
   tree.mc_truth_tth_h0Wl1_id = DEFVAL;
   tree.mc_truth_tth_h0Wnu1_id = DEFVAL;
   tree.mc_truth_tth_h0Wtau1_id = DEFVAL;
   tree.mc_truth_tth_h0Wnutau1_id = DEFVAL;
   tree.mc_truth_tth_h0Wtaul1_id = DEFVAL;
   tree.mc_truth_tth_h0Wtaunu1_id = DEFVAL;
   tree.mc_truth_tth_h0Wtaunutau1_id = DEFVAL;
   tree.mc_truth_tth_h0Wl2_id = DEFVAL;
   tree.mc_truth_tth_h0Wnu2_id = DEFVAL;
   tree.mc_truth_tth_h0Wtau2_id = DEFVAL;
   tree.mc_truth_tth_h0Wnutau2_id = DEFVAL;
   tree.mc_truth_tth_h0Wtaul2_id = DEFVAL;
   tree.mc_truth_tth_h0Wtaunu2_id = DEFVAL;
   tree.mc_truth_tth_h0Wtaunutau2_id = DEFVAL;
   tree.mc_truth_tth_h0Wq11_id = DEFVAL;
   tree.mc_truth_tth_h0Wq21_id = DEFVAL;
   tree.mc_truth_tth_h0Wq12_id = DEFVAL;
   tree.mc_truth_tth_h0Wq22_id = DEFVAL;
   
   tree.mc_truth_tth_h0Z1_id = DEFVAL;
   tree.mc_truth_tth_h0Z2_id = DEFVAL;
   tree.mc_truth_tth_h0Zl11_id = DEFVAL;
   tree.mc_truth_tth_h0Zl21_id = DEFVAL;
   tree.mc_truth_tth_h0Ztau11_id = DEFVAL;
   tree.mc_truth_tth_h0Ztau21_id = DEFVAL;
   tree.mc_truth_tth_h0Ztaul11_id = DEFVAL;
   tree.mc_truth_tth_h0Ztaul21_id = DEFVAL;
   tree.mc_truth_tth_h0Ztaunu11_id = DEFVAL;
   tree.mc_truth_tth_h0Ztaunu21_id = DEFVAL;
   tree.mc_truth_tth_h0Ztaunutau11_id = DEFVAL;
   tree.mc_truth_tth_h0Ztaunutau21_id = DEFVAL;
   tree.mc_truth_tth_h0Zq11_id = DEFVAL;
   tree.mc_truth_tth_h0Zq21_id = DEFVAL;
   tree.mc_truth_tth_h0Zl12_id = DEFVAL;
   tree.mc_truth_tth_h0Zl22_id = DEFVAL;
   tree.mc_truth_tth_h0Ztau12_id = DEFVAL;
   tree.mc_truth_tth_h0Ztau22_id = DEFVAL;
   tree.mc_truth_tth_h0Ztaul12_id = DEFVAL;
   tree.mc_truth_tth_h0Ztaul22_id = DEFVAL;
   tree.mc_truth_tth_h0Ztaunu12_id = DEFVAL;
   tree.mc_truth_tth_h0Ztaunu22_id = DEFVAL;
   tree.mc_truth_tth_h0Ztaunutau12_id = DEFVAL;
   tree.mc_truth_tth_h0Ztaunutau22_id = DEFVAL;
   tree.mc_truth_tth_h0Zq12_id = DEFVAL;
   tree.mc_truth_tth_h0Zq22_id = DEFVAL;
   tree.mc_truth_tth_h0Znu11_id = DEFVAL;
   tree.mc_truth_tth_h0Znu21_id = DEFVAL;
   tree.mc_truth_tth_h0Znu12_id = DEFVAL;
   tree.mc_truth_tth_h0Znu22_id = DEFVAL;
   
   tree.mc_truth_tth_h0tau1_id = DEFVAL;
   tree.mc_truth_tth_h0tau2_id = DEFVAL;
   tree.mc_truth_tth_h0taul1_id = DEFVAL;
   tree.mc_truth_tth_h0taunutau1_id = DEFVAL;
   tree.mc_truth_tth_h0taunu1_id = DEFVAL;
   tree.mc_truth_tth_h0taul2_id = DEFVAL;
   tree.mc_truth_tth_h0taunutau2_id = DEFVAL;
   tree.mc_truth_tth_h0taunu2_id = DEFVAL;
   
   tree.mc_truth_tth_t1_id = DEFVAL;
   tree.mc_truth_tth_t2_id = DEFVAL;
   tree.mc_truth_tth_tb1_id = DEFVAL;
   tree.mc_truth_tth_tb2_id = DEFVAL;
   
   tree.mc_truth_tth_tW1_id = DEFVAL;
   tree.mc_truth_tth_tWnu1_id = DEFVAL;
   tree.mc_truth_tth_tWnutau1_id = DEFVAL;
   tree.mc_truth_tth_tWl1_id = DEFVAL;
   tree.mc_truth_tth_tWtau1_id = DEFVAL;
   tree.mc_truth_tth_tWtaunu1_id = DEFVAL;
   tree.mc_truth_tth_tWtaunutau1_id = DEFVAL;
   tree.mc_truth_tth_tWtaul1_id = DEFVAL;
   tree.mc_truth_tth_tWq11_id = DEFVAL;
   tree.mc_truth_tth_tWq21_id = DEFVAL;

   tree.mc_truth_tth_tW2_id = DEFVAL;
   tree.mc_truth_tth_tWnu2_id = DEFVAL;
   tree.mc_truth_tth_tWnutau2_id = DEFVAL;
   tree.mc_truth_tth_tWl2_id = DEFVAL;
   tree.mc_truth_tth_tWtau2_id = DEFVAL;
   tree.mc_truth_tth_tWtaunu2_id = DEFVAL;
   tree.mc_truth_tth_tWtaunutau2_id = DEFVAL;
   tree.mc_truth_tth_tWtaul2_id = DEFVAL;
   tree.mc_truth_tth_tWq12_id = DEFVAL;
   tree.mc_truth_tth_tWq22_id = DEFVAL;

   tree.mc_truth_tth_j1_id = DEFVAL;
   tree.mc_truth_tth_j2_id = DEFVAL;
   tree.mc_truth_tth_j3_id = DEFVAL;

   // status

   tree.mc_truth_tth_h0_status = DEFVAL;

   tree.mc_truth_tth_h0W1_status = DEFVAL;
   tree.mc_truth_tth_h0W2_status = DEFVAL;
   tree.mc_truth_tth_h0Wl1_status = DEFVAL;
   tree.mc_truth_tth_h0Wnu1_status = DEFVAL;
   tree.mc_truth_tth_h0Wtau1_status = DEFVAL;
   tree.mc_truth_tth_h0Wnutau1_status = DEFVAL;
   tree.mc_truth_tth_h0Wtaul1_status = DEFVAL;
   tree.mc_truth_tth_h0Wtaunu1_status = DEFVAL;
   tree.mc_truth_tth_h0Wtaunutau1_status = DEFVAL;
   tree.mc_truth_tth_h0Wl2_status = DEFVAL;
   tree.mc_truth_tth_h0Wnu2_status = DEFVAL;
   tree.mc_truth_tth_h0Wtau2_status = DEFVAL;
   tree.mc_truth_tth_h0Wnutau2_status = DEFVAL;
   tree.mc_truth_tth_h0Wtaul2_status = DEFVAL;
   tree.mc_truth_tth_h0Wtaunu2_status = DEFVAL;
   tree.mc_truth_tth_h0Wtaunutau2_status = DEFVAL;
   tree.mc_truth_tth_h0Wq11_status = DEFVAL;
   tree.mc_truth_tth_h0Wq21_status = DEFVAL;
   tree.mc_truth_tth_h0Wq12_status = DEFVAL;
   tree.mc_truth_tth_h0Wq22_status = DEFVAL;
   
   tree.mc_truth_tth_h0Z1_status = DEFVAL;
   tree.mc_truth_tth_h0Z2_status = DEFVAL;
   tree.mc_truth_tth_h0Zl11_status = DEFVAL;
   tree.mc_truth_tth_h0Zl21_status = DEFVAL;
   tree.mc_truth_tth_h0Ztau11_status = DEFVAL;
   tree.mc_truth_tth_h0Ztau21_status = DEFVAL;
   tree.mc_truth_tth_h0Ztaul11_status = DEFVAL;
   tree.mc_truth_tth_h0Ztaul21_status = DEFVAL;
   tree.mc_truth_tth_h0Ztaunu11_status = DEFVAL;
   tree.mc_truth_tth_h0Ztaunu21_status = DEFVAL;
   tree.mc_truth_tth_h0Ztaunutau11_status = DEFVAL;
   tree.mc_truth_tth_h0Ztaunutau21_status = DEFVAL;
   tree.mc_truth_tth_h0Zq11_status = DEFVAL;
   tree.mc_truth_tth_h0Zq21_status = DEFVAL;
   tree.mc_truth_tth_h0Zl12_status = DEFVAL;
   tree.mc_truth_tth_h0Zl22_status = DEFVAL;
   tree.mc_truth_tth_h0Ztau12_status = DEFVAL;
   tree.mc_truth_tth_h0Ztau22_status = DEFVAL;
   tree.mc_truth_tth_h0Ztaul12_status = DEFVAL;
   tree.mc_truth_tth_h0Ztaul22_status = DEFVAL;
   tree.mc_truth_tth_h0Ztaunu12_status = DEFVAL;
   tree.mc_truth_tth_h0Ztaunu22_status = DEFVAL;
   tree.mc_truth_tth_h0Ztaunutau12_status = DEFVAL;
   tree.mc_truth_tth_h0Ztaunutau22_status = DEFVAL;
   tree.mc_truth_tth_h0Zq12_status = DEFVAL;
   tree.mc_truth_tth_h0Zq22_status = DEFVAL;
   tree.mc_truth_tth_h0Znu11_status = DEFVAL;
   tree.mc_truth_tth_h0Znu21_status = DEFVAL;
   tree.mc_truth_tth_h0Znu12_status = DEFVAL;
   tree.mc_truth_tth_h0Znu22_status = DEFVAL;
   
   tree.mc_truth_tth_h0tau1_status = DEFVAL;
   tree.mc_truth_tth_h0tau2_status = DEFVAL;
   tree.mc_truth_tth_h0taul1_status = DEFVAL;
   tree.mc_truth_tth_h0taunutau1_status = DEFVAL;
   tree.mc_truth_tth_h0taunu1_status = DEFVAL;
   tree.mc_truth_tth_h0taul2_status = DEFVAL;
   tree.mc_truth_tth_h0taunutau2_status = DEFVAL;
   tree.mc_truth_tth_h0taunu2_status = DEFVAL;
   
   tree.mc_truth_tth_t1_status = DEFVAL;
   tree.mc_truth_tth_t2_status = DEFVAL;
   tree.mc_truth_tth_tb1_status = DEFVAL;
   tree.mc_truth_tth_tb2_status = DEFVAL;
   
   tree.mc_truth_tth_tW1_status = DEFVAL;
   tree.mc_truth_tth_tWnu1_status = DEFVAL;
   tree.mc_truth_tth_tWnutau1_status = DEFVAL;
   tree.mc_truth_tth_tWl1_status = DEFVAL;
   tree.mc_truth_tth_tWtau1_status = DEFVAL;
   tree.mc_truth_tth_tWtaunu1_status = DEFVAL;
   tree.mc_truth_tth_tWtaunutau1_status = DEFVAL;
   tree.mc_truth_tth_tWtaul1_status = DEFVAL;
   tree.mc_truth_tth_tWq11_status = DEFVAL;
   tree.mc_truth_tth_tWq21_status = DEFVAL;

   tree.mc_truth_tth_tW2_status = DEFVAL;
   tree.mc_truth_tth_tWnu2_status = DEFVAL;
   tree.mc_truth_tth_tWnutau2_status = DEFVAL;
   tree.mc_truth_tth_tWl2_status = DEFVAL;
   tree.mc_truth_tth_tWtau2_status = DEFVAL;
   tree.mc_truth_tth_tWtaunu2_status = DEFVAL;
   tree.mc_truth_tth_tWtaunutau2_status = DEFVAL;
   tree.mc_truth_tth_tWtaul2_status = DEFVAL;
   tree.mc_truth_tth_tWq12_status = DEFVAL;
   tree.mc_truth_tth_tWq22_status = DEFVAL;

   tree.mc_truth_tth_j1_status = DEFVAL;
   tree.mc_truth_tth_j2_status = DEFVAL;
   tree.mc_truth_tth_j3_status = DEFVAL;

   // tZq
   tree.mc_truth_tzq_channel = DEFVAL;

   // TLV

   tree.mc_truth_tzq_Z_p4.Clear();
   tree.mc_truth_tzq_Zl1_p4.Clear();
   tree.mc_truth_tzq_Zl2_p4.Clear();
   tree.mc_truth_tzq_Ztau1_p4.Clear();
   tree.mc_truth_tzq_Ztau2_p4.Clear();
   tree.mc_truth_tzq_Ztaul1_p4.Clear();
   tree.mc_truth_tzq_Ztaul2_p4.Clear();
   tree.mc_truth_tzq_Ztaunu1_p4.Clear();
   tree.mc_truth_tzq_Ztaunu2_p4.Clear();
   tree.mc_truth_tzq_Ztaunutau1_p4.Clear();
   tree.mc_truth_tzq_Ztaunutau2_p4.Clear();
   
   tree.mc_truth_tzq_t_p4.Clear();
   tree.mc_truth_tzq_tb_p4.Clear();
   tree.mc_truth_tzq_tW_p4.Clear();
   tree.mc_truth_tzq_tWnu_p4.Clear();
   tree.mc_truth_tzq_tWnutau_p4.Clear();
   tree.mc_truth_tzq_tWl_p4.Clear();
   tree.mc_truth_tzq_tWtau_p4.Clear();
   tree.mc_truth_tzq_tWtaunu_p4.Clear();
   tree.mc_truth_tzq_tWtaunutau_p4.Clear();
   tree.mc_truth_tzq_tWtaul_p4.Clear();
   tree.mc_truth_tzq_tWq1_p4.Clear();
   tree.mc_truth_tzq_tWq2_p4.Clear();

   tree.mc_truth_tzq_j1_p4.Clear();
   tree.mc_truth_tzq_j2_p4.Clear();
   tree.mc_truth_tzq_j3_p4.Clear();

   // pdgId

   tree.mc_truth_tzq_Z_id = DEFVAL;
   tree.mc_truth_tzq_Zl1_id = DEFVAL;
   tree.mc_truth_tzq_Zl2_id = DEFVAL;
   tree.mc_truth_tzq_Ztau1_id = DEFVAL;
   tree.mc_truth_tzq_Ztau2_id = DEFVAL;
   tree.mc_truth_tzq_Ztaul1_id = DEFVAL;
   tree.mc_truth_tzq_Ztaul2_id = DEFVAL;
   tree.mc_truth_tzq_Ztaunu1_id = DEFVAL;
   tree.mc_truth_tzq_Ztaunu2_id = DEFVAL;
   tree.mc_truth_tzq_Ztaunutau1_id = DEFVAL;
   tree.mc_truth_tzq_Ztaunutau2_id = DEFVAL;
   
   tree.mc_truth_tzq_t_id = DEFVAL;
   tree.mc_truth_tzq_tb_id = DEFVAL;
   tree.mc_truth_tzq_tW_id = DEFVAL;
   tree.mc_truth_tzq_tWnu_id = DEFVAL;
   tree.mc_truth_tzq_tWnutau_id = DEFVAL;
   tree.mc_truth_tzq_tWl_id = DEFVAL;
   tree.mc_truth_tzq_tWtau_id = DEFVAL;
   tree.mc_truth_tzq_tWtaunu_id = DEFVAL;
   tree.mc_truth_tzq_tWtaunutau_id = DEFVAL;
   tree.mc_truth_tzq_tWtaul_id = DEFVAL;
   tree.mc_truth_tzq_tWq1_id = DEFVAL;
   tree.mc_truth_tzq_tWq2_id = DEFVAL;
   
   tree.mc_truth_tzq_j1_id = DEFVAL;
   tree.mc_truth_tzq_j2_id = DEFVAL;
   tree.mc_truth_tzq_j3_id = DEFVAL;

   // status
   
   tree.mc_truth_tzq_Z_status = DEFVAL;
   tree.mc_truth_tzq_Zl1_status = DEFVAL;
   tree.mc_truth_tzq_Zl2_status = DEFVAL;
   tree.mc_truth_tzq_Ztau1_status = DEFVAL;
   tree.mc_truth_tzq_Ztau2_status = DEFVAL;
   tree.mc_truth_tzq_Ztaul1_status = DEFVAL;
   tree.mc_truth_tzq_Ztaul2_status = DEFVAL;
   tree.mc_truth_tzq_Ztaunu1_status = DEFVAL;
   tree.mc_truth_tzq_Ztaunu2_status = DEFVAL;
   tree.mc_truth_tzq_Ztaunutau1_status = DEFVAL;
   tree.mc_truth_tzq_Ztaunutau2_status = DEFVAL;
   
   tree.mc_truth_tzq_t_status = DEFVAL;
   tree.mc_truth_tzq_tb_status = DEFVAL;
   tree.mc_truth_tzq_tW_status = DEFVAL;
   tree.mc_truth_tzq_tWnu_status = DEFVAL;
   tree.mc_truth_tzq_tWnutau_status = DEFVAL;
   tree.mc_truth_tzq_tWl_status = DEFVAL;
   tree.mc_truth_tzq_tWtau_status = DEFVAL;
   tree.mc_truth_tzq_tWtaunu_status = DEFVAL;
   tree.mc_truth_tzq_tWtaunutau_status = DEFVAL;
   tree.mc_truth_tzq_tWtaul_status = DEFVAL;
   tree.mc_truth_tzq_tWq1_status = DEFVAL;
   tree.mc_truth_tzq_tWq2_status = DEFVAL;
   
   tree.mc_truth_tzq_j1_status = DEFVAL;
   tree.mc_truth_tzq_j2_status = DEFVAL;
   tree.mc_truth_tzq_j3_status = DEFVAL;
}

void MCTruth::fillTTHSignalGenParticles(const edm::Event& iEvent,
					const edm::EventSetup& iSetup,
					FlatTree& tree,
					const edm::Handle<std::vector<reco::GenParticle> >& GenParticles)
{
   reco::GenParticle *h0 = 0;
   
   reco::GenParticle *h0W1 = 0;
   reco::GenParticle *h0W2 = 0;
   reco::GenParticle *h0Wl1 = 0;
   reco::GenParticle *h0Wnu1 = 0;
   reco::GenParticle *h0Wtau1 = 0;
   reco::GenParticle *h0Wnutau1 = 0;
   reco::GenParticle *h0Wtaul1 = 0;
   reco::GenParticle *h0Wtaunu1 = 0;
   reco::GenParticle *h0Wtaunutau1 = 0;
   reco::GenParticle *h0Wl2 = 0;
   reco::GenParticle *h0Wnu2 = 0;
   reco::GenParticle *h0Wtau2 = 0;
   reco::GenParticle *h0Wnutau2 = 0;
   reco::GenParticle *h0Wtaul2 = 0;
   reco::GenParticle *h0Wtaunu2 = 0;
   reco::GenParticle *h0Wtaunutau2 = 0;
   reco::GenParticle *h0Wq11 = 0;
   reco::GenParticle *h0Wq21 = 0;
   reco::GenParticle *h0Wq12 = 0;
   reco::GenParticle *h0Wq22 = 0;

   reco::GenParticle *h0Z1 = 0;
   reco::GenParticle *h0Z2 = 0;
   reco::GenParticle *h0Zl11 = 0;
   reco::GenParticle *h0Zl21 = 0;
   reco::GenParticle *h0Ztau11 = 0;
   reco::GenParticle *h0Ztau21 = 0;
   reco::GenParticle *h0Ztaul11 = 0;
   reco::GenParticle *h0Ztaul21 = 0;
   reco::GenParticle *h0Ztaunu11 = 0;
   reco::GenParticle *h0Ztaunu21 = 0;
   reco::GenParticle *h0Ztaunutau11 = 0;
   reco::GenParticle *h0Ztaunutau21 = 0;
   reco::GenParticle *h0Zq11 = 0;
   reco::GenParticle *h0Zq21 = 0;
   reco::GenParticle *h0Zl12 = 0;
   reco::GenParticle *h0Zl22 = 0;
   reco::GenParticle *h0Ztau12 = 0;
   reco::GenParticle *h0Ztau22 = 0;
   reco::GenParticle *h0Ztaul12 = 0;
   reco::GenParticle *h0Ztaul22 = 0;
   reco::GenParticle *h0Ztaunu12 = 0;
   reco::GenParticle *h0Ztaunu22 = 0;
   reco::GenParticle *h0Ztaunutau12 = 0;
   reco::GenParticle *h0Ztaunutau22 = 0;
   reco::GenParticle *h0Zq12 = 0;
   reco::GenParticle *h0Zq22 = 0;
   reco::GenParticle *h0Znu11 = 0;
   reco::GenParticle *h0Znu21 = 0;
   reco::GenParticle *h0Znu12 = 0;
   reco::GenParticle *h0Znu22 = 0;
   
   reco::GenParticle *h0tau1 = 0;
   reco::GenParticle *h0tau2 = 0;
   reco::GenParticle *h0taul1 = 0;
   reco::GenParticle *h0taunutau1 = 0;
   reco::GenParticle *h0taunu1 = 0;
   reco::GenParticle *h0taul2 = 0;
   reco::GenParticle *h0taunutau2 = 0;
   reco::GenParticle *h0taunu2 = 0;
   
   reco::GenParticle *t1 = 0;
   reco::GenParticle *t2 = 0;   

   reco::GenParticle *tb1 = 0;
   reco::GenParticle *tb2 = 0;
   
   reco::GenParticle *tW1 = 0;
   reco::GenParticle *tWnu1 = 0;
   reco::GenParticle *tWnutau1 = 0;
   reco::GenParticle *tWl1 = 0;
   reco::GenParticle *tWtau1 = 0;
   reco::GenParticle *tWtaunu1 = 0;
   reco::GenParticle *tWtaunutau1 = 0;
   reco::GenParticle *tWtaul1 = 0;
   reco::GenParticle *tWq11 = 0;
   reco::GenParticle *tWq21 = 0;
   
   reco::GenParticle *tW2 = 0;
   reco::GenParticle *tWnu2 = 0;
   reco::GenParticle *tWnutau2 = 0;
   reco::GenParticle *tWl2 = 0;
   reco::GenParticle *tWtau2 = 0;
   reco::GenParticle *tWtaunu2 = 0;
   reco::GenParticle *tWtaunutau2 = 0;
   reco::GenParticle *tWtaul2 = 0;
   reco::GenParticle *tWq12 = 0;
   reco::GenParticle *tWq22 = 0;

   reco::GenParticle *j1 = 0;
   reco::GenParticle *j2 = 0;
   reco::GenParticle *j3 = 0;
   
   int chan = -666;

   // 0   = (t->bW,W->lnu)(t->bW,W->lnu)
   // 1   = (t->bW,W->qq)(t->bW,W->qq)
   // 2   = (t->bW,W->qq)(t->bW,W->lnu)
   // 3   = (t->bW,W->lnu)(t->bW,W->qq)
   // 4   = (t->bW,W->tauLnu)(t->bW,W->tauLnu)
   // 5   = (t->bW,W->qq)(t->bW,W->tauLnu)
   // 6   = (t->bW,W->tauLnu)(t->bW,W->qq)
   // 7   = (t->bW,W->tauHnu)(t->bW,W->tauHnu)
   // 8   = (t->bW,W->qq)(t->bW,W->tauHnu)
   // 9   = (t->bW,W->tauHnu)(t->bW,W->qq)
   // 10  = (t->bW,W->tauHnu)(t->bW,W->tauLnu)
   // 11  = (t->bW,W->tauLnu)(t->bW,W->tauHnu)
   // 12  = (t->bW,W->lnu)(t->bW,W->tauLnu)
   // 13  = (t->bW,W->lnu)(t->bW,W->tauHnu)
   // 14  = (t->bW,W->tauLnu)(t->bW,W->lnu)
   // 15  = (t->bW,W->tauHnu)(t->bW,W->lnu)
   
   // (H->WW:W->lnu,W->lnu)                 +0
   // (H->WW:W->tauLtaunu,W->tauLtaunu)     +20
   // (H->WW:W->tauHtaunu,W->tauHtaunu)     +40
   // (H->WW:W->tauHtaunu,W->tauLtaunu)     +60
   // (H->WW:W->tauLtaunu,W->tauHtaunu)     +80
   // (H->WW:W->qq,W->qq)                   +100
   // (H->WW:W->lnu,W->qq)                  +120
   // (H->WW:W->tauLnutau,W->qq)            +140
   // (H->WW:W->tauHnutau,W->qq)            +160
   // (H->WW:W->qq,W->lnu)                  +180
   // (H->WW:W->qq,W->tauLtaunu)            +200
   // (H->WW:W->qq,W->tauHtaunu)            +220
   // (H->WW:W->lnu,W->tauLtaunu)           +240
   // (H->WW:W->lnu,W->tauHtaunu)           +260
   // (H->WW:W->tauLtaunu,W->lnu)           +280
   // (H->WW:W->tauHtaunu,W->lnu)           +300
   
   // +1000
   // (H->ZZ:Z->ll,Z->ll)                   +0
   // (H->ZZ:Z->tauLtauL,Z->tauLtauL)       +20
   // (H->ZZ:Z->tauLtauL,Z->tauHtauH)       +40
   // (H->ZZ:Z->tauLtauL,Z->tauLtauH)       +60
   // (H->ZZ:Z->tauLtauL,Z->tauHtauL)       +80
   // (H->ZZ:Z->tauHtauH,Z->tauLtauL)       +100
   // (H->ZZ:Z->tauHtauH,Z->tauHtauH)       +120
   // (H->ZZ:Z->tauHtauH,Z->tauLtauH)       +140
   // (H->ZZ:Z->tauHtauH,Z->tauHtauL)       +160
   // (H->ZZ:Z->tauLtauH,Z->tauLtauL)       +180
   // (H->ZZ:Z->tauLtauH,Z->tauHtauH)       +200
   // (H->ZZ:Z->tauLtauH,Z->tauLtauH)       +220
   // (H->ZZ:Z->tauLtauH,Z->tauHtauL)       +240
   // (H->ZZ:Z->tauHtauL,Z->tauLtauL)       +260
   // (H->ZZ:Z->tauHtauL,Z->tauHtauH)       +280
   // (H->ZZ:Z->tauHtauL,Z->tauLtauH)       +300
   // (H->ZZ:Z->tauHtauL,Z->tauHtauL)       +320
   // (H->ZZ:Z->qq,Z->qq)                   +340
   // (H->ZZ:Z->ll,Z->qq)                   +360
   // (H->ZZ:Z->tauHtauH,Z->qq)             +380
   // (H->ZZ:Z->tauLtauL,Z->qq)             +400
   // (H->ZZ:Z->tauHtauL,Z->qq)             +420
   // (H->ZZ:Z->tauLtauH,Z->qq)             +440
   // (H->ZZ:Z->qq,Z->ll)                   +460
   // (H->ZZ:Z->qq,Z->tauHtauH)             +480
   // (H->ZZ:Z->qq,Z->tauLtauL)             +500
   // (H->ZZ:Z->qq,Z->tauHtauL)             +520
   // (H->ZZ:Z->qq,Z->tauLtauH)             +540
   // (H->ZZ:Z->ll,Z->nunu)                 +560
   // (H->ZZ:Z->tauHtauH,Z->nunu)           +580
   // (H->ZZ:Z->tauLtauL,Z->nunu)           +600
   // (H->ZZ:Z->tauHtauL,Z->nunu)           +620
   // (H->ZZ:Z->tauLtauH,Z->nunu)           +640
   // (H->ZZ:Z->qq,Z->nunu)                 +660
   // (H->ZZ:Z->nunu,Z->qq)                 +680
   // (H->ZZ:Z->nunu,Z->ll)                 +700
   // (H->ZZ:Z->nunu,Z->tauHtauH)           +720
   // (H->ZZ:Z->nunu,Z->tauLtauL)           +740
   // (H->ZZ:Z->nunu,Z->tauHtauL)           +760
   // (H->ZZ:Z->nunu,Z->tauLtauH)           +780
   // (H->ZZ:Z->nunu,Z->nunu)               +800
   // (H->ZZ:Z->tauHtauH,Z->ll)             +820
   // (H->ZZ:Z->tauLtauL,Z->ll)             +840
   // (H->ZZ:Z->tauHtauL,Z->ll)             +860
   // (H->ZZ:Z->tauLtauH,Z->ll)             +880
   // (H->ZZ:Z->ll,Z->tauHtauH)             +900
   // (H->ZZ:Z->ll,Z->tauLtauL)             +920
   // (H->ZZ:Z->ll,Z->tauHtauL)             +940
   // (H->ZZ:Z->ll,Z->tauLtauH)             +960
   
   // +2000
   // (H->tautau:tauL,tauL)     +2020
   // (H->tautau:tauH,tauH)     +2040
   // (H->tautau:tauL,tauH)     +2060
   // (H->tautau:tauH,tauL)     +2080

   reco::GenParticleCollection genParticlesCollection = *GenParticles;
   reco::GenParticleCollection::const_iterator genParticleSrc;
   
   int ipart = 0;
   
   for(genParticleSrc = genParticlesCollection.begin();
       genParticleSrc != genParticlesCollection.end(); 
       genParticleSrc++)
     {
	reco::GenParticle *mcp = &(const_cast<reco::GenParticle&>(*genParticleSrc));

	int barcode = ipart; // in CMSSW barcode is the index of genParticle in the event
	// https://twiki.cern.ch/twiki/bin/view/CMS/GenParticles2HepMCConverter
	ipart++;
	
	// Additional partons (up to three)
	if( (fabs(mcp->pdgId()) <= 6 || fabs(mcp->pdgId()) == 21) &&
	    mcp->status() == 23 && barcode == 8 )
	  {
	     if( !j1 )
	       j1 = mcp;
	     else if( !j2 )
	       j2 = mcp;
	     else if( !j3 )
	       j3 = mcp;
	  }	

	// Higgs decays
	if( fabs(mcp->pdgId()) == 25 )
	  {
	     if( (mcp->status() == 62) ||
		 (mcp->status() == 3) )
	       {
		  h0 = const_cast<reco::GenParticle*>(mcp);
		  
		  const reco::GenParticleRefVector& daughterRefs = mcp->daughterRefVector();
		  for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr) 
		    {
		       if( idr->isAvailable() ) 
			 {		       
			    const reco::GenParticleRef& genParticle = (*idr);
			    const reco::GenParticle *d = genParticle.get();
			    reco::GenParticle *pf = getUnique(d,0);

			    // h0 -> bbbar
			    if( fabs(pf->pdgId()) == 5 ) chan = 10000;
			    // h0 -> ee/mumu
			    if( fabs(pf->pdgId()) == 11 || fabs(pf->pdgId()) == 13 ) chan = 10001;
			    // h0 -> gg
			    if( fabs(pf->pdgId()) == 21 ) chan = 10002;
			    // h0 -> gammagamma
			    if( fabs(pf->pdgId()) == 22 ) chan = 10003;
			    // h0 -> qqbar (non-b)
			    if( fabs(pf->pdgId()) == 6 || fabs(pf->pdgId()) <= 4 ) chan = 10004;
			    
			    // h0 -> WW
			    if( fabs(pf->pdgId()) == 24 )
			      {
				 if( h0W1 && !h0W2 ) {h0W2 = pf;}
				 if( !h0W1 ) {h0W1 = pf;}				 
				 
				 if( h0W1 && !h0W2 )
				   {
				      const reco::GenParticleRefVector& daughterRefs = h0W1->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator h0W1_idr = daughterRefs.begin(); 
					  h0W1_idr!= daughterRefs.end(); ++h0W1_idr) 
					{
					   if( h0W1_idr->isAvailable() ) 
					     {		       
						const reco::GenParticleRef& genParticle = (*h0W1_idr);
						const reco::GenParticle *h0W1_d = genParticle.get();
						reco::GenParticle *pff = getUnique(h0W1_d,0);
						
						if( fabs(pff->pdgId()) == 12 ||
						    fabs(pff->pdgId()) == 14 ) // nu
						  {
						     h0Wnu1 = pff;
						  }
						if( fabs(pff->pdgId()) == 16 ) // nutau
						  {
						     h0Wnutau1 = pff;
						  }		
						if( fabs(pff->pdgId()) == 11 ||
						    fabs(pff->pdgId()) == 13 ) // l
						  {
						     h0Wl1 = pff;
						  }		
						if( fabs(pff->pdgId()) == 15 ) // tau
						  {
						     h0Wtau1 = pff;
						     
						     const reco::GenParticleRefVector& daughterRefs = h0Wtau1->daughterRefVector();
						     for(reco::GenParticleRefVector::const_iterator h0Wtau1_idr = daughterRefs.begin();
							 h0Wtau1_idr!= daughterRefs.end(); ++h0Wtau1_idr) 
						       {
							  if( h0Wtau1_idr->isAvailable() ) 
							    {		       
							       const reco::GenParticleRef& genParticle = (*h0Wtau1_idr);
							       const reco::GenParticle *h0Wtau1_d = genParticle.get();
							       reco::GenParticle *pfff = getUnique(h0Wtau1_d,0);
							       
							       if( fabs(pfff->pdgId()) == 12 ||
								   fabs(pfff->pdgId()) == 14 ) // nu
								 {
								    h0Wtaunu1 = pfff;
								 }
							       if( fabs(pfff->pdgId()) == 16 ) // nutau
								 {
								    h0Wtaunutau1 = pfff;
								 }		
							       if( fabs(pfff->pdgId()) == 11 ||
								   fabs(pfff->pdgId()) == 13 ) // l
								 {
								    h0Wtaul1 = pfff;  
								 }
							    }
						       }
						  }						
						if( fabs(pff->pdgId()) <= 6 )
						  {
						     if( h0Wq11 && !h0Wq21 ) {h0Wq21 = pff;}
						     if( !h0Wq11 ) {h0Wq11 = pff;}
						  }
					     }
					}				      
				   }
				 if( h0W2 )
				   {
				      const reco::GenParticleRefVector& daughterRefs = h0W2->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator h0W2_idr = daughterRefs.begin(); 
					  h0W2_idr!= daughterRefs.end(); ++h0W2_idr) 
					{
					   if( h0W2_idr->isAvailable() ) 
					     {		       
						const reco::GenParticleRef& genParticle = (*h0W2_idr);
						const reco::GenParticle *h0W2_d = genParticle.get();
						reco::GenParticle *pff = getUnique(h0W2_d,0);
						
						if( fabs(pff->pdgId()) == 12 ||
						    fabs(pff->pdgId()) == 14 ) // nu
						  {
						     h0Wnu2 = pff;
						  }
						if( fabs(pff->pdgId()) == 16 ) // nutau
						  {
						     h0Wnutau2 = pff;
						  }		
						if( fabs(pff->pdgId()) == 11 ||
						    fabs(pff->pdgId()) == 13 ) // l
						  {
						     h0Wl2 = pff;
						  }		
						if( fabs(pff->pdgId()) == 15 ) // tau
						  {
						     h0Wtau2 = pff;
						     
						     const reco::GenParticleRefVector& daughterRefs = h0Wtau2->daughterRefVector();
						     for(reco::GenParticleRefVector::const_iterator h0Wtau2_idr = daughterRefs.begin();
							 h0Wtau2_idr!= daughterRefs.end(); ++h0Wtau2_idr) 
						       {
							  if( h0Wtau2_idr->isAvailable() ) 
							    {		       
							       const reco::GenParticleRef& genParticle = (*h0Wtau2_idr);
							       const reco::GenParticle *h0Wtau2_d = genParticle.get();
							       reco::GenParticle *pfff = getUnique(h0Wtau2_d,0);
							       
							       if( fabs(pfff->pdgId()) == 12 ||
								   fabs(pfff->pdgId()) == 14 ) // nu
								 {
								    h0Wtaunu2 = pfff;
								 }
							       if( fabs(pfff->pdgId()) == 16 ) // nutau
								 {
								    h0Wtaunutau2 = pfff;
								 }		
							       if( fabs(pfff->pdgId()) == 11 ||
								   fabs(pfff->pdgId()) == 13 ) // l
								 {
								    h0Wtaul2 = pfff;  
								 }
							    }
						       }
						  }						
						if( fabs(pff->pdgId()) <= 6 )
						  {
						     if( h0Wq12 && !h0Wq22 ) {h0Wq22 = pff;}
						     if( !h0Wq12 ) {h0Wq12 = pff;}
						  }
					     }
					}				      				      
				   }				 
			      }

			    // h0 -> ZZ
			    if( fabs(pf->pdgId()) == 23 )
			      {
				 if( h0Z1 && !h0Z2 ) {h0Z2 = pf;}
				 if( !h0Z1 ) {h0Z1 = pf;}
				 
				 if( h0Z1 && !h0Z2 )
				   {
				      const reco::GenParticleRefVector& daughterRefs = h0Z1->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator h0Z1_idr = daughterRefs.begin(); 
					  h0Z1_idr!= daughterRefs.end(); ++h0Z1_idr) 
					{
					   if( h0Z1_idr->isAvailable() ) 
					     {		       
						const reco::GenParticleRef& genParticle = (*h0Z1_idr);
						const reco::GenParticle *h0Z1_d = genParticle.get();
						reco::GenParticle *pff = getUnique(h0Z1_d,0);
						
						if( fabs(pff->pdgId()) == 11 ||
						    fabs(pff->pdgId()) == 13 ) // l
						  {
						     if( h0Zl11 && !h0Zl21 ) {h0Zl21 = pff;}
						     if( !h0Zl11 ) {h0Zl11 = pff;}
						  }				
						if( fabs(pff->pdgId()) == 15 ) // tau
						  {
						     if( h0Ztau11 && !h0Ztau21 )
						       {
							  h0Ztau21 = pff;
							  
							  const reco::GenParticleRefVector& daughterRefs = h0Ztau21->daughterRefVector();
							  for(reco::GenParticleRefVector::const_iterator h0Ztau21_idr = daughterRefs.begin(); 
							      h0Ztau21_idr!= daughterRefs.end(); ++h0Ztau21_idr) 
							    {
							       if( h0Ztau21_idr->isAvailable() ) 
								 {		       
								    const reco::GenParticleRef& genParticle = (*h0Ztau21_idr);
								    const reco::GenParticle *h0Ztau21_d = genParticle.get();
								    reco::GenParticle *pfff = getUnique(h0Ztau21_d,0);
								    
								    if( fabs(pfff->pdgId()) == 12 ||
									fabs(pfff->pdgId()) == 14 ) // nu
								      {
									 h0Ztaunu21 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 16 ) // nutau
								      {
									 h0Ztaunutau21 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 11 ||
									fabs(pfff->pdgId()) == 13 ) // l
								      {
									 h0Ztaul21 = pfff;
								      }
								 }
							    }							  
						       }
						     if( !h0Ztau11 )
						       {
							  h0Ztau11 = pff;

							  const reco::GenParticleRefVector& daughterRefs = h0Ztau11->daughterRefVector();
							  for(reco::GenParticleRefVector::const_iterator h0Ztau11_idr = daughterRefs.begin(); 
							      h0Ztau11_idr!= daughterRefs.end(); ++h0Ztau11_idr) 
							    {
							       if( h0Ztau11_idr->isAvailable() ) 
								 {		       
								    const reco::GenParticleRef& genParticle = (*h0Ztau11_idr);
								    const reco::GenParticle *h0Ztau11_d = genParticle.get();
								    reco::GenParticle *pfff = getUnique(h0Ztau11_d,0);
								    
								    if( fabs(pfff->pdgId()) == 12 ||
									fabs(pfff->pdgId()) == 14 ) // nu
								      {
									 h0Ztaunu11 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 16 ) // nutau
								      {
									 h0Ztaunutau11 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 11 ||
									fabs(pfff->pdgId()) == 13 ) // l
								      {
									 h0Ztaul11 = pfff;
								      }
								 }
							    }							  
						       }						     
						  }
						if( fabs(pff->pdgId()) <= 6 ) // q
						  {
						     if( h0Zq11 && !h0Zq21 ) {h0Zq21 = pff;}
						     if( !h0Zq11 ) {h0Zq11 = pff;}
						  }				
						if( fabs(pff->pdgId()) == 12 ||
						    fabs(pff->pdgId()) == 14 ||
						    fabs(pff->pdgId()) == 16 ) // nu
						  {
						     if( h0Znu11 && !h0Znu21 ) {h0Znu21 = pff;}
						     if( !h0Znu11 ) {h0Znu11 = pff;}
						  }						
					     }					   
					}				      
				   }
				 if( h0Z2 )
				   {
				      const reco::GenParticleRefVector& daughterRefs = h0Z2->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator h0Z2_idr = daughterRefs.begin();
					  h0Z2_idr!= daughterRefs.end(); ++h0Z2_idr) 
					{
					   if( h0Z2_idr->isAvailable() ) 
					     {		       
						const reco::GenParticleRef& genParticle = (*h0Z2_idr);
						const reco::GenParticle *h0Z2_d = genParticle.get();
						reco::GenParticle *pff = getUnique(h0Z2_d,0);
						
						if( fabs(pff->pdgId()) == 11 ||
						    fabs(pff->pdgId()) == 13 ) // l
						  {
						     if( h0Zl12 && !h0Zl22 ) {h0Zl22 = pff;}
						     if( !h0Zl12 ) {h0Zl12 = pff;}
						  }				
						if( fabs(pff->pdgId()) == 15 ) // tau
						  {
						     if( h0Ztau12 && !h0Ztau22 )
						       {
							  h0Ztau22 = pff;
							  
							  const reco::GenParticleRefVector& daughterRefs = h0Ztau22->daughterRefVector();
							  for(reco::GenParticleRefVector::const_iterator h0Ztau22_idr = daughterRefs.begin(); 
							      h0Ztau22_idr!= daughterRefs.end(); ++h0Ztau22_idr) 
							    {
							       if( h0Ztau22_idr->isAvailable() ) 
								 {		       
								    const reco::GenParticleRef& genParticle = (*h0Ztau22_idr);
								    const reco::GenParticle *h0Ztau22_d = genParticle.get();
								    reco::GenParticle *pfff = getUnique(h0Ztau22_d,0);
								    
								    if( fabs(pfff->pdgId()) == 12 ||
									fabs(pfff->pdgId()) == 14 ) // nu
								      {
									 h0Ztaunu22 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 16 ) // nutau
								      {
									 h0Ztaunutau22 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 11 ||
									fabs(pfff->pdgId()) == 13 ) // l
								      {
									 h0Ztaul22 = pfff;
								      }
								 }
							    }							  
						       }
						     if( !h0Ztau12 )
						       {
							  h0Ztau12 = pff;

							  const reco::GenParticleRefVector& daughterRefs = h0Ztau12->daughterRefVector();
							  for(reco::GenParticleRefVector::const_iterator h0Ztau12_idr = daughterRefs.begin(); 
							      h0Ztau12_idr!= daughterRefs.end(); ++h0Ztau12_idr)
							    {
							       if( h0Ztau12_idr->isAvailable() ) 
								 {		       
								    const reco::GenParticleRef& genParticle = (*h0Ztau12_idr);
								    const reco::GenParticle *h0Ztau12_d = genParticle.get();
								    reco::GenParticle *pfff = getUnique(h0Ztau12_d,0);
								    
								    if( fabs(pfff->pdgId()) == 12 ||
									fabs(pfff->pdgId()) == 14 ) // nu
								      {
									 h0Ztaunu12 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 16 ) // nutau
								      {
									 h0Ztaunutau12 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 11 ||
									fabs(pfff->pdgId()) == 13 ) // l
								      {
									 h0Ztaul12 = pfff;
								      }
								 }
							    }							  
						       }						     
						  }
						if( fabs(pff->pdgId()) <= 6 ) // q
						  {
						     if( h0Zq12 && !h0Zq22 ) {h0Zq22 = pff;}
						     if( !h0Zq12 ) {h0Zq12 = pff;}
						  }				
						if( fabs(pff->pdgId()) == 12 ||
						    fabs(pff->pdgId()) == 14 ||
						    fabs(pff->pdgId()) == 16 ) // nu
						  {
						     if( h0Znu12 && !h0Znu22 ) {h0Znu22 = pff;}
						     if( !h0Znu12 ) {h0Znu12 = pff;}
						  }						
					     }					   
					}				      				      
				   }				 
			      }			    
			    
			    // h0 -> tautau
			    if( fabs(pf->pdgId()) == 15 && pf->status() == 2 ) // tau to decay
			      {
				 if( h0tau1 && !h0tau2 ) {h0tau2 = pf;}
				 if( !h0tau1 ) {h0tau1 = pf;}
				 
				 if( h0tau1 && !h0tau2 )
				   {
				      const reco::GenParticleRefVector& daughterRefs = h0tau1->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator h0tau1_idr = daughterRefs.begin();
					  h0tau1_idr!= daughterRefs.end(); ++h0tau1_idr) 
					{
					   if( h0tau1_idr->isAvailable() ) 
					     {		       
						const reco::GenParticleRef& genParticle = (*h0tau1_idr);
						const reco::GenParticle *h0tau1_d = genParticle.get();
						reco::GenParticle *pff = getUnique(h0tau1_d,0);
						
						if( fabs(pff->pdgId()) == 11 ||
						    fabs(pff->pdgId()) == 13 ||
						    fabs(pff->pdgId()) == 15 ) // l
						  {
						     h0taul1 = pff;
						  }		
						if( fabs(pff->pdgId()) == 16 ) // nu_tau
						  {
						     h0taunutau1 = pff;
						  }		
						if( fabs(pff->pdgId()) == 12 ||
						    fabs(pff->pdgId()) == 14 ) // nu_e or nu_mu
						  {
						     h0taunu1 = pff;
						  }						
					     }
					}				      				      
				   }
				 if( h0tau2 )
				   {
				      const reco::GenParticleRefVector& daughterRefs = h0tau2->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator h0tau2_idr = daughterRefs.begin();
					  h0tau2_idr!= daughterRefs.end(); ++h0tau2_idr) 
					{
					   if( h0tau2_idr->isAvailable() )
					     {		       
						const reco::GenParticleRef& genParticle = (*h0tau2_idr);
						const reco::GenParticle *h0tau2_d = genParticle.get();
						reco::GenParticle *pff = getUnique(h0tau2_d,0);
						
						if( fabs(pff->pdgId()) == 11 ||
						    fabs(pff->pdgId()) == 13 ||
						    fabs(pff->pdgId()) == 15 ) // l
						  {
						     h0taul2 = pff;
						  }		
						if( fabs(pff->pdgId()) == 16 ) // nu_tau
						  {
						     h0taunutau2 = pff;
						  }		
						if( fabs(pff->pdgId()) == 12 ||
						    fabs(pff->pdgId()) == 14 ) // nu_e or nu_mu
						  {
						     h0taunu2 = pff;
						  }						
					     }
					}				      				      
				   }				 
			      }			    
			 }		       
		    }		  
	       }	     
	  }	
	
	// top decays
	if( fabs(mcp->pdgId()) == 6
	    && ( (mcp->status() == 62) || 
		 (mcp->status() == 3)
	       ) )
	  {
	     if( t1 && !t2 ) {t2 = const_cast<reco::GenParticle*>(mcp);}
	     if( !t1 ) {t1 = const_cast<reco::GenParticle*>(mcp);}

	     const reco::GenParticleRefVector& daughterRefs = mcp->daughterRefVector();
	     for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr) 
	       {
		  if( idr->isAvailable() ) 
		    {		       
		       const reco::GenParticleRef& genParticle = (*idr);
		       const reco::GenParticle *d = genParticle.get();
		       reco::GenParticle *pf = getUnique(d,0);
		       
//		       if( pf->status() != 3 && pf->status() != 62 ) continue;
		       
		       if( fabs(pf->pdgId()) == 5 || fabs(pf->pdgId()) == 3 || fabs(pf->pdgId()) == 1 ) // b or s or d
			 {
			    if( tb1 && !tb2 ) {tb2 = pf;}
			    if( !tb1 ) {tb1 = pf;}			    
			 }		       
		       
		       if( fabs(pf->pdgId()) == 24 ) // W
			 {
			    if( tW1 && !tW2 )
			      {
				 tW2 = pf;
				 const reco::GenParticleRefVector& tW2_daughterRefs = tW2->daughterRefVector();
				 for(reco::GenParticleRefVector::const_iterator tW2_idr = tW2_daughterRefs.begin();
				     tW2_idr!= tW2_daughterRefs.end(); ++tW2_idr) 
				   {
				      if( tW2_idr->isAvailable() ) 
					{		       
					   const reco::GenParticleRef& tW2_genParticle = (*tW2_idr);
					   const reco::GenParticle *tW2_d = tW2_genParticle.get();
					   reco::GenParticle *pff = getUnique(tW2_d,0);
					   
					   if( fabs(pff->pdgId()) == 12 ||
					       fabs(pff->pdgId()) == 14 ) // nu
					     {
						tWnu2 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 16 ) // nu_tau
					     {
						tWnutau2 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 11 ||
					       fabs(pff->pdgId()) == 13 ) // l
					     {
						tWl2 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 15 ) // tau
					     {
						tWtau2 = pff;
						
						const reco::GenParticleRefVector& tWtau2_daughterRefs = tWtau2->daughterRefVector();
						for(reco::GenParticleRefVector::const_iterator tWtau2_idr = tWtau2_daughterRefs.begin();
						    tWtau2_idr!= tWtau2_daughterRefs.end(); ++tWtau2_idr) 
						  {
						     if( tWtau2_idr->isAvailable() ) 
						       {		       
							  const reco::GenParticleRef& tWtau2_genParticle = (*tWtau2_idr);
							  const reco::GenParticle *tWtau2_d = tWtau2_genParticle.get();
							  reco::GenParticle *pfff = getUnique(tWtau2_d,0);
							  
							  if( fabs(pfff->pdgId()) == 12 ||
							      fabs(pfff->pdgId()) == 14 ) // nu
							    {
							       tWtaunu2 = pfff;
							    }		
							  if( fabs(pfff->pdgId()) == 16 ) // nu_tau
							    {
							       tWtaunutau2 = pfff;
							    }			
							  if( fabs(pfff->pdgId()) == 11 ||
							      fabs(pfff->pdgId()) == 13 ) // l
							    {
							       tWtaul2 = pfff;
							    }							  
						       }
						  }						
					     }
					   if( fabs(pff->pdgId()) <= 6 ) // q
					     {
						if( tWq12 && !tWq22 ) {tWq22 = pff;}
						if( !tWq12 ) {tWq12 = pff;}
					     }					   					   
					}
				   }				
			      }
			    
			    if( !tW1 )
			      {
				 tW1 = pf;
				 const reco::GenParticleRefVector& tW1_daughterRefs = tW1->daughterRefVector();
				 for(reco::GenParticleRefVector::const_iterator tW1_idr = tW1_daughterRefs.begin();
				     tW1_idr!= tW1_daughterRefs.end(); ++tW1_idr) 
				   {
				      if( tW1_idr->isAvailable() ) 
					{		       
					   const reco::GenParticleRef& tW1_genParticle = (*tW1_idr);
					   const reco::GenParticle *tW1_d = tW1_genParticle.get();
					   reco::GenParticle *pff = getUnique(tW1_d,0);
					   
					   if( fabs(pff->pdgId()) == 12 ||
					       fabs(pff->pdgId()) == 14 ) // nu
					     {
						tWnu1 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 16 ) // nu_tau
					     {
						tWnutau1 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 11 ||
					       fabs(pff->pdgId()) == 13 ) // l
					     {
						tWl1 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 15 ) // tau
					     {
						tWtau1 = pff;
						
						const reco::GenParticleRefVector& tWtau1_daughterRefs = tWtau1->daughterRefVector();
						for(reco::GenParticleRefVector::const_iterator tWtau1_idr = tWtau1_daughterRefs.begin();
						    tWtau1_idr!= tWtau1_daughterRefs.end(); ++tWtau1_idr) 
						  {
						     if( tWtau1_idr->isAvailable() ) 
						       {		       
							  const reco::GenParticleRef& tWtau1_genParticle = (*tWtau1_idr);
							  const reco::GenParticle *tWtau1_d = tWtau1_genParticle.get();
							  reco::GenParticle *pfff = getUnique(tWtau1_d,0);
							  
							  if( fabs(pfff->pdgId()) == 12 ||
							      fabs(pfff->pdgId()) == 14 ) // nu
							    {
							       tWtaunu1 = pfff;
							    }		
							  if( fabs(pfff->pdgId()) == 16 ) // nu_tau
							    {
							       tWtaunutau1 = pfff;
							    }			
							  if( fabs(pfff->pdgId()) == 11 ||
							      fabs(pfff->pdgId()) == 13 ) // l
							    {
							       tWtaul1 = pfff;
							    }							  
						       }
						  }						
					     }
					   if( fabs(pff->pdgId()) <= 6 ) // q
					     {
						if( tWq11 && !tWq21 ) {tWq21 = pff;}
						if( !tWq11 ) {tWq11 = pff;}
					     }					   					   
					}
				   }				
			      }			    
			 }		       
		    }
	       }	     
	  }
     }

   bool doCheck = 0;

   if( h0 && t1 && t2 && tb1 && tb2 && tW1 && tW2 )
     {	
	int tchan = -666;
	if( tWl1 && tWl2 )   tchan = 0;
	if( tWq11 && tWq12 ) tchan = 1;
	if( tWq11 && tWl2 )  tchan = 2;
	if( tWl1 && tWq12 )  tchan = 3;
	if( tWtaul1 && tWtaul2 )  tchan = 4;
	if( tWq11 && tWtaul2 )  tchan = 5;
	if( tWtaul1 && tWq12 )  tchan = 6;
	if( tWtaunutau1 && tWtaunutau2 && !tWtaul1 && !tWtaul2 )  tchan = 7;
	if( tWq11 && tWtaunutau2 && !tWtaul2 )  tchan = 8;
	if( tWtaunutau1 && !tWtaul1 && tWq12 )  tchan = 9;
	if( tWtaunutau1 && tWtaul2 && !tWtaul1 )  tchan = 10;
	if( tWtaul1 && tWtaunutau2 && !tWtaul2 )  tchan = 11;
	if( tWl1 && tWtaul2 )  tchan = 12;
	if( tWl1 && !tWtaul2 && tWtaunutau2 )  tchan = 13;
	if( tWl2 && tWtaul1 )  tchan = 14;
	if( tWl2 && !tWtaul1 && tWtaunutau1 )  tchan = 15;
	
	if( tchan < 0 && doCheck )
	  {	     
	     std::cout << "Failed to identify top-quark decay chain" << std::endl;
	     
	     std::cout << "t1 = " << bool(t1) << std::endl;
	     std::cout << "t1->W = " << bool(tW1) << std::endl;
	     std::cout << "t1->W->l = " << bool(tWl1) << std::endl;
	     std::cout << "t1->W->nu = " << bool(tWnu1) << std::endl;
	     std::cout << "t1->W->nutau = " << bool(tWnutau1) << std::endl;
	     std::cout << "t1->W->tau = " << bool(tWtau1) << std::endl;
	     std::cout << "t1->W->tau->l = " << bool(tWtaul1) << std::endl;
	     std::cout << "t1->W->tau->nu = " << bool(tWtaunu1) << std::endl;
	     std::cout << "t1->W->tau->nutau = " << bool(tWtaunutau1) << std::endl;
	     std::cout << "t1->W->q = " << bool(tWq11) << std::endl;
	             
	     std::cout << "t2 = " << bool(t2) << std::endl;
	     std::cout << "t2->W = " << bool(tW2) << std::endl;
	     std::cout << "t2->W->l = " << bool(tWl2) << std::endl;
	     std::cout << "t2->W->nu = " << bool(tWnu2) << std::endl;
	     std::cout << "t2->W->nutau = " << bool(tWnutau2) << std::endl;
	     std::cout << "t2->W->tau = " << bool(tWtau2) << std::endl;
	     std::cout << "t2->W->tau->l = " << bool(tWtaul2) << std::endl;
	     std::cout << "t2->W->tau->nu = " << bool(tWtaunu2) << std::endl;
	     std::cout << "t2->W->tau->nutau = " << bool(tWtaunutau2) << std::endl;
	     std::cout << "t2->W->q = " << bool(tWq12) << std::endl;
	             
	     exit(1);
	  }
	
	if( h0W1 && h0W2 )
	  {	              
	     int chan0 = 0;
	     if( h0Wl1 && h0Wl2 ) chan = chan0 + 0 + tchan;
	     if( h0Wtaul1 && h0Wtaul2 ) chan = chan0 + 20 + tchan;
	     if( h0Wtaunutau1 && !h0Wtaul1 && h0Wtaunutau2 && !h0Wtaul2 ) chan = chan0 + 40 + tchan;
	     if( h0Wtaunutau1 && !h0Wtaul1 && h0Wtaul2 ) chan = chan0 + 60 + tchan;
	     if( h0Wtaul1 && h0Wtaunutau2 && !h0Wtaul2 ) chan = chan0 + 80 + tchan;
	     if( h0Wq11 && h0Wq12 ) chan = chan0 + 100 + tchan;
	     if( h0Wl1 && h0Wq12 ) chan = chan0 + 120 + tchan;
	     if( h0Wtaul1 && h0Wq12 ) chan = chan0 + 140 + tchan;
	     if( h0Wtaunutau1 && !h0Wtaul1 && h0Wq12 ) chan = chan0 + 160 + tchan;
	     if( h0Wq11 && h0Wl2 ) chan = chan0 + 180 + tchan;
	     if( h0Wq11 && h0Wtaul2 ) chan = chan0 + 200 + tchan;
	     if( h0Wq11 && h0Wtaunutau2 && !h0Wtaul2 ) chan = chan0 + 220 + tchan;
	     if( h0Wl1 && h0Wtaul2 ) chan = chan0 + 240 + tchan;
	     if( h0Wl1 && !h0Wtaul2 && h0Wtaunutau2 ) chan = chan0 + 260 + tchan;
	     if( h0Wl2 && h0Wtaul1 ) chan = chan0 + 280 + tchan;
	     if( h0Wl2 && !h0Wtaul1 && h0Wtaunutau1 ) chan = chan0 + 300 + tchan;
	  }

	if( h0Z1 && h0Z2 )
	  {
	     int chan0 = 1000;
	     if( h0Zl11 && h0Zl12 ) chan = chan0 + 0 + tchan;
	     if( h0Ztaul11 && h0Ztaul21 && h0Ztaul12 && h0Ztaul22 ) chan = chan0 + 20 + tchan;
	     if( h0Ztaul11 && h0Ztaul21 && !h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau12 && h0Ztaunutau22 ) chan = chan0 + 40 + tchan;
	     if( h0Ztaul11 && h0Ztaul21 && h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau22 ) chan = chan0 + 60 + tchan;
	     if( h0Ztaul11 && h0Ztaul21 && h0Ztaul22 && !h0Ztaul12 && h0Ztaunutau12 ) chan = chan0 + 80 + tchan;
	     if( !h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau11 && h0Ztaunutau21 && h0Ztaul12 && h0Ztaul22 ) chan = chan0 + 100 + tchan;
	     if( !h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau11 && h0Ztaunutau21 && !h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau12 && h0Ztaunutau22 ) chan = chan0 + 120 + tchan;
	     if( !h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau11 && h0Ztaunutau21 && h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau22 ) chan = chan0 + 140 + tchan;
	     if( !h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau11 && h0Ztaunutau21 && !h0Ztaul12 && h0Ztaul22 && h0Ztaunutau12 ) chan = chan0 + 160 + tchan;
	     if( h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau21 && h0Ztaul12 && h0Ztaul22 ) chan = chan0 + 180 + tchan;
	     if( h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau21 && !h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau12 && h0Ztaunutau22 ) chan = chan0 + 200 + tchan;
	     if( h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau21 && h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau22 ) chan = chan0 + 220 + tchan;
	     if( h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau21 && !h0Ztaul12 && h0Ztaul22 && h0Ztaunutau12 ) chan = chan0 + 240 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && h0Ztaul21 && h0Ztaul12 && h0Ztaul22 ) chan = chan0 + 260 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && h0Ztaul21 && !h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau12 && h0Ztaunutau22 ) chan = chan0 + 280 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && h0Ztaul21 && h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau22 ) chan = chan0 + 300 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && h0Ztaul21 && !h0Ztaul12 && h0Ztaul22 && h0Ztaunutau12 ) chan = chan0 + 320 + tchan;
	     if( h0Zq11 && h0Zq12 ) chan = chan0 + 340 + tchan;
	     if( h0Zl11 && h0Zq12 ) chan = chan0 + 360 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && !h0Ztaul21 && h0Ztaunutau21 && h0Zq12 ) chan = chan0 + 380 + tchan;
	     if( h0Ztaul11 && h0Ztaul21 && h0Zq12 ) chan = chan0 + 400 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && h0Ztaul21 && h0Zq12 ) chan = chan0 + 420 + tchan;
	     if( h0Ztaul11 && h0Ztaunutau21 && !h0Ztaul21 && h0Zq12 ) chan = chan0 + 440 + tchan;
	     if( h0Zq11 && h0Zl12 ) chan = chan0 + 460 + tchan;
	     if( h0Zq11 && !h0Ztaul12 && h0Ztaunutau12 && !h0Ztaul22 && h0Ztaunutau22 ) chan = chan0 + 480 + tchan;
	     if( h0Zq11 && h0Ztaul12 && h0Ztaul22 ) chan = chan0 + 500 + tchan;
	     if( h0Zq11 && !h0Ztaul12 && h0Ztaunutau12 && h0Ztaul22 ) chan = chan0 + 520 + tchan;
	     if( h0Zq11 && h0Ztaul12 && h0Ztaunutau22 && !h0Ztaul22 ) chan = chan0 + 540 + tchan;
	     if( h0Zl11 && h0Znu12 ) chan = chan0 + 560 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && !h0Ztaul21 && h0Ztaunutau21 && h0Znu12 ) chan = chan0 + 580 + tchan;
	     if( h0Ztaul11 && h0Ztaul21 && h0Znu12 ) chan = chan0 + 600 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && h0Ztaul21 && h0Znu12 ) chan = chan0 + 620 + tchan;
	     if( h0Ztaul11 && h0Ztaunutau21 && !h0Ztaul21 && h0Znu12 ) chan = chan0 + 640 + tchan;
	     if( h0Zq11 && h0Znu12 ) chan = chan0 + 660 + tchan;
	     if( h0Znu11 && h0Zq12 ) chan = chan0 + 680 + tchan;
	     if( h0Znu11 && h0Zl12 ) chan = chan0 + 700 + tchan;
	     if( h0Znu11 && !h0Ztaul12 && h0Ztaunutau12 && !h0Ztaul22 && h0Ztaunutau22 ) chan = chan0 + 720 + tchan;
	     if( h0Znu11 && h0Ztaul12 && h0Ztaul22 ) chan = chan0 + 740 + tchan;
	     if( h0Znu11 && !h0Ztaul12 && h0Ztaunutau12 && h0Ztaul22 ) chan = chan0 + 760 + tchan;
	     if( h0Znu11 && h0Ztaul12 && h0Ztaunutau22 && !h0Ztaul22 ) chan = chan0 + 780 + tchan;
	     if( h0Znu11 && h0Znu12 ) chan = chan0 + 800 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && !h0Ztaul21 && h0Ztaunutau21 && h0Zl12 ) chan = chan0 + 820 + tchan;
	     if( h0Ztaul11 && h0Ztaul21 && h0Zl12 ) chan = chan0 + 840 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && h0Ztaul21 && h0Zl12 ) chan = chan0 + 860 + tchan;
	     if( h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau21 && h0Zl12 ) chan = chan0 + 880 + tchan;
	     if( h0Zl11 && !h0Ztaul12 && h0Ztaunutau12 && !h0Ztaul22 && h0Ztaunutau22 ) chan = chan0 + 900 + tchan;
	     if( h0Zl11 && h0Ztaul12 && h0Ztaul22 ) chan = chan0 + 920 + tchan;
	     if( h0Zl11 && !h0Ztaul12 && h0Ztaunutau12 && h0Ztaul22 ) chan = chan0 + 940 + tchan;
	     if( h0Zl11 && h0Ztaul12 && h0Ztaunutau22 && !h0Ztaul22 ) chan = chan0 + 960 + tchan;	     
	  }	
	
	if( h0tau1 && h0tau2 )
	  {	     
	     int chan0 = 2000;
	     if( h0taul1 && h0taul2 && h0taunu1 && h0taunu2 ) chan = chan0 + 20 + tchan;
	     if( !h0taul1 && !h0taul2 && h0taunutau1 && h0taunutau2 ) chan = chan0 + 40 + tchan;
	     if( h0taul1 && !h0taul2 && h0taunutau2 && h0taunu1 ) chan = chan0 + 60 + tchan;
	     if( h0taul2 && !h0taul1 && h0taunutau1 && h0taunu2 ) chan = chan0 + 80 + tchan;
	  }	

	if( chan < 0 && doCheck )
	  {
	     std::cout << "Unknown channel found" << std::endl;
	     std::cout << "chan = " << chan << std::endl;

	     std::cout << "j1 = " << bool(j1) << std::endl;
	     std::cout << "j2 = " << bool(j2) << std::endl;
	     std::cout << "j3 = " << bool(j3) << std::endl;
	     
	     std::cout << "h0 = " << bool(h0) << std::endl;
	     std::cout << "h0->W1 = " << bool(h0W1) << std::endl;
	     std::cout << "h0->W1->l = " << bool(h0Wl1) << std::endl;
	     std::cout << "h0->W1->nu = " << bool(h0Wnu1) << std::endl;
	     std::cout << "h0->W1->tau = " << bool(h0Wtau1) << std::endl;
	     std::cout << "h0->W1->nutau = " << bool(h0Wnutau1) << std::endl;
	     std::cout << "h0->W1->tau->l = " << bool(h0Wtaul1) << std::endl;
	     std::cout << "h0->W1->tau->nu = " << bool(h0Wtaunu1) << std::endl;
	     std::cout << "h0->W1->tau->nutau = " << bool(h0Wtaunutau1) << std::endl;
	     std::cout << "h0->W1->q1 = " << bool(h0Wq11) << std::endl;
	     std::cout << "h0->W1->q2 = " << bool(h0Wq21) << std::endl;
	     std::cout << "h0->W2 = " << bool(h0W2) << std::endl;
	     std::cout << "h0->W2->l = " << bool(h0Wl2) << std::endl;
	     std::cout << "h0->W2->nu = " << bool(h0Wnu2) << std::endl;
	     std::cout << "h0->W2->tau = " << bool(h0Wtau2) << std::endl;
	     std::cout << "h0->W2->nutau = " << bool(h0Wnutau2) << std::endl;
	     std::cout << "h0->W2->tau->l = " << bool(h0Wtaul2) << std::endl;
	     std::cout << "h0->W2->tau->nu = " << bool(h0Wtaunu2) << std::endl;
	     std::cout << "h0->W2->tau->nutau = " << bool(h0Wtaunutau2) << std::endl;
	     std::cout << "h0->W2->q1 = " << bool(h0Wq12) << std::endl;
	     std::cout << "h0->W2->q2 = " << bool(h0Wq22) << std::endl;
	     
	     std::cout << "h0->Z1 = " << bool(h0Z1) << std::endl;
	     std::cout << "h0->Z1->l1 = " << bool(h0Zl11) << std::endl;
	     std::cout << "h0->Z1->l2 = " << bool(h0Zl21) << std::endl;
	     std::cout << "h0->Z1->tau1 = " << bool(h0Ztau11) << std::endl;
	     std::cout << "h0->Z1->tau1->l = " << bool(h0Ztaul11) << std::endl;
	     std::cout << "h0->Z1->tau1->nu = " << bool(h0Ztaunu11) << std::endl;
	     std::cout << "h0->Z1->tau1->nutau = " << bool(h0Ztaunutau11) << std::endl;
	     std::cout << "h0->Z1->tau2 = " << bool(h0Ztau21) << std::endl;
	     std::cout << "h0->Z1->tau2->l = " << bool(h0Ztaul21) << std::endl;
	     std::cout << "h0->Z1->tau2->nu = " << bool(h0Ztaunu21) << std::endl;
	     std::cout << "h0->Z1->tau2->nutau = " << bool(h0Ztaunutau21) << std::endl;
	     std::cout << "h0->Z1->q1 = " << bool(h0Zq11) << std::endl;
	     std::cout << "h0->Z1->q2 = " << bool(h0Zq21) << std::endl;
	     std::cout << "h0->Z1->nu1 = " << bool(h0Znu11) << std::endl;
	     std::cout << "h0->Z1->nu2 = " << bool(h0Znu21) << std::endl;
	     
	     std::cout << "h0->Z2 = " << bool(h0Z2) << std::endl;
	     std::cout << "h0->Z2->l1 = " << bool(h0Zl12) << std::endl;
	     std::cout << "h0->Z2->l2 = " << bool(h0Zl22) << std::endl;
	     std::cout << "h0->Z2->tau1 = " << bool(h0Ztau12) << std::endl;
	     std::cout << "h0->Z2->tau1->l = " << bool(h0Ztaul12) << std::endl;
	     std::cout << "h0->Z2->tau1->nu = " << bool(h0Ztaunu12) << std::endl;
	     std::cout << "h0->Z2->tau1->nutau = " << bool(h0Ztaunutau12) << std::endl;
	     std::cout << "h0->Z2->tau2 = " << bool(h0Ztau22) << std::endl;
	     std::cout << "h0->Z2->tau2->l = " << bool(h0Ztaul22) << std::endl;
	     std::cout << "h0->Z2->tau2->nu = " << bool(h0Ztaunu22) << std::endl;
	     std::cout << "h0->Z2->tau2->nutau = " << bool(h0Ztaunutau22) << std::endl;
	     std::cout << "h0->Z2->q1 = " << bool(h0Zq12) << std::endl;
	     std::cout << "h0->Z2->q2 = " << bool(h0Zq22) << std::endl;
	     std::cout << "h0->Z2->nu1 = " << bool(h0Znu12) << std::endl;
	     std::cout << "h0->Z2->nu2 = " << bool(h0Znu22) << std::endl;
	     
	     std::cout << "h0->tau1 = " << bool(h0tau1) << std::endl;
	     std::cout << "h0->tau1->l = " << bool(h0taul1) << std::endl;
	     std::cout << "h0->tau1->nu_tau = " << bool(h0taunutau1) << std::endl;
	     std::cout << "h0->tau1->nu = " << bool(h0taunu1) << std::endl;
	             
	     std::cout << "h0->tau2 = " << bool(h0tau2) << std::endl;
	     std::cout << "h0->tau2->l = " << bool(h0taul2) << std::endl;
	     std::cout << "h0->tau2->nu_tau = " << bool(h0taunutau2) << std::endl;
	     std::cout << "h0->tau2->nu = " << bool(h0taunu2) << std::endl;
	     
	     exit(1);
	  }
     }

   tree.mc_truth_tth_channel = chan;

   // TLV

   if( h0 ) p4toTLV(h0->p4(),tree.mc_truth_tth_h0_p4);      
   
   if( h0W1 ) p4toTLV(h0W1->p4(),tree.mc_truth_tth_h0W1_p4);
   if( h0W2 ) p4toTLV(h0W2->p4(),tree.mc_truth_tth_h0W2_p4);
   if( h0Wl1 ) p4toTLV(h0Wl1->p4(),tree.mc_truth_tth_h0Wl1_p4);
   if( h0Wnu1 ) p4toTLV(h0Wnu1->p4(),tree.mc_truth_tth_h0Wnu1_p4);
   if( h0Wtau1 ) p4toTLV(h0Wtau1->p4(),tree.mc_truth_tth_h0Wtau1_p4);
   if( h0Wnutau1 ) p4toTLV(h0Wnutau1->p4(),tree.mc_truth_tth_h0Wnutau1_p4);
   if( h0Wtaul1 ) p4toTLV(h0Wtaul1->p4(),tree.mc_truth_tth_h0Wtaul1_p4);
   if( h0Wtaunu1 ) p4toTLV(h0Wtaunu1->p4(),tree.mc_truth_tth_h0Wtaunu1_p4);
   if( h0Wtaunutau1 ) p4toTLV(h0Wtaunutau1->p4(),tree.mc_truth_tth_h0Wtaunutau1_p4);
   if( h0Wl2 ) p4toTLV(h0Wl2->p4(),tree.mc_truth_tth_h0Wl2_p4);
   if( h0Wnu2 ) p4toTLV(h0Wnu2->p4(),tree.mc_truth_tth_h0Wnu2_p4);
   if( h0Wtau2 ) p4toTLV(h0Wtau2->p4(),tree.mc_truth_tth_h0Wtau2_p4);
   if( h0Wnutau2 ) p4toTLV(h0Wnutau2->p4(),tree.mc_truth_tth_h0Wnutau2_p4);
   if( h0Wtaul2 ) p4toTLV(h0Wtaul2->p4(),tree.mc_truth_tth_h0Wtaul2_p4);
   if( h0Wtaunu2 ) p4toTLV(h0Wtaunu2->p4(),tree.mc_truth_tth_h0Wtaunu2_p4);
   if( h0Wtaunutau2 ) p4toTLV(h0Wtaunutau2->p4(),tree.mc_truth_tth_h0Wtaunutau2_p4);
   if( h0Wq11 ) p4toTLV(h0Wq11->p4(),tree.mc_truth_tth_h0Wq11_p4);
   if( h0Wq21 ) p4toTLV(h0Wq21->p4(),tree.mc_truth_tth_h0Wq21_p4);
   if( h0Wq12 ) p4toTLV(h0Wq12->p4(),tree.mc_truth_tth_h0Wq12_p4);
   if( h0Wq22 ) p4toTLV(h0Wq22->p4(),tree.mc_truth_tth_h0Wq22_p4);
   
   if( h0Z1 ) p4toTLV(h0Z1->p4(),tree.mc_truth_tth_h0Z1_p4);
   if( h0Z2 ) p4toTLV(h0Z2->p4(),tree.mc_truth_tth_h0Z2_p4);
   if( h0Zl11 ) p4toTLV(h0Zl11->p4(),tree.mc_truth_tth_h0Zl11_p4);
   if( h0Zl21 ) p4toTLV(h0Zl21->p4(),tree.mc_truth_tth_h0Zl21_p4);
   if( h0Ztau11 ) p4toTLV(h0Ztau11->p4(),tree.mc_truth_tth_h0Ztau11_p4);
   if( h0Ztau21 ) p4toTLV(h0Ztau21->p4(),tree.mc_truth_tth_h0Ztau21_p4);
   if( h0Ztaul11 ) p4toTLV(h0Ztaul11->p4(),tree.mc_truth_tth_h0Ztaul11_p4);
   if( h0Ztaul21 ) p4toTLV(h0Ztaul21->p4(),tree.mc_truth_tth_h0Ztaul21_p4);
   if( h0Ztaunu11 ) p4toTLV(h0Ztaunu11->p4(),tree.mc_truth_tth_h0Ztaunu11_p4);
   if( h0Ztaunu21 ) p4toTLV(h0Ztaunu21->p4(),tree.mc_truth_tth_h0Ztaunu21_p4);
   if( h0Ztaunutau11 ) p4toTLV(h0Ztaunutau11->p4(),tree.mc_truth_tth_h0Ztaunutau11_p4);
   if( h0Ztaunutau21 ) p4toTLV(h0Ztaunutau21->p4(),tree.mc_truth_tth_h0Ztaunutau21_p4);
   if( h0Zq11 ) p4toTLV(h0Zq11->p4(),tree.mc_truth_tth_h0Zq11_p4);
   if( h0Zq21 ) p4toTLV(h0Zq21->p4(),tree.mc_truth_tth_h0Zq21_p4);
   if( h0Zl12 ) p4toTLV(h0Zl12->p4(),tree.mc_truth_tth_h0Zl12_p4);
   if( h0Zl22 ) p4toTLV(h0Zl22->p4(),tree.mc_truth_tth_h0Zl22_p4);
   if( h0Ztau12 ) p4toTLV(h0Ztau12->p4(),tree.mc_truth_tth_h0Ztau12_p4);
   if( h0Ztau22 ) p4toTLV(h0Ztau22->p4(),tree.mc_truth_tth_h0Ztau22_p4);
   if( h0Ztaul12 ) p4toTLV(h0Ztaul12->p4(),tree.mc_truth_tth_h0Ztaul12_p4);
   if( h0Ztaul22 ) p4toTLV(h0Ztaul22->p4(),tree.mc_truth_tth_h0Ztaul22_p4);
   if( h0Ztaunu12 ) p4toTLV(h0Ztaunu12->p4(),tree.mc_truth_tth_h0Ztaunu12_p4);
   if( h0Ztaunu22 ) p4toTLV(h0Ztaunu22->p4(),tree.mc_truth_tth_h0Ztaunu22_p4);
   if( h0Ztaunutau12 ) p4toTLV(h0Ztaunutau12->p4(),tree.mc_truth_tth_h0Ztaunutau12_p4);
   if( h0Ztaunutau22 ) p4toTLV(h0Ztaunutau22->p4(),tree.mc_truth_tth_h0Ztaunutau22_p4);
   if( h0Zq12 ) p4toTLV(h0Zq12->p4(),tree.mc_truth_tth_h0Zq12_p4);
   if( h0Zq22 ) p4toTLV(h0Zq22->p4(),tree.mc_truth_tth_h0Zq22_p4);
   if( h0Znu11 ) p4toTLV(h0Znu11->p4(),tree.mc_truth_tth_h0Znu11_p4);
   if( h0Znu21 ) p4toTLV(h0Znu21->p4(),tree.mc_truth_tth_h0Znu21_p4);
   if( h0Znu12 ) p4toTLV(h0Znu12->p4(),tree.mc_truth_tth_h0Znu12_p4);
   if( h0Znu22 ) p4toTLV(h0Znu22->p4(),tree.mc_truth_tth_h0Znu22_p4);
   
   if( h0tau1 ) p4toTLV(h0tau1->p4(),tree.mc_truth_tth_h0tau1_p4);
   if( h0tau2 ) p4toTLV(h0tau2->p4(),tree.mc_truth_tth_h0tau2_p4);
   if( h0taul1 ) p4toTLV(h0taul1->p4(),tree.mc_truth_tth_h0taul1_p4);
   if( h0taunutau1 ) p4toTLV(h0taunutau1->p4(),tree.mc_truth_tth_h0taunutau1_p4);
   if( h0taunu1 ) p4toTLV(h0taunu1->p4(),tree.mc_truth_tth_h0taunu1_p4);
   if( h0taul2 ) p4toTLV(h0taul2->p4(),tree.mc_truth_tth_h0taul2_p4);
   if( h0taunutau2 ) p4toTLV(h0taunutau2->p4(),tree.mc_truth_tth_h0taunutau2_p4);
   if( h0taunu2 ) p4toTLV(h0taunu2->p4(),tree.mc_truth_tth_h0taunu2_p4);
   
   if( t1 ) p4toTLV(t1->p4(),tree.mc_truth_tth_t1_p4);
   if( t2 ) p4toTLV(t2->p4(),tree.mc_truth_tth_t2_p4);
   if( tb1 ) p4toTLV(tb1->p4(),tree.mc_truth_tth_tb1_p4);
   if( tb2 ) p4toTLV(tb2->p4(),tree.mc_truth_tth_tb2_p4);
   
   if( tW1 ) p4toTLV(tW1->p4(),tree.mc_truth_tth_tW1_p4);
   if( tWnu1 ) p4toTLV(tWnu1->p4(),tree.mc_truth_tth_tWnu1_p4);
   if( tWnutau1 ) p4toTLV(tWnutau1->p4(),tree.mc_truth_tth_tWnutau1_p4);
   if( tWl1 ) p4toTLV(tWl1->p4(),tree.mc_truth_tth_tWl1_p4);
   if( tWtau1 ) p4toTLV(tWtau1->p4(),tree.mc_truth_tth_tWtau1_p4);
   if( tWtaunu1 ) p4toTLV(tWtaunu1->p4(),tree.mc_truth_tth_tWtaunu1_p4);
   if( tWtaunutau1 ) p4toTLV(tWtaunutau1->p4(),tree.mc_truth_tth_tWtaunutau1_p4);
   if( tWtaul1 ) p4toTLV(tWtaul1->p4(),tree.mc_truth_tth_tWtaul1_p4);
   if( tWq11 ) p4toTLV(tWq11->p4(),tree.mc_truth_tth_tWq11_p4);
   if( tWq21 ) p4toTLV(tWq21->p4(),tree.mc_truth_tth_tWq21_p4);

   if( tW2 ) p4toTLV(tW2->p4(),tree.mc_truth_tth_tW2_p4);
   if( tWnu2 ) p4toTLV(tWnu2->p4(),tree.mc_truth_tth_tWnu2_p4);
   if( tWnutau2 ) p4toTLV(tWnutau2->p4(),tree.mc_truth_tth_tWnutau2_p4);
   if( tWl2 ) p4toTLV(tWl2->p4(),tree.mc_truth_tth_tWl2_p4);
   if( tWtau2 ) p4toTLV(tWtau2->p4(),tree.mc_truth_tth_tWtau2_p4);
   if( tWtaunu2 ) p4toTLV(tWtaunu2->p4(),tree.mc_truth_tth_tWtaunu2_p4);
   if( tWtaunutau2 ) p4toTLV(tWtaunutau2->p4(),tree.mc_truth_tth_tWtaunutau2_p4);
   if( tWtaul2 ) p4toTLV(tWtaul2->p4(),tree.mc_truth_tth_tWtaul2_p4);
   if( tWq12 ) p4toTLV(tWq12->p4(),tree.mc_truth_tth_tWq12_p4);
   if( tWq22 ) p4toTLV(tWq22->p4(),tree.mc_truth_tth_tWq22_p4);

   if( j1 ) p4toTLV(j1->p4(),tree.mc_truth_tth_j1_p4);
   if( j2 ) p4toTLV(j2->p4(),tree.mc_truth_tth_j2_p4);
   if( j3 ) p4toTLV(j3->p4(),tree.mc_truth_tth_j3_p4);
   
   // pdgId
   
   if( h0 ) tree.mc_truth_tth_h0_id = h0->pdgId();

   if( h0W1 ) tree.mc_truth_tth_h0W1_id = h0W1->pdgId();
   if( h0W2 ) tree.mc_truth_tth_h0W2_id = h0W2->pdgId();
   if( h0Wl1 ) tree.mc_truth_tth_h0Wl1_id = h0Wl1->pdgId();
   if( h0Wnu1 ) tree.mc_truth_tth_h0Wnu1_id = h0Wnu1->pdgId();
   if( h0Wtau1 ) tree.mc_truth_tth_h0Wtau1_id = h0Wtau1->pdgId();
   if( h0Wnutau1 ) tree.mc_truth_tth_h0Wnutau1_id = h0Wnutau1->pdgId();
   if( h0Wtaul1 ) tree.mc_truth_tth_h0Wtaul1_id = h0Wtaul1->pdgId();
   if( h0Wtaunu1 ) tree.mc_truth_tth_h0Wtaunu1_id = h0Wtaunu1->pdgId();
   if( h0Wtaunutau1 ) tree.mc_truth_tth_h0Wtaunutau1_id = h0Wtaunutau1->pdgId();
   if( h0Wl2 ) tree.mc_truth_tth_h0Wl2_id = h0Wl2->pdgId();
   if( h0Wnu2 ) tree.mc_truth_tth_h0Wnu2_id = h0Wnu2->pdgId();
   if( h0Wtau2 ) tree.mc_truth_tth_h0Wtau2_id = h0Wtau2->pdgId();
   if( h0Wnutau2 ) tree.mc_truth_tth_h0Wnutau2_id = h0Wnutau2->pdgId();
   if( h0Wtaul2 ) tree.mc_truth_tth_h0Wtaul2_id = h0Wtaul2->pdgId();
   if( h0Wtaunu2 ) tree.mc_truth_tth_h0Wtaunu2_id = h0Wtaunu2->pdgId();
   if( h0Wtaunutau2 ) tree.mc_truth_tth_h0Wtaunutau2_id = h0Wtaunutau2->pdgId();
   if( h0Wq11 ) tree.mc_truth_tth_h0Wq11_id = h0Wq11->pdgId();
   if( h0Wq21 ) tree.mc_truth_tth_h0Wq21_id = h0Wq21->pdgId();
   if( h0Wq12 ) tree.mc_truth_tth_h0Wq12_id = h0Wq12->pdgId();
   if( h0Wq22 ) tree.mc_truth_tth_h0Wq22_id = h0Wq22->pdgId();
   
   if( h0Z1 ) tree.mc_truth_tth_h0Z1_id = h0Z1->pdgId();
   if( h0Z2 ) tree.mc_truth_tth_h0Z2_id = h0Z2->pdgId();
   if( h0Zl11 ) tree.mc_truth_tth_h0Zl11_id = h0Zl11->pdgId();
   if( h0Zl21 ) tree.mc_truth_tth_h0Zl21_id = h0Zl21->pdgId();
   if( h0Ztau11 ) tree.mc_truth_tth_h0Ztau11_id = h0Ztau11->pdgId();
   if( h0Ztau21 ) tree.mc_truth_tth_h0Ztau21_id = h0Ztau21->pdgId();
   if( h0Ztaul11 ) tree.mc_truth_tth_h0Ztaul11_id = h0Ztaul11->pdgId();
   if( h0Ztaul21 ) tree.mc_truth_tth_h0Ztaul21_id = h0Ztaul21->pdgId();
   if( h0Ztaunu11 ) tree.mc_truth_tth_h0Ztaunu11_id = h0Ztaunu11->pdgId();
   if( h0Ztaunu21 ) tree.mc_truth_tth_h0Ztaunu21_id = h0Ztaunu21->pdgId();
   if( h0Ztaunutau11 ) tree.mc_truth_tth_h0Ztaunutau11_id = h0Ztaunutau11->pdgId();
   if( h0Ztaunutau21 ) tree.mc_truth_tth_h0Ztaunutau21_id = h0Ztaunutau21->pdgId();
   if( h0Zq11 ) tree.mc_truth_tth_h0Zq11_id = h0Zq11->pdgId();
   if( h0Zq21 ) tree.mc_truth_tth_h0Zq21_id = h0Zq21->pdgId();
   if( h0Zq12 ) tree.mc_truth_tth_h0Zl12_id = h0Zl12->pdgId();
   if( h0Zq22 ) tree.mc_truth_tth_h0Zl22_id = h0Zl22->pdgId();
   if( h0Ztau12 ) tree.mc_truth_tth_h0Ztau12_id = h0Ztau12->pdgId();
   if( h0Ztau22 ) tree.mc_truth_tth_h0Ztau22_id = h0Ztau22->pdgId();
   if( h0Ztaul12 ) tree.mc_truth_tth_h0Ztaul12_id = h0Ztaul12->pdgId();
   if( h0Ztaul22 ) tree.mc_truth_tth_h0Ztaul22_id = h0Ztaul22->pdgId();
   if( h0Ztaunu12 ) tree.mc_truth_tth_h0Ztaunu12_id = h0Ztaunu12->pdgId();
   if( h0Ztaunu22 ) tree.mc_truth_tth_h0Ztaunu22_id = h0Ztaunu22->pdgId();
   if( h0Ztaunutau12 ) tree.mc_truth_tth_h0Ztaunutau12_id = h0Ztaunutau12->pdgId();
   if( h0Ztaunutau22 ) tree.mc_truth_tth_h0Ztaunutau22_id = h0Ztaunutau22->pdgId();
   if( h0Zq12 ) tree.mc_truth_tth_h0Zq12_id = h0Zq12->pdgId();
   if( h0Zq22 ) tree.mc_truth_tth_h0Zq22_id = h0Zq22->pdgId();
   if( h0Znu11 ) tree.mc_truth_tth_h0Znu11_id = h0Znu11->pdgId();
   if( h0Znu21 ) tree.mc_truth_tth_h0Znu21_id = h0Znu21->pdgId();
   if( h0Znu12 ) tree.mc_truth_tth_h0Znu12_id = h0Znu12->pdgId();
   if( h0Znu22 ) tree.mc_truth_tth_h0Znu22_id = h0Znu22->pdgId();
   
   if( h0tau1 ) tree.mc_truth_tth_h0tau1_id = h0tau1->pdgId();
   if( h0tau2 ) tree.mc_truth_tth_h0tau2_id = h0tau2->pdgId();
   if( h0taul1 ) tree.mc_truth_tth_h0taul1_id = h0taul1->pdgId();
   if( h0taunutau1 ) tree.mc_truth_tth_h0taunutau1_id = h0taunutau1->pdgId();
   if( h0taunu1 ) tree.mc_truth_tth_h0taunu1_id = h0taunu1->pdgId();
   if( h0taul2 ) tree.mc_truth_tth_h0taul2_id = h0taul2->pdgId();
   if( h0taunutau2 ) tree.mc_truth_tth_h0taunutau2_id = h0taunutau2->pdgId();
   if( h0taunu2 ) tree.mc_truth_tth_h0taunu2_id = h0taunu2->pdgId();
   
   if( t1 ) tree.mc_truth_tth_t1_id = t1->pdgId();
   if( t2 ) tree.mc_truth_tth_t2_id = t2->pdgId();
   if( tb1 ) tree.mc_truth_tth_tb1_id = tb1->pdgId();
   if( tb2 ) tree.mc_truth_tth_tb2_id = tb2->pdgId();
   
   if( tW1 ) tree.mc_truth_tth_tW1_id = tW1->pdgId();
   if( tWnu1 ) tree.mc_truth_tth_tWnu1_id = tWnu1->pdgId();
   if( tWnutau1 ) tree.mc_truth_tth_tWnutau1_id = tWnutau1->pdgId();
   if( tWl1 ) tree.mc_truth_tth_tWl1_id = tWl1->pdgId();
   if( tWtau1 ) tree.mc_truth_tth_tWtau1_id = tWtau1->pdgId();
   if( tWtaunu1 ) tree.mc_truth_tth_tWtaunu1_id = tWtaunu1->pdgId();
   if( tWtaunutau1 ) tree.mc_truth_tth_tWtaunutau1_id = tWtaunutau1->pdgId();
   if( tWtaul1 ) tree.mc_truth_tth_tWtaul1_id = tWtaul1->pdgId();
   if( tWq11 ) tree.mc_truth_tth_tWq11_id = tWq11->pdgId();
   if( tWq21 ) tree.mc_truth_tth_tWq21_id = tWq21->pdgId();

   if( tW2 ) tree.mc_truth_tth_tW2_id = tW2->pdgId();
   if( tWnu2 ) tree.mc_truth_tth_tWnu2_id = tWnu2->pdgId();
   if( tWnutau2 ) tree.mc_truth_tth_tWnutau2_id = tWnutau2->pdgId();
   if( tWl2 ) tree.mc_truth_tth_tWl2_id = tWl2->pdgId();
   if( tWtau2 ) tree.mc_truth_tth_tWtau2_id = tWtau2->pdgId();
   if( tWtaunu2 ) tree.mc_truth_tth_tWtaunu2_id = tWtaunu2->pdgId();
   if( tWtaunutau2 ) tree.mc_truth_tth_tWtaunutau2_id = tWtaunutau2->pdgId();
   if( tWtaul2 ) tree.mc_truth_tth_tWtaul2_id = tWtaul2->pdgId();
   if( tWq12 ) tree.mc_truth_tth_tWq12_id = tWq12->pdgId();
   if( tWq22 ) tree.mc_truth_tth_tWq22_id = tWq22->pdgId();

   if( j1 ) tree.mc_truth_tth_j1_id = j1->pdgId();
   if( j2 ) tree.mc_truth_tth_j2_id = j2->pdgId();
   if( j3 ) tree.mc_truth_tth_j3_id = j3->pdgId();

   // status

   if( h0 ) tree.mc_truth_tth_h0_status = h0->status();

   if( h0W1 ) tree.mc_truth_tth_h0W1_status = h0W1->status();
   if( h0W2 ) tree.mc_truth_tth_h0W2_status = h0W2->status();
   if( h0Wl1 ) tree.mc_truth_tth_h0Wl1_status = h0Wl1->status();
   if( h0Wnu1 ) tree.mc_truth_tth_h0Wnu1_status = h0Wnu1->status();
   if( h0Wtau1 ) tree.mc_truth_tth_h0Wtau1_status = h0Wtau1->status();
   if( h0Wnutau1 ) tree.mc_truth_tth_h0Wnutau1_status = h0Wnutau1->status();
   if( h0Wtaul1 ) tree.mc_truth_tth_h0Wtaul1_status = h0Wtaul1->status();
   if( h0Wtaunu1 ) tree.mc_truth_tth_h0Wtaunu1_status = h0Wtaunu1->status();
   if( h0Wtaunutau1 ) tree.mc_truth_tth_h0Wtaunutau1_status = h0Wtaunutau1->status();
   if( h0Wl2 ) tree.mc_truth_tth_h0Wl2_status = h0Wl2->status();
   if( h0Wnu2 ) tree.mc_truth_tth_h0Wnu2_status = h0Wnu2->status();
   if( h0Wtau2 ) tree.mc_truth_tth_h0Wtau2_status = h0Wtau2->status();
   if( h0Wnutau2 ) tree.mc_truth_tth_h0Wnutau2_status = h0Wnutau2->status();
   if( h0Wtaul2 ) tree.mc_truth_tth_h0Wtaul2_status = h0Wtaul2->status();
   if( h0Wtaunu2 ) tree.mc_truth_tth_h0Wtaunu2_status = h0Wtaunu2->status();
   if( h0Wtaunutau2 ) tree.mc_truth_tth_h0Wtaunutau2_status = h0Wtaunutau2->status();
   if( h0Wq11 ) tree.mc_truth_tth_h0Wq11_status = h0Wq11->status();
   if( h0Wq21 ) tree.mc_truth_tth_h0Wq21_status = h0Wq21->status();
   if( h0Wq12 ) tree.mc_truth_tth_h0Wq12_status = h0Wq12->status();
   if( h0Wq22 ) tree.mc_truth_tth_h0Wq22_status = h0Wq22->status();
   
   if( h0Z1 ) tree.mc_truth_tth_h0Z1_status = h0Z1->status();
   if( h0Z2 ) tree.mc_truth_tth_h0Z2_status = h0Z2->status();
   if( h0Zl11 ) tree.mc_truth_tth_h0Zl11_status = h0Zl11->status();
   if( h0Zl21 ) tree.mc_truth_tth_h0Zl21_status = h0Zl21->status();
   if( h0Ztau11 ) tree.mc_truth_tth_h0Ztau11_status = h0Ztau11->status();
   if( h0Ztau21 ) tree.mc_truth_tth_h0Ztau21_status = h0Ztau21->status();
   if( h0Ztaul11 ) tree.mc_truth_tth_h0Ztaul11_status = h0Ztaul11->status();
   if( h0Ztaul21 ) tree.mc_truth_tth_h0Ztaul21_status = h0Ztaul21->status();
   if( h0Ztaunu11 ) tree.mc_truth_tth_h0Ztaunu11_status = h0Ztaunu11->status();
   if( h0Ztaunu21 ) tree.mc_truth_tth_h0Ztaunu21_status = h0Ztaunu21->status();
   if( h0Ztaunutau11 ) tree.mc_truth_tth_h0Ztaunutau11_status = h0Ztaunutau11->status();
   if( h0Ztaunutau21 ) tree.mc_truth_tth_h0Ztaunutau21_status = h0Ztaunutau21->status();
   if( h0Zq11 ) tree.mc_truth_tth_h0Zq11_status = h0Zq11->status();
   if( h0Zq21 ) tree.mc_truth_tth_h0Zq21_status = h0Zq21->status();
   if( h0Zl12 ) tree.mc_truth_tth_h0Zl12_status = h0Zl12->status();
   if( h0Zl22 ) tree.mc_truth_tth_h0Zl22_status = h0Zl22->status();
   if( h0Ztau12 ) tree.mc_truth_tth_h0Ztau12_status = h0Ztau12->status();
   if( h0Ztau22 ) tree.mc_truth_tth_h0Ztau22_status = h0Ztau22->status();
   if( h0Ztaul12 ) tree.mc_truth_tth_h0Ztaul12_status = h0Ztaul12->status();
   if( h0Ztaul22 ) tree.mc_truth_tth_h0Ztaul22_status = h0Ztaul22->status();
   if( h0Ztaunu12 ) tree.mc_truth_tth_h0Ztaunu12_status = h0Ztaunu12->status();
   if( h0Ztaunu22 ) tree.mc_truth_tth_h0Ztaunu22_status = h0Ztaunu22->status();
   if( h0Ztaunutau12 ) tree.mc_truth_tth_h0Ztaunutau12_status = h0Ztaunutau12->status();
   if( h0Ztaunutau22 ) tree.mc_truth_tth_h0Ztaunutau22_status = h0Ztaunutau22->status();
   if( h0Zq12 ) tree.mc_truth_tth_h0Zq12_status = h0Zq12->status();
   if( h0Zq22 ) tree.mc_truth_tth_h0Zq22_status = h0Zq22->status();
   if( h0Znu11 ) tree.mc_truth_tth_h0Znu11_status = h0Znu11->status();
   if( h0Znu21 ) tree.mc_truth_tth_h0Znu21_status = h0Znu21->status();
   if( h0Znu12 ) tree.mc_truth_tth_h0Znu12_status = h0Znu12->status();
   if( h0Znu22 ) tree.mc_truth_tth_h0Znu22_status = h0Znu22->status();
   
   if( h0tau1 ) tree.mc_truth_tth_h0tau1_status = h0tau1->status();
   if( h0tau2 ) tree.mc_truth_tth_h0tau2_status = h0tau2->status();
   if( h0taul1 ) tree.mc_truth_tth_h0taul1_status = h0taul1->status();
   if( h0taunutau1 ) tree.mc_truth_tth_h0taunutau1_status = h0taunutau1->status();
   if( h0taunu1 ) tree.mc_truth_tth_h0taunu1_status = h0taunu1->status();
   if( h0taul2 ) tree.mc_truth_tth_h0taul2_status = h0taul2->status();
   if( h0taunutau2 ) tree.mc_truth_tth_h0taunutau2_status = h0taunutau2->status();
   if( h0taunu2 ) tree.mc_truth_tth_h0taunu2_status = h0taunu2->status();
   
   if( t1 ) tree.mc_truth_tth_t1_status = t1->status();
   if( t2 ) tree.mc_truth_tth_t2_status = t2->status();
   if( tb1 ) tree.mc_truth_tth_tb1_status = tb1->status();
   if( tb2 ) tree.mc_truth_tth_tb2_status = tb2->status();
   
   if( tW1 ) tree.mc_truth_tth_tW1_status = tW1->status();
   if( tWnu1 ) tree.mc_truth_tth_tWnu1_status = tWnu1->status();
   if( tWnutau1 ) tree.mc_truth_tth_tWnutau1_status = tWnutau1->status();
   if( tWl1 ) tree.mc_truth_tth_tWl1_status = tWl1->status();
   if( tWtau1 ) tree.mc_truth_tth_tWtau1_status = tWtau1->status();
   if( tWtaunu1 ) tree.mc_truth_tth_tWtaunu1_status = tWtaunu1->status();
   if( tWtaunutau1 ) tree.mc_truth_tth_tWtaunutau1_status = tWtaunutau1->status();
   if( tWtaul1 ) tree.mc_truth_tth_tWtaul1_status = tWtaul1->status();
   if( tWq11 ) tree.mc_truth_tth_tWq11_status = tWq11->status();
   if( tWq21 ) tree.mc_truth_tth_tWq21_status = tWq21->status();

   if( tW2 ) tree.mc_truth_tth_tW2_status = tW2->status();
   if( tWnu2 ) tree.mc_truth_tth_tWnu2_status = tWnu2->status();
   if( tWnutau2 ) tree.mc_truth_tth_tWnutau2_status = tWnutau2->status();
   if( tWl2 ) tree.mc_truth_tth_tWl2_status = tWl2->status();
   if( tWtau2 ) tree.mc_truth_tth_tWtau2_status = tWtau2->status();
   if( tWtaunu2 ) tree.mc_truth_tth_tWtaunu2_status = tWtaunu2->status();
   if( tWtaunutau2 ) tree.mc_truth_tth_tWtaunutau2_status = tWtaunutau2->status();
   if( tWtaul2 ) tree.mc_truth_tth_tWtaul2_status = tWtaul2->status();
   if( tWq12 ) tree.mc_truth_tth_tWq12_status = tWq12->status();
   if( tWq22 ) tree.mc_truth_tth_tWq22_status = tWq22->status();

   if( j1 ) tree.mc_truth_tth_j1_status = j1->status();
   if( j2 ) tree.mc_truth_tth_j2_status = j2->status();
   if( j3 ) tree.mc_truth_tth_j3_status = j3->status();
}

// tZq MC analyzer
void MCTruth::fillTZQSignalGenParticles(const edm::Event& iEvent,
					const edm::EventSetup& iSetup,
					FlatTree& tree,
					const edm::Handle<std::vector<reco::GenParticle> >& GenParticles)
{
   reco::GenParticle *Z = 0;
   
   reco::GenParticle *Zl1 = 0;
   reco::GenParticle *Zl2 = 0;
   reco::GenParticle *Ztau1 = 0;
   reco::GenParticle *Ztau2 = 0;
   reco::GenParticle *Ztaul1 = 0;
   reco::GenParticle *Ztaul2 = 0;
   reco::GenParticle *Ztaunu1 = 0;
   reco::GenParticle *Ztaunu2 = 0;
   reco::GenParticle *Ztaunutau1 = 0;
   reco::GenParticle *Ztaunutau2 = 0;
   
   reco::GenParticle *t = 0;
   reco::GenParticle *tb = 0;
   reco::GenParticle *tW = 0;
   reco::GenParticle *tWnu = 0;
   reco::GenParticle *tWnutau = 0;
   reco::GenParticle *tWl = 0;
   reco::GenParticle *tWtau = 0;
   reco::GenParticle *tWtaunu = 0;
   reco::GenParticle *tWtaunutau = 0;
   reco::GenParticle *tWtaul = 0;
   reco::GenParticle *tWq1 = 0;
   reco::GenParticle *tWq2 = 0;

   reco::GenParticle *j1 = 0;
   reco::GenParticle *j2 = 0;
   reco::GenParticle *j3 = 0;
   
   int chan = -666;

   // 0   = (t->bW,W->lnu)
   // 1   = (t->bW,W->qq)
   // 2   = (t->bW,W->tauLnu)
   // 3   = (t->bW,W->tauHnu)
   
   // (Z->ll)             +0
   // (Z->tauLtauL)       +20
   // (Z->tauHtauH)       +40
   // (Z->tauLtauH)       +60
   // (Z->tauHtauL)       +80
         
   reco::GenParticleCollection genParticlesCollection = *GenParticles;
   reco::GenParticleCollection::const_iterator genParticleSrc;
   
   int ipart = 0;

//   std::cout << "event" << std::endl;
   
   for(genParticleSrc = genParticlesCollection.begin();
       genParticleSrc != genParticlesCollection.end(); 
       genParticleSrc++)
     {
	reco::GenParticle *mcp = &(const_cast<reco::GenParticle&>(*genParticleSrc));

	int barcode = ipart; // in CMSSW barcode is the index of genParticle in the event
	// https://twiki.cern.ch/twiki/bin/view/CMS/GenParticles2HepMCConverter
	ipart++;
	
//	if( fabs(mcp->pdgId()) <= 6 )
//	  {	     
//	     std::cout << "pdgId=" << mcp->pdgId() << " pt=" <<
//	       mcp->pt() << " status=" << mcp->status() << 
//	       " motherPdg=" << mcp->mother()->pdgId() << std::endl;
//	  }	
	
	// Additional partons (up to three)
	if( (fabs(mcp->pdgId()) <= 6 || fabs(mcp->pdgId()) == 21) &&
	    mcp->status() == 23 && barcode == 8 )
	  {
	     if( !j1 )
	       j1 = mcp;
	     else if( !j2 )
	       j2 = mcp;
	     else if( !j3 )
	       j3 = mcp;
	  }	

//	if( fabs(mcp->pdgId()) == 15 )
//	  {	     
//	     std::cout << mcp->pdgId() << " " << mcp->status() << " " << mcp->numberOfDaughters() << std::endl;
//	     
//	     const reco::GenParticleRefVector& daughterRefs = mcp->daughterRefVector();
//	     for(reco::GenParticleRefVector::const_iterator mcp_idr = daughterRefs.begin(); 
//		 mcp_idr!= daughterRefs.end(); ++mcp_idr)
//	       {
//		  if( mcp_idr->isAvailable() ) 
//		    {
//		       const reco::GenParticleRef& genParticle = (*mcp_idr);
//		       const reco::GenParticle *mcp_d = genParticle.get();
//		       std::cout << "daughter " << mcp_d->pdgId() << " " << mcp_d->status() << std::endl;
//		    }
//	       }	     	     
//	  }

	if( ((fabs(mcp->pdgId()) == 11 ||
	      fabs(mcp->pdgId()) == 13) &&
	     mcp->status() == 3) ||
	    (fabs(mcp->pdgId()) == 15 &&
		mcp->status() == 2) )
	  {
	     if( fabs(mcp->pdgId()) == 11 ||
		 fabs(mcp->pdgId()) == 13 ) // l
	       {
		  if( Zl1 && !Zl2 ) {Zl2 = mcp;}
		  if( !Zl1 ) {Zl1 = mcp;}
	       }			    
	     if( fabs(mcp->pdgId()) == 15 ) // tau
	       {
		  if( Ztau1 )
		    {
		       Ztau2 = mcp;
				      
		       const reco::GenParticleRefVector& daughterRefs = Ztau2->daughterRefVector();
		       for(reco::GenParticleRefVector::const_iterator Ztau2_idr = daughterRefs.begin(); 
			   Ztau2_idr!= daughterRefs.end(); ++Ztau2_idr)
			 {
			    if( Ztau2_idr->isAvailable() ) 
			      {		       
				 const reco::GenParticleRef& genParticle = (*Ztau2_idr);
				 const reco::GenParticle *Ztau2_d = genParticle.get();
				 reco::GenParticle *pfff = getUnique(Ztau2_d,0);
				 
				 if( fabs(pfff->pdgId()) == 12 ||
				     fabs(pfff->pdgId()) == 14 ) // nu
				   {
				      Ztaunu2 = pfff;
				   }
				 if( fabs(pfff->pdgId()) == 16 ) // nutau
				   {
				      Ztaunutau2 = pfff;
				   }
				 if( fabs(pfff->pdgId()) == 11 ||
				     fabs(pfff->pdgId()) == 13 ) // l
				   {
				      Ztaul2 = pfff;
				   }
			      }
			 }							  
		    }
		  if( !Ztau1 )
		    {
		       Ztau1 = mcp;
		       
		       const reco::GenParticleRefVector& daughterRefs = Ztau1->daughterRefVector();
		       for(reco::GenParticleRefVector::const_iterator Ztau1_idr = daughterRefs.begin(); 
			   Ztau1_idr!= daughterRefs.end(); ++Ztau1_idr)
			 {
			    if( Ztau1_idr->isAvailable() ) 
			      {		       
				 const reco::GenParticleRef& genParticle = (*Ztau1_idr);
				 const reco::GenParticle *Ztau1_d = genParticle.get();
				 reco::GenParticle *pfff = getUnique(Ztau1_d,0);
				 
				 if( fabs(pfff->pdgId()) == 12 ||
				     fabs(pfff->pdgId()) == 14 ) // nu
				   {
				      Ztaunu1 = pfff;
				   }
				 if( fabs(pfff->pdgId()) == 16 ) // nutau
				   {
				      Ztaunutau1 = pfff;
				   }
				 if( fabs(pfff->pdgId()) == 11 ||
				     fabs(pfff->pdgId()) == 13 ) // l
				   {
				      Ztaul1 = pfff;
				   }
			      }
			 }							  
		    }
	       }	     
	  }		  
	
	// Z
	if( fabs(mcp->pdgId()) == 23 )
	  {
	     if( !Z ) {Z = mcp;}
	     
	     if( Z )
	       {
		  const reco::GenParticleRefVector& daughterRefs = Z->daughterRefVector();
		  for(reco::GenParticleRefVector::const_iterator Z_idr = daughterRefs.begin(); 
		      Z_idr!= daughterRefs.end(); ++Z_idr)
		    {
		       if( Z_idr->isAvailable() ) 
			 {		       
			    const reco::GenParticleRef& genParticle = (*Z_idr);
			    const reco::GenParticle *Z_d = genParticle.get();
			    reco::GenParticle *pff = getUnique(Z_d,0);
			    
			    if( fabs(pff->pdgId()) == 11 ||
				fabs(pff->pdgId()) == 13 ) // l
			      {
				 if( Zl1 && !Zl2 ) {Zl2 = pff;}
				 if( !Zl1 ) {Zl1 = pff;}
			      }			    
			    if( fabs(pff->pdgId()) == 15 ) // tau
			      {
				 if( Ztau1 )
				   {
				      Ztau2 = pff;
				      
				      const reco::GenParticleRefVector& daughterRefs = Ztau2->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator Ztau2_idr = daughterRefs.begin(); 
					  Ztau2_idr!= daughterRefs.end(); ++Ztau2_idr)
					{
					   if( Ztau2_idr->isAvailable() ) 
					     {		       
						const reco::GenParticleRef& genParticle = (*Ztau2_idr);
						const reco::GenParticle *Ztau2_d = genParticle.get();
						reco::GenParticle *pfff = getUnique(Ztau2_d,0);
						
						if( fabs(pfff->pdgId()) == 12 ||
						    fabs(pfff->pdgId()) == 14 ) // nu
						  {
						     Ztaunu2 = pfff;
						  }
						if( fabs(pfff->pdgId()) == 16 ) // nutau
						  {
						     Ztaunutau2 = pfff;
						  }
						if( fabs(pfff->pdgId()) == 11 ||
						    fabs(pfff->pdgId()) == 13 ) // l
						  {
						     Ztaul2 = pfff;
						  }
					     }
					}							  
				   }
				 if( !Ztau1 )
				   {
				      Ztau1 = pff;
				      
				      const reco::GenParticleRefVector& daughterRefs = Ztau1->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator Ztau1_idr = daughterRefs.begin(); 
					  Ztau1_idr!= daughterRefs.end(); ++Ztau1_idr)
					{
					   if( Ztau1_idr->isAvailable() ) 
					     {		       
						const reco::GenParticleRef& genParticle = (*Ztau1_idr);
						const reco::GenParticle *Ztau1_d = genParticle.get();
						reco::GenParticle *pfff = getUnique(Ztau1_d,0);
						
						if( fabs(pfff->pdgId()) == 12 ||
						    fabs(pfff->pdgId()) == 14 ) // nu
						  {
						     Ztaunu1 = pfff;
						  }
						if( fabs(pfff->pdgId()) == 16 ) // nutau
						  {
						     Ztaunutau1 = pfff;
						  }
						if( fabs(pfff->pdgId()) == 11 ||
						    fabs(pfff->pdgId()) == 13 ) // l
						  {
						     Ztaul1 = pfff;
						  }
					     }
					}							  
				   }						     
			      }
			 }					   
		    }				      
	       }
	  }	
	
	// top decays
	if( fabs(mcp->pdgId()) == 6
	    && ( (mcp->status() == 62) || 
		 (mcp->status() == 3)
	       ) )
	  {
	     if( !t ) {t = const_cast<reco::GenParticle*>(mcp);}

	     const reco::GenParticleRefVector& daughterRefs = mcp->daughterRefVector();
	     for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr) 
	       {
		  if( idr->isAvailable() ) 
		    {		       
		       const reco::GenParticleRef& genParticle = (*idr);
		       const reco::GenParticle *d = genParticle.get();
		       reco::GenParticle *pf = getUnique(d,0);
		       
//		       if( pf->status() != 3 && pf->status() != 62 ) continue;
		       
		       if( fabs(pf->pdgId()) == 5 || fabs(pf->pdgId()) == 3 || fabs(pf->pdgId()) == 1 ) // b or s or d
			 {
			    if( !tb ) {tb = pf;}
			 }		       
		       
		       if( fabs(pf->pdgId()) == 24 ) // W
			 {
			    if( !tW )
			      {
				 tW = pf;
				 const reco::GenParticleRefVector& tW_daughterRefs = tW->daughterRefVector();
				 for(reco::GenParticleRefVector::const_iterator tW_idr = tW_daughterRefs.begin();
				     tW_idr!= tW_daughterRefs.end(); ++tW_idr)
				   {
				      if( tW_idr->isAvailable() ) 
					{		       
					   const reco::GenParticleRef& tW_genParticle = (*tW_idr);
					   const reco::GenParticle *tW_d = tW_genParticle.get();
					   reco::GenParticle *pff = getUnique(tW_d,0);
					   
					   if( fabs(pff->pdgId()) == 12 ||
					       fabs(pff->pdgId()) == 14 ) // nu
					     {
						tWnu = pff;
					     }		
					   if( fabs(pff->pdgId()) == 16 ) // nu_tau
					     {
						tWnutau = pff;
					     }		
					   if( fabs(pff->pdgId()) == 11 ||
					       fabs(pff->pdgId()) == 13 ) // l
					     {
						tWl = pff;
					     }		
					   if( fabs(pff->pdgId()) == 15 ) // tau
					     {
						tWtau = pff;
						
						const reco::GenParticleRefVector& tWtau_daughterRefs = tWtau->daughterRefVector();
						for(reco::GenParticleRefVector::const_iterator tWtau_idr = tWtau_daughterRefs.begin();
						    tWtau_idr!= tWtau_daughterRefs.end(); ++tWtau_idr)
						  {
						     if( tWtau_idr->isAvailable() ) 
						       {		       
							  const reco::GenParticleRef& tWtau_genParticle = (*tWtau_idr);
							  const reco::GenParticle *tWtau_d = tWtau_genParticle.get();
							  reco::GenParticle *pfff = getUnique(tWtau_d,0);
							  
							  if( fabs(pfff->pdgId()) == 12 ||
							      fabs(pfff->pdgId()) == 14 ) // nu
							    {
							       tWtaunu = pfff;
							    }		
							  if( fabs(pfff->pdgId()) == 16 ) // nu_tau
							    {
							       tWtaunutau = pfff;
							    }			
							  if( fabs(pfff->pdgId()) == 11 ||
							      fabs(pfff->pdgId()) == 13 ) // l
							    {
							       tWtaul = pfff;
							    }							  
						       }
						  }						
					     }
					   if( fabs(pff->pdgId()) <= 6 ) // q
					     {
						if( tWq1 && !tWq2 ) {tWq2 = pff;}
						if( !tWq1 ) {tWq1 = pff;}
					     }					   					   
					}
				   }				
			      }			    
			 }		       
		    }
	       }	     
	  }
     }

   bool doCheck = 0;
   
   if( t && tb && tW )
     {	
	int tchan = -666;
	if( tWl )   tchan = 0;
	if( tWq1 && tWq2 ) tchan = 1;
	if( tWtaul )  tchan = 2;
	if( tWtaunutau && !tWtaul )  tchan = 3;
	
	if( tchan < 0 && doCheck )
	  {	     
	     std::cout << "Failed to identify top-quark decay chain" << std::endl;
	     
	     std::cout << "t = " << bool(t) << std::endl;
	     std::cout << "t->W = " << bool(tW) << std::endl;
	     std::cout << "t->W->l = " << bool(tWl) << std::endl;
	     std::cout << "t->W->nu = " << bool(tWnu) << std::endl;
	     std::cout << "t->W->nutau = " << bool(tWnutau) << std::endl;
	     std::cout << "t->W->tau = " << bool(tWtau) << std::endl;
	     std::cout << "t->W->tau->l = " << bool(tWtaul) << std::endl;
	     std::cout << "t->W->tau->nu = " << bool(tWtaunu) << std::endl;
	     std::cout << "t->W->tau->nutau = " << bool(tWtaunutau) << std::endl;
	     std::cout << "t->W->q = " << bool(tWq1) << std::endl;

	     exit(1);
	  }
	
	if( Zl1 || Ztaul1 || Ztaunutau1 )
	  {
	     int chan0 = 0;
	     if( Zl1 ) chan = chan0 + 0 + tchan;
	     if( Ztaul1 && Ztaul2 ) chan = chan0 + 20 + tchan;
	     if( Ztaunutau1 && ! Ztaul1 && Ztaunutau2 && ! Ztaul2 ) chan = chan0 + 40 + tchan;
	     if( Ztaul1 && Ztaunutau2 && ! Ztaul2 ) chan = chan0 + 60 + tchan;
	     if( Ztaul2 && Ztaunutau1 && ! Ztaul1 ) chan = chan0 + 80 + tchan;
	  }	

	if( chan < 0 && doCheck )
	  {
	     std::cout << "Unknown channel found" << std::endl;
	     std::cout << "chan = " << chan << std::endl;

	     std::cout << "j1 = " << bool(j1) << std::endl;
	     std::cout << "j2 = " << bool(j2) << std::endl;
	     std::cout << "j3 = " << bool(j3) << std::endl;
	     
	     std::cout << "Z = " << bool(Z) << std::endl;
	     std::cout << "Z->l1 = " << bool(Zl1) << std::endl;
	     std::cout << "Z->l2 = " << bool(Zl2) << std::endl;
	     std::cout << "Z->tau1 = " << bool(Ztau1) << std::endl;
	     std::cout << "Z->tau1->l = " << bool(Ztaul1) << std::endl;
	     std::cout << "Z->tau1->nu = " << bool(Ztaunu1) << std::endl;
	     std::cout << "Z->tau1->nutau = " << bool(Ztaunutau1) << std::endl;
	     std::cout << "Z->tau2 = " << bool(Ztau2) << std::endl;
	     std::cout << "Z->tau2->l = " << bool(Ztaul2) << std::endl;
	     std::cout << "Z->tau2->nu = " << bool(Ztaunu2) << std::endl;
	     std::cout << "Z->tau2->nutau = " << bool(Ztaunutau2) << std::endl;
	     
	     exit(1);
	  }
     }

   tree.mc_truth_tzq_channel = chan;

   // TLV

   if( Z ) p4toTLV(Z->p4(),tree.mc_truth_tzq_Z_p4);
   if( Zl1 ) p4toTLV(Zl1->p4(),tree.mc_truth_tzq_Zl1_p4);
   if( Zl2 ) p4toTLV(Zl2->p4(),tree.mc_truth_tzq_Zl2_p4);
   if( Ztau1 ) p4toTLV(Ztau1->p4(),tree.mc_truth_tzq_Ztau1_p4);
   if( Ztau2 ) p4toTLV(Ztau2->p4(),tree.mc_truth_tzq_Ztau2_p4);
   if( Ztaul1 ) p4toTLV(Ztaul1->p4(),tree.mc_truth_tzq_Ztaul1_p4);
   if( Ztaul2 ) p4toTLV(Ztaul2->p4(),tree.mc_truth_tzq_Ztaul2_p4);
   if( Ztaunu1 ) p4toTLV(Ztaunu1->p4(),tree.mc_truth_tzq_Ztaunu1_p4);
   if( Ztaunu2 ) p4toTLV(Ztaunu2->p4(),tree.mc_truth_tzq_Ztaunu2_p4);
   if( Ztaunutau1 ) p4toTLV(Ztaunutau1->p4(),tree.mc_truth_tzq_Ztaunutau1_p4);
   if( Ztaunutau2 ) p4toTLV(Ztaunutau2->p4(),tree.mc_truth_tzq_Ztaunutau2_p4);
   
   if( t ) p4toTLV(t->p4(),tree.mc_truth_tzq_t_p4);
   if( tb ) p4toTLV(tb->p4(),tree.mc_truth_tzq_tb_p4);
   if( tW ) p4toTLV(tW->p4(),tree.mc_truth_tzq_tW_p4);
   if( tWnu ) p4toTLV(tWnu->p4(),tree.mc_truth_tzq_tWnu_p4);
   if( tWnutau ) p4toTLV(tWnutau->p4(),tree.mc_truth_tzq_tWnutau_p4);
   if( tWl ) p4toTLV(tWl->p4(),tree.mc_truth_tzq_tWl_p4);
   if( tWtau ) p4toTLV(tWtau->p4(),tree.mc_truth_tzq_tWtau_p4);
   if( tWtaunu ) p4toTLV(tWtaunu->p4(),tree.mc_truth_tzq_tWtaunu_p4);
   if( tWtaunutau ) p4toTLV(tWtaunutau->p4(),tree.mc_truth_tzq_tWtaunutau_p4);
   if( tWtaul ) p4toTLV(tWtaul->p4(),tree.mc_truth_tzq_tWtaul_p4);
   if( tWq1 ) p4toTLV(tWq1->p4(),tree.mc_truth_tzq_tWq1_p4);
   if( tWq2 ) p4toTLV(tWq2->p4(),tree.mc_truth_tzq_tWq2_p4);

   if( j1 ) p4toTLV(j1->p4(),tree.mc_truth_tzq_j1_p4);
   if( j2 ) p4toTLV(j2->p4(),tree.mc_truth_tzq_j2_p4);
   if( j3 ) p4toTLV(j3->p4(),tree.mc_truth_tzq_j3_p4);

   // pdgId

   if( Z ) tree.mc_truth_tzq_Z_id = Z->pdgId();
   if( Zl1 ) tree.mc_truth_tzq_Zl1_id = Zl1->pdgId();
   if( Zl2 ) tree.mc_truth_tzq_Zl2_id = Zl2->pdgId();
   if( Ztau1 ) tree.mc_truth_tzq_Ztau1_id = Ztau1->pdgId();
   if( Ztau2 ) tree.mc_truth_tzq_Ztau2_id = Ztau2->pdgId();
   if( Ztaul1 ) tree.mc_truth_tzq_Ztaul1_id = Ztaul1->pdgId();
   if( Ztaul2 ) tree.mc_truth_tzq_Ztaul2_id = Ztaul2->pdgId();
   if( Ztaunu1 ) tree.mc_truth_tzq_Ztaunu1_id = Ztaunu1->pdgId();
   if( Ztaunu2 ) tree.mc_truth_tzq_Ztaunu2_id = Ztaunu2->pdgId();
   if( Ztaunutau1 ) tree.mc_truth_tzq_Ztaunutau1_id = Ztaunutau1->pdgId();
   if( Ztaunutau2 ) tree.mc_truth_tzq_Ztaunutau2_id = Ztaunutau2->pdgId();
   
   if( t ) tree.mc_truth_tzq_t_id = t->pdgId();
   if( tb ) tree.mc_truth_tzq_tb_id = tb->pdgId();
   if( tW ) tree.mc_truth_tzq_tW_id = tW->pdgId();
   if( tWnu ) tree.mc_truth_tzq_tWnu_id = tWnu->pdgId();
   if( tWnutau ) tree.mc_truth_tzq_tWnutau_id = tWnutau->pdgId();
   if( tWl ) tree.mc_truth_tzq_tWl_id = tWl->pdgId();
   if( tWtau ) tree.mc_truth_tzq_tWtau_id = tWtau->pdgId();
   if( tWtaunu ) tree.mc_truth_tzq_tWtaunu_id = tWtaunu->pdgId();
   if( tWtaunutau ) tree.mc_truth_tzq_tWtaunutau_id = tWtaunutau->pdgId();
   if( tWtaul ) tree.mc_truth_tzq_tWtaul_id = tWtaul->pdgId();
   if( tWq1 ) tree.mc_truth_tzq_tWq1_id = tWq1->pdgId();
   if( tWq2 ) tree.mc_truth_tzq_tWq2_id = tWq2->pdgId();
   
   if( j1 ) tree.mc_truth_tzq_j1_id = j1->pdgId();
   if( j2 ) tree.mc_truth_tzq_j2_id = j2->pdgId();
   if( j3 ) tree.mc_truth_tzq_j3_id = j3->pdgId();

   // status
   
   if( Z ) tree.mc_truth_tzq_Z_status = Z->status();
   if( Zl1 ) tree.mc_truth_tzq_Zl1_status = Zl1->status();
   if( Zl2 ) tree.mc_truth_tzq_Zl2_status = Zl2->status();
   if( Ztau1 ) tree.mc_truth_tzq_Ztau1_status = Ztau1->status();
   if( Ztau2 ) tree.mc_truth_tzq_Ztau2_status = Ztau2->status();
   if( Ztaul1 ) tree.mc_truth_tzq_Ztaul1_status = Ztaul1->status();
   if( Ztaul2 ) tree.mc_truth_tzq_Ztaul2_status = Ztaul2->status();
   if( Ztaunu1 ) tree.mc_truth_tzq_Ztaunu1_status = Ztaunu1->status();
   if( Ztaunu2 ) tree.mc_truth_tzq_Ztaunu2_status = Ztaunu2->status();
   if( Ztaunutau1 ) tree.mc_truth_tzq_Ztaunutau1_status = Ztaunutau1->status();
   if( Ztaunutau2 ) tree.mc_truth_tzq_Ztaunutau2_status = Ztaunutau2->status();
   
   if( t ) tree.mc_truth_tzq_t_status = t->status();
   if( tb ) tree.mc_truth_tzq_tb_status = tb->status();
   if( tW ) tree.mc_truth_tzq_tW_status = tW->status();
   if( tWnu ) tree.mc_truth_tzq_tWnu_status = tWnu->status();
   if( tWnutau ) tree.mc_truth_tzq_tWnutau_status = tWnutau->status();
   if( tWl ) tree.mc_truth_tzq_tWl_status = tWl->status();
   if( tWtau ) tree.mc_truth_tzq_tWtau_status = tWtau->status();
   if( tWtaunu ) tree.mc_truth_tzq_tWtaunu_status = tWtaunu->status();
   if( tWtaunutau ) tree.mc_truth_tzq_tWtaunutau_status = tWtaunutau->status();
   if( tWtaul ) tree.mc_truth_tzq_tWtaul_status = tWtaul->status();
   if( tWq1 ) tree.mc_truth_tzq_tWq1_status = tWq1->status();
   if( tWq2 ) tree.mc_truth_tzq_tWq2_status = tWq2->status();
   
   if( j1 ) tree.mc_truth_tzq_j1_status = j1->status();
   if( j2 ) tree.mc_truth_tzq_j2_status = j2->status();
   if( j3 ) tree.mc_truth_tzq_j3_status = j3->status();
}

reco::GenParticle* MCTruth::getUnique(const reco::GenParticle* p,bool verbose)
{
   reco::GenParticle *pcur = const_cast<reco::GenParticle*>(p);
   
   if( verbose )
     {	
	std::cout << "---------b--------" << std::endl;
	std::cout << "INITIAL = " << pcur->pdgId() << " " << pcur->status() << std::endl;
     }
   
   while( 1 )
     {
	bool foundDupl = false;

	const reco::GenParticleRefVector& daughterRefs = pcur->daughterRefVector();
	for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr)
	  {	     
	     if( idr->isAvailable() )
	       {		  
		  const reco::GenParticleRef& genParticle = (*idr);
		  const reco::GenParticle *d = genParticle.get();

		  if( d )
		    {			       
		       if( fabs(d->pdgId()) != 15 && d->status() == 2 ) continue;
		       if( verbose )
			 {		  
			    std::cout << "current: " << d->pdgId() << " " << d->status() << std::endl;
			    std::cout << "pcur: " << pcur->pdgId() << " " << pcur->status() << std::endl;
			 }
		       
		       if( d->pdgId() == pcur->pdgId() )
			 {
			    pcur = const_cast<reco::GenParticle*>(d);
			    foundDupl = true;
		       
			    if( verbose )
			      {		       
				 std::cout << "Found duplicate, switch to it" << std::endl;
				 std::cout << "Number of children = " << pcur->numberOfDaughters() << std::endl;
			      }		  
			 }
		    }
		  else break; // the world is fcked up in this case
	       }
	     else break;
	  }
	
	if( !foundDupl ) break;
     }   
   
   if( verbose )
     {   
	std::cout << "FINAL = " << pcur->pdgId() << " " << pcur->status() << std::endl;
	std::cout << "---------e--------" << std::endl;
     }
      
   return pcur;
}

void MCTruth::p4toTLV(reco::Particle::LorentzVector vp4,TLorentzVector& tlv)
{
   return tlv.SetPxPyPzE(vp4.px(),vp4.py(),vp4.pz(),vp4.energy());
}
