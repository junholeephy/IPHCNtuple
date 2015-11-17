#include "../include/TTbarHiggsMultileptonAnalysis.h"
#include "TSystem.h"

TTbarHiggsMultileptonAnalysis::TTbarHiggsMultileptonAnalysis() 
{
   _printLHCO_MC = false;
   _processLHCO_MC = -1;
   
   _printLHCO_RECO = false;
   _processLHCO_RECO = -1;

}

void TTbarHiggsMultileptonAnalysis::InitLHCO(int process_MC, int process_RECO) 
{  
   _printLHCO_MC = true;
   _processLHCO_MC = process_MC;
   fout_MC.open("LHCO_MC.txt");
   
   _printLHCO_RECO = true;
   _processLHCO_RECO = process_RECO;
   fout_RECO.open("LHCO_RECO.txt");
   
   fline00 = "#   typ	  eta	 phi	   pt  jmass  ntrk  btag   had/em  dummy dummy";
   del = "    ";
   trig = "8";
   
}

TTbarHiggsMultileptonAnalysis::TTbarHiggsMultileptonAnalysis(TString inputfilename, TTree *tree, TString theSampleName) 
{    
    gSystem->Load("libNtuple.so");

    if (tree == 0) 
    {
        TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(inputfilename.Data());
        if (!f || !f->IsOpen()) 
        {
            f = TFile::Open(inputfilename.Data(), "READ");
        }
        f->GetObject("Nt",tree);
    }
    Init(tree);

    theHistoManager = new HistoManager();

    sampleName = theSampleName;

    std::string foutlog = "output.txt";
    fevc = fopen(foutlog.c_str(),"w");

   outputfile = new TFile("output.root", "recreate");

}

void TTbarHiggsMultileptonAnalysis::createHistograms()
{    
 
   outputfile->cd();

    theHistoManager->addHisto("CutFlow",      "noSel", "emu", sampleName.Data(), 10, 0, 10);

    theHistoManager->addHisto("MuonPt",      "noSel", "emu", sampleName.Data(), 100, 0, 200);
    theHistoManager->addHisto("MuonEta",     "noSel", "emu", sampleName.Data(), 100, -3, 3);

    theHistoManager->addHisto("ElectronPt",  "noSel", "emu", sampleName.Data(), 100, 0, 200);
    theHistoManager->addHisto("ElectronEta", "noSel", "emu", sampleName.Data(), 100, -3, 3);

    theHistoManager->addHisto("JetPt",       "noSel", "emu", sampleName.Data(), 100, 0, 200);
    theHistoManager->addHisto("JetEta",      "noSel", "emu", sampleName.Data(), 100, -3, 3);

    initializeOutputTree();
}

void TTbarHiggsMultileptonAnalysis::writeHistograms()
{  
    outputfile->cd();

    std::vector<TH1F*> the1DHisto =  theHistoManager->getHisto1D_list();
    std::vector<TH2F*> the2DHisto =  theHistoManager->getHisto2D_list();

    for(unsigned int i=0; i<the1DHisto.size(); i++)  the1DHisto[i]->Write();
    for(unsigned int i=0; i<the2DHisto.size(); i++)  the2DHisto[i]->Write();

    tOutput->Write();
         
}

void TTbarHiggsMultileptonAnalysis::initializeOutputTree()
{

  outputfile->cd();
  tOutput = new TTree("Tree", "Tree");

  tOutput->Branch("mc_event",&mc_event,"mc_event/I");
  tOutput->Branch("mc_weight",&mc_weight,"mc_weight/F");
  tOutput->Branch("mc_3l_category",&mc_3l_category,"mc_3l_category/I");
  tOutput->Branch("mc_ttbar_decay",&mc_ttbar_decay,"mc_ttbar_decay/I");
  tOutput->Branch("mc_boson_decay",&mc_boson_decay,"mc_boson_decay/I");
  tOutput->Branch("mc_ttZhypAllowed",&mc_ttZhypAllowed,"mc_ttZhypAllowed/I");
  tOutput->Branch("mc_nJets25",&mc_nJets25,"mc_nJets25/I");
  tOutput->Branch("mc_nBtagJets25",&mc_nBtagJets25,"mc_nBtagJets25/I");
  tOutput->Branch("mc_nNonBtagJets25",&mc_nNonBtagJets25,"mc_nNonBtagJets25/I");

  tOutput->Branch("multilepton_Bjet1_Id",&multilepton_Bjet1_Id,"multilepton_Bjet1_Id/I");
  tOutput->Branch("multilepton_Bjet1_P4","TLorentzVector",&multilepton_Bjet1_P4);
  tOutput->Branch("multilepton_Bjet2_Id",&multilepton_Bjet2_Id,"multilepton_Bjet2_Id/I");
  tOutput->Branch("multilepton_Bjet2_P4","TLorentzVector",&multilepton_Bjet2_P4);
  tOutput->Branch("multilepton_Lepton1_Id",&multilepton_Lepton1_Id,"multilepton_Lepton1_Id/I");
  tOutput->Branch("multilepton_Lepton1_P4","TLorentzVector",&multilepton_Lepton1_P4);
  tOutput->Branch("multilepton_Lepton2_Id",&multilepton_Lepton2_Id,"multilepton_Lepton2_Id/I");
  tOutput->Branch("multilepton_Lepton2_P4","TLorentzVector",&multilepton_Lepton2_P4);
  tOutput->Branch("multilepton_Lepton3_Id",&multilepton_Lepton3_Id,"multilepton_Lepton3_Id/I");
  tOutput->Branch("multilepton_Lepton3_P4","TLorentzVector",&multilepton_Lepton3_P4);
  tOutput->Branch("multilepton_JetHighestPt1_Id",&multilepton_JetHighestPt1_Id,"multilepton_JetHighestPt1_Id/I");
  tOutput->Branch("multilepton_JetHighestPt1_P4","TLorentzVector",&multilepton_JetHighestPt1_P4);
  tOutput->Branch("multilepton_JetHighestPt2_Id",&multilepton_JetHighestPt2_Id,"multilepton_JetHighestPt2_Id/I");
  tOutput->Branch("multilepton_JetHighestPt2_P4","TLorentzVector",&multilepton_JetHighestPt2_P4);
  tOutput->Branch("multilepton_JetClosestMw1_Id",&multilepton_JetClosestMw1_Id,"multilepton_JetClosestMw1_Id/I");
  tOutput->Branch("multilepton_JetClosestMw1_P4","TLorentzVector",&multilepton_JetClosestMw1_P4);
  tOutput->Branch("multilepton_JetClosestMw2_Id",&multilepton_JetClosestMw2_Id,"multilepton_JetClosestMw2_Id/I");
  tOutput->Branch("multilepton_JetClosestMw2_P4","TLorentzVector",&multilepton_JetClosestMw2_P4);
  tOutput->Branch("multilepton_JetLowestMjj1_Id",&multilepton_JetLowestMjj1_Id,"multilepton_JetLowestMjj1_Id/I");
  tOutput->Branch("multilepton_JetLowestMjj1_P4","TLorentzVector",&multilepton_JetLowestMjj1_P4);
  tOutput->Branch("multilepton_JetLowestMjj2_Id",&multilepton_JetLowestMjj2_Id,"multilepton_JetLowestMjj2_Id/I");
  tOutput->Branch("multilepton_JetLowestMjj2_P4","TLorentzVector",&multilepton_JetLowestMjj2_P4);
  tOutput->Branch("multilepton_mET","TLorentzVector",&multilepton_mET);
  tOutput->Branch("multilepton_Ptot","TLorentzVector",&multilepton_Ptot);
 
  return;
}

void TTbarHiggsMultileptonAnalysis::selectBjets(std::string BjetSel, int* ibsel1, int* ibsel2){

  //Assumes there are at least 2 b-tagged jets in the event
  int ib1=-1, ib2=-1;

  if (BjetSel=="HighestBtagDiscrim"){
    float btag_max=-999, btag_max2=-999;
    for (unsigned int ib=0; ib<vSelectedJets.size(); ib++){
      if (vSelectedJets.at(ib).CSVv2()<0.423) continue;
      if (vSelectedJets.at(ib).CSVv2()>btag_max){
        btag_max2 = btag_max;
        ib2 = ib1;
        btag_max = vSelectedJets.at(ib).CSVv2();
        ib1 = ib;
      }
      if (vSelectedJets.at(ib).CSVv2()<btag_max && vSelectedJets.at(ib).CSVv2()>btag_max2){
        btag_max2 = vSelectedJets.at(ib).CSVv2();
        ib2 = ib;
      }
    }
  }
  if (BjetSel=="BtagHighestPt"){
    float pt_max=0, pt_max2=0;
    for (unsigned int ib=0; ib<vSelectedJets.size(); ib++){
      if (vSelectedJets.at(ib).CSVv2()<0.423) continue;
      if (vSelectedJets.at(ib).pt()>pt_max){
        pt_max2 = pt_max;
        ib2 = ib1;
        pt_max = vSelectedJets.at(ib).pt();
        ib1 = ib;
      }
      if (vSelectedJets.at(ib).pt()<pt_max && vSelectedJets.at(ib).pt()>pt_max2){
        pt_max2 = vSelectedJets.at(ib).pt();
        ib2 = ib;
      }
    }
  }

  *ibsel1 = ib1;
  *ibsel2 = ib2;

}

void TTbarHiggsMultileptonAnalysis::fillOutputTree(){

  if (vSelectedLeptons.size()!=3 || vSelectedBTagJets.size()<2 || vSelectedJets.size()<4 ) return; 

  mc_weight = vEvent->at(0).mc_weight();
  //std::cout << "fillOutputTree mc_weight="<<mc_weight<<std::endl;

  multilepton_Lepton1_P4 = vSelectedLeptons.at(0).p4();
  multilepton_Lepton1_Id = vSelectedLeptons.at(0).id();
  multilepton_Lepton2_P4 = vSelectedLeptons.at(1).p4();
  multilepton_Lepton2_Id = vSelectedLeptons.at(1).id();
  multilepton_Lepton3_P4 = vSelectedLeptons.at(2).p4();
  multilepton_Lepton3_Id = vSelectedLeptons.at(2).id();

  //Choosing 2 b-jets
  TLorentzVector Bjet1, Bjet2; 
  int ib1=-1, ib2=-1;
  selectBjets("HighestBtagDiscrim", &ib1, &ib2);
  Bjet1.SetPtEtaPhiE(vSelectedJets.at(ib1).pt(), vSelectedJets.at(ib1).eta(), vSelectedJets.at(ib1).phi(), vSelectedJets.at(ib1).E());
  Bjet2.SetPtEtaPhiE(vSelectedJets.at(ib2).pt(), vSelectedJets.at(ib2).eta(), vSelectedJets.at(ib2).phi(), vSelectedJets.at(ib2).E());

  multilepton_Bjet1_P4 = Bjet1;
  multilepton_Bjet1_Id = 5;
  multilepton_Bjet2_P4 = Bjet2;
  multilepton_Bjet2_Id = 5;

  //Choose 2 jets
      TLorentzVector Pjet1, Pjet2;
      float pt_max=0, pt_max2=0; int ij1=-1, ij2=-1;
      float diffmass_min = 10000, mass_min = 10000; int ik1=-1, ik2=-1, il1=-1, il2=-1;
      for (unsigned int ij=0; ij<vSelectedJets.size(); ij++){
	if (ij==ib1 || ij==ib2) continue;
        if (vSelectedJets.at(ij).pt() > pt_max ) {
           pt_max2 = pt_max;
           ij2 = ij1;
           pt_max = vSelectedJets.at(ij).pt();
           ij1 = ij;
         } 
         if (vSelectedJets.at(ij).pt() < pt_max && vSelectedJets.at(ij).pt() > pt_max2){
           pt_max2 = vSelectedJets.at(ij).pt(); 
           ij2 = ij; 
         } 
         for (unsigned int ik=0; ik<vSelectedJets.size(); ik++){
           if (ik==ij) continue;
	   if (ik==ib1 || ik==ib2) continue;
	   Pjet1.SetPtEtaPhiE(vSelectedJets.at(ij).pt(), vSelectedJets.at(ij).eta(), vSelectedJets.at(ij).phi(), vSelectedJets.at(ij).E());
	   Pjet2.SetPtEtaPhiE(vSelectedJets.at(ik).pt(), vSelectedJets.at(ik).eta(), vSelectedJets.at(ik).phi(), vSelectedJets.at(ik).E()); 
           if (TMath::Abs((Pjet1+Pjet2).M()-80.419)<diffmass_min){
             ik1=ij;
             ik2=ik;
             diffmass_min = TMath::Abs((Pjet1+Pjet2).M()-80.419);
           } 
           if ((Pjet1+Pjet2).M()<mass_min){
             il1=ij;
             il2=ik;
             mass_min = (Pjet1+Pjet2).M();
           } 
         } 
      }  
      if (ij1!=-1 && ij2!=-1) {
	multilepton_JetHighestPt1_Id = 1;
	multilepton_JetHighestPt2_Id = 1;
	multilepton_JetHighestPt1_P4.SetPtEtaPhiE(vSelectedJets.at(ij1).pt(), vSelectedJets.at(ij1).eta(), vSelectedJets.at(ij1).phi(), vSelectedJets.at(ij1).E());
        multilepton_JetHighestPt2_P4.SetPtEtaPhiE(vSelectedJets.at(ij2).pt(), vSelectedJets.at(ij2).eta(), vSelectedJets.at(ij2).phi(), vSelectedJets.at(ij2).E());
      }
      if (ik1!=-1 && ik2!=-1){
	multilepton_JetClosestMw1_Id = 2;
	multilepton_JetClosestMw2_Id = 2;
	multilepton_JetClosestMw1_P4.SetPtEtaPhiE(vSelectedJets.at(ik1).pt(), vSelectedJets.at(ik1).eta(), vSelectedJets.at(ik1).phi(), vSelectedJets.at(ik1).E());
	multilepton_JetClosestMw2_P4.SetPtEtaPhiE(vSelectedJets.at(ik2).pt(), vSelectedJets.at(ik2).eta(), vSelectedJets.at(ik2).phi(), vSelectedJets.at(ik2).E());
      }
      if (il1!=-1 && il2!=-1){
	multilepton_JetLowestMjj1_Id = 3;
	multilepton_JetLowestMjj2_Id = 3;
	multilepton_JetLowestMjj1_P4.SetPtEtaPhiE(vSelectedJets.at(il1).pt(), vSelectedJets.at(il1).eta(), vSelectedJets.at(il1).phi(), vSelectedJets.at(il1).E());
        multilepton_JetLowestMjj2_P4.SetPtEtaPhiE(vSelectedJets.at(il2).pt(), vSelectedJets.at(il2).eta(), vSelectedJets.at(il2).phi(), vSelectedJets.at(il2).E());
      }

  multilepton_mET.SetPtEtaPhiE(vEvent->at(0).metpt(), 0, vEvent->at(0).metphi(), vEvent->at(0).metpt());

  if ( vSelectedLeptons.at(0).charge()==vSelectedLeptons.at(1).charge() && vSelectedLeptons.at(1).charge()==vSelectedLeptons.at(2).charge() ) mc_ttZhypAllowed =-1;
  else if (  ( vSelectedLeptons.at(0).id() == -vSelectedLeptons.at(1).id() ) 
          || ( vSelectedLeptons.at(0).id() == -vSelectedLeptons.at(2).id() ) 
          || ( vSelectedLeptons.at(1).id() == -vSelectedLeptons.at(2).id() ))
	mc_ttZhypAllowed = 1;
  else mc_ttZhypAllowed = 0; 

  mc_nJets25 = vSelectedJets.size();
  mc_nBtagJets25 = vSelectedBTagJets.size();
  mc_nNonBtagJets25 = vSelectedNonBTagJets.size();

  //std::cout << "mc_ttZhypAllowed="<<mc_ttZhypAllowed<<" mc_nJets25="<<mc_nJets25<<" mc_nBtagJets25="<<mc_nBtagJets25<<" mc_nNonBtagJets25="<<mc_nNonBtagJets25<<std::endl;

  tOutput->Fill();
  
  //if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);

}

void TTbarHiggsMultileptonAnalysis::Init(TTree *tree)
{
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;

    vElectron = new std::vector<Electron>();
    vMuon     = new std::vector<Muon>();
    vEvent    = new std::vector<Event>();
    vJet      = new std::vector<Jet>();
    vTruth    = new std::vector<Truth>();

    fChain->SetBranchAddress("Electron", &vElectron);
    fChain->SetBranchAddress("Muon",     &vMuon);
    fChain->SetBranchAddress("Jet",      &vJet);
    fChain->SetBranchAddress("Event",    &vEvent);
    fChain->SetBranchAddress("Truth",    &vTruth);
}

void TTbarHiggsMultileptonAnalysis::Loop()
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntries();

    std::cout << "Number of input events = " << nentries << std::endl;

    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
        Long64_t ientry = fChain->LoadTree(jentry);
        if (ientry < 0) break;

        if(jentry%10000 == 0) std::cout << "number of processed events " << jentry << std::endl;
        //std::cout << "number of processed events " << jentry << std::endl;

        //if(jentry > 100000) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        //int nPart = vTruth->at(0).mc_truth_n();
        //std::vector<int> labTruth = vTruth->at(0).mc_truth_label();


        float theweight = vEvent->at(0).mc_weight();
	
        //std::vector<Muon>     vSelectedMuons;
        //std::vector<Electron> vSelectedElectrons;
        //std::vector<Lepton>   vSelectedLeptons;
        //std::vector<Jet>      vSelectedNonBTagJets;
        //std::vector<Jet>      vSelectedBTagJets;
        //std::vector<Jet>      vSelectedJets;
	vSelectedMuons.clear();
	vSelectedElectrons.clear();
	vSelectedLeptons.clear();
	vSelectedNonBTagJets.clear();
	vSelectedBTagJets.clear();
	vSelectedJets.clear();

        //std::cout << "Jusqu'ici tout va bien 2" << std::endl;

        for(unsigned int imuon=0; imuon < vMuon->size() ; imuon++)
        {
            // ADD CONDITION ON LEP MVA!

            //if( !vMuon->at(imuon).passPtEta() ) continue;
            //if( !vMuon->at(imuon).isLoose() )   continue;

            //	     theHistoManager->fillHisto("MuonPt",  "noSel", "emu", sampleName.Data(),  vMuon->at(imuon).pt(), theweight);
            //	     theHistoManager->fillHisto("MuonEta", "noSel", "emu", sampleName.Data(),  vMuon->at(imuon).eta(), theweight);

            vSelectedMuons.push_back(vMuon->at(imuon));

            Lepton l; l.setLepton(&vMuon->at(imuon),imuon,0);
            vSelectedLeptons.push_back(l);
        }     

        //std::cout << "Jusqu'ici tout va bien 3" << std::endl;

        for(unsigned int ielectron=0; ielectron < vElectron->size() ; ielectron++)
        {
            // ADD CONDITION ON LEP MVA!

            //if( !vElectron->at(ielectron).passPtEta() ) continue;
            //if( !vElectron->at(ielectron).isLoose() )   continue;

            //	     theHistoManager->fillHisto("ElectronPt",  "noSel", "emu", sampleName.Data(),  vElectron->at(ielectron).pt(), theweight);
            //	     theHistoManager->fillHisto("ElectronEta", "noSel", "emu", sampleName.Data(),  vElectron->at(ielectron).eta(), theweight);

            vSelectedElectrons.push_back(vElectron->at(ielectron));	     

            Lepton l; l.setLepton(&vElectron->at(ielectron),ielectron,1);
            vSelectedLeptons.push_back(l);
        }  

        //std::cout << "Jusqu'ici tout va bien 4" << std::endl;

        int nLooseBJets  = 0;
        int nMediumBJets = 0;

        for(unsigned int ijet=0; ijet < vJet->size() ; ijet++)
        {
            //if(  vJet->at(ijet).pt() < 25   )         continue;
            //if(  fabs(vJet->at(ijet).eta()) > 2.4 )   continue;

            //	     theHistoManager->fillHisto("JetPt",  "noSel", "emu", sampleName.Data(),  vJet->at(ijet).pt(), theweight);
            //	     theHistoManager->fillHisto("JetEta", "noSel", "emu", sampleName.Data(),  vJet->at(ijet).eta(), theweight);

            if( vJet->at(ijet).CSVv2() > 0.423 ) nLooseBJets  = nLooseBJets  + 1;
            if( vJet->at(ijet).CSVv2() > 0.814 ) nMediumBJets = nMediumBJets + 1;

            if(vJet->at(ijet).CSVv2() >= 0.423 ) vSelectedBTagJets.push_back(vJet->at(ijet));
            else                                 vSelectedNonBTagJets.push_back(vJet->at(ijet));

            vSelectedJets.push_back(vJet->at(ijet));
        }

        //	theHistoManager->fillHisto("CutFlow", "noSel", "emu", sampleName.Data(),  0, theweight);

        //bool nLep = (vSelectedLeptons.size() >= 2);

       //std::cout << "Jusqu'ici tout va bien 5" << std::endl; 

        // #################################
        // # Three leptons event selection #
        // #################################

        bool nLep        = ( vSelectedLeptons.size()                                  == 3 );
        bool nJets       = ( (vSelectedBTagJets.size() + vSelectedNonBTagJets.size()) >= 4 );
        bool nLooseBtag  = ( nLooseBJets                                              >= 2 );
        bool nMediumBtag = ( nMediumBJets                                             >= 1 );

        //if( nLep )
        //{	 
            // PrintEventList(vSelectedLeptons,vSelectedJets);
            // theHistoManager->fillHisto("CutFlow", "noSel", "emu", sampleName.Data(),  1, theweight);

        //    if( (vSelectedBTagJets.size() + vSelectedNonBTagJets.size()) >= 4)
        //    {
        //        theHistoManager->fillHisto("CutFlow", "noSel", "emu", sampleName.Data(),  2, theweight);

        //        if(vSelectedBTagJets.size() >=1)
        //        {		     
        //            theHistoManager->fillHisto("CutFlow", "noSel", "emu", sampleName.Data(),  3, theweight);
        //        }//NBjet selection
        //    }//end NJet selection
        //}//end nlepton selection

        ThreeLeptonSelection(vSelectedLeptons, vSelectedJets, vSelectedBTagJets, vSelectedNonBTagJets, nLooseBJets, nMediumBJets, jentry);

        /*std::cout << "Number of selected leptons:    " << vSelectedLeptons.size()     << std::endl;
        std::cout << "Number of selected jets   :    " << (vSelectedBTagJets.size() + vSelectedNonBTagJets.size()) << std::endl;
        std::cout << "Number of loose b jets    :    " << nLooseBJets                 << std::endl;
        std::cout << "Number of medium b jets   :    " << nMediumBJets                << std::endl;

        // #################################
        // # Three leptons event selection #
        // #################################

        if ( nLep && nJets && (nLooseBtag || nMediumBtag) )
        {

            // #############################
            // # FILLING MULTILEPTON CLASS #
            // #############################

            std::cout << "******** SELECTED EVENT OF THE THREE LEPTON CATEGORY! ********" << std::endl;
            std::cout << "Number of selected leptons:    " << vSelectedLeptons.size()     << std::endl;
            std::cout << "Number of selected jets   :    " << (vSelectedBTagJets.size() + vSelectedNonBTagJets.size()) << std::endl;
            std::cout << "Number of loose b jets    :    " << nLooseBJets                 << std::endl;
            std::cout << "Number of medium b jets   :    " << nMediumBJets                << std::endl;
        }*/
        
	if (_printLHCO_MC && ThreeLeptonSelection_MC()) PrintLHCOforMadweight_MC(jentry);
	
    }

}

void TTbarHiggsMultileptonAnalysis::PrintEventList(std::vector<Lepton> leptons, std::vector<Jet> jets)
{
    int run  = vEvent->at(0).run();
    int id   = vEvent->at(0).id();
    int lumi = vEvent->at(0).lumi();

    float metpt  = vEvent->at(0).metpt();
    float metphi = vEvent->at(0).metphi();

    int   l1id  = leptons.at(0).id();
    float l1pt  = leptons.at(0).pt();
    float l1eta = leptons.at(0).eta();
    float l1phi = leptons.at(0).phi();

    int   l2id  = leptons.at(1).id();
    float l2pt  = leptons.at(1).pt();
    float l2eta = leptons.at(1).eta();
    float l2phi = leptons.at(1).phi();

    int   njets = jets.size();

    fprintf(fevc,"%6d %6d %10d  %+2d  %6.2f %+4.2f %+4.2f   %+2d  %6.2f %+4.2f %+4.2f    %6.1f  %+4.2f    %d \n",
            run, lumi, id,
            l1id, l1pt, l1eta, l1phi,
            l2id, l2pt, l2eta, l2phi,
            metpt, metphi,
            njets);
}

void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection(std::vector<Lepton> vSelectedLeptons, 
                                                         std::vector<Jet>    vSelectedJets, 
                                                         std::vector<Jet>    vSelectedBTagJets, 
                                                         std::vector<Jet>    vSelectedNonBTagJets,
                                                         int                 nLooseBJets, 
                                                         int                 nMediumBJets, 
							 int evt)
{
    //std::cout << "Number of selected leptons:    " << vSelectedLeptons.size()       << std::endl;
    //std::cout << "Number of selected jets   :    " << (vSelectedBTagJets.size() 
    //                                                 + vSelectedNonBTagJets.size()) << std::endl;
    //std::cout << "Number of loose b jets    :    " << nLooseBJets                   << std::endl;
    //std::cout << "Number of medium b jets   :    " << nMediumBJets                  << std::endl;

    // #################################
    // # Three leptons event selection #
    // #################################

    // three lepton selection
    bool nLep        = ( vSelectedLeptons.size()                                  == 3 );
    bool nJets       = ( (vSelectedBTagJets.size() + vSelectedNonBTagJets.size()) >= 4 );


    bool pass_OSSF = true;
    if (nLep) 
    { 
        if ( ( ( vSelectedLeptons.at(0).id() == -vSelectedLeptons.at(1).id() ) && ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M() - 91.188 < 10 ) )
          || ( ( vSelectedLeptons.at(0).id() == -vSelectedLeptons.at(2).id() ) && ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(2).p4() ).M() - 91.188 < 10 ) )
          || ( ( vSelectedLeptons.at(1).id() == -vSelectedLeptons.at(2).id() ) && ( ( vSelectedLeptons.at(1).p4() + vSelectedLeptons.at(2).p4() ).M() - 91.188 < 10 ) ) )
        { pass_OSSF = false ;}
    }

    // common selection
    bool leading_lep_pt = 0;
    if (nLep) leading_lep_pt = ( vSelectedLeptons.at(0).pt() > 20 );
    bool following_lep_pt = 0;
    if (nLep) following_lep_pt = ( vSelectedLeptons.at(1).pt() > 10 );

    bool passMll12Gt12;
    if (nLep) passMll12Gt12  = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M()  > 12
                              && ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(2).p4() ).M()  > 12
                              && ( vSelectedLeptons.at(1).p4() + vSelectedLeptons.at(2).p4() ).M()  > 12 );

    bool nLooseBtag  = ( nLooseBJets                                              >= 2 );
    bool nMediumBtag = ( nMediumBJets                                             >= 1 );
/*
    std::cout << "nLep: "              << nLep 
              << " nJets: "            << nJets 
              << " pass_OSSF: "        << pass_OSSF 
              << " leading_lep_pt: "   << leading_lep_pt 
              << " following_lep_pt: " << following_lep_pt
              << " passMll12Gt12: "    << passMll12Gt12 
              << " nLooseBtag: "       << nLooseBtag 
              << " nMediumBtag: "      << nMediumBtag << std::endl; 
*/
    bool passCharge = true;
    if (vSelectedLeptons.size()>=3 && vSelectedLeptons.at(0).charge()==vSelectedLeptons.at(1).charge() && vSelectedLeptons.at(1).charge()==vSelectedLeptons.at(2).charge() ) passCharge = false;

    if ( nLep && nJets && passCharge && pass_OSSF && leading_lep_pt && following_lep_pt && passMll12Gt12 && (nLooseBtag || nMediumBtag) )
    {

        // #############################
        // # FILLING MULTILEPTON CLASS #
        // #############################
/*
        std::cout << "******** SELECTED EVENT OF THE THREE LEPTON CATEGORY! ********" << std::endl;
        std::cout << "Number of selected leptons:    " << vSelectedLeptons.size()     << std::endl;
        std::cout << "Number of selected jets   :    " << (vSelectedBTagJets.size() + vSelectedNonBTagJets.size()) << std::endl;
        std::cout << "Number of loose b jets    :    " << nLooseBJets                 << std::endl;
        std::cout << "Number of medium b jets   :    " << nMediumBJets                << std::endl;
*/

      fillOutputTree();
      
      if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);
    
    }
}

bool TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_MC() 
{ 
  bool sel_MC = true;
  
  //std::cout <<" sel_MC 1"<< std::endl;
   
  //Check decays and presence of genjets
  if (!((vTruth->at(0).ttbar_decay()==1 && vTruth->at(0).boson_decay()==0) || (vTruth->at(0).ttbar_decay()==2 && vTruth->at(0).boson_decay()==1))) 
  { 
    sel_MC = false; 
    return sel_MC;}
    
  if (vTruth->at(0).JetsHighestPt_phi().size()<2 || vTruth->at(0).JetsClosestMw_phi().size()<2 || vTruth->at(0).JetsLowestMjj_phi().size()<2) 
  { 
    sel_MC = false; 
    return sel_MC;}
  
  //std::cout <<" sel_MC 212 "<<vTruth->at(0).Leptons_id().size()<< std::endl;
  
  //SFOS 
  if (!(vTruth->at(0).Leptons_id().at(0)==-vTruth->at(0).Leptons_id().at(1) || 
        vTruth->at(0).Leptons_id().at(1)==-vTruth->at(0).Leptons_id().at(2) ||
	vTruth->at(0).Leptons_id().at(0)==-vTruth->at(0).Leptons_id().at(2)   )) sel_MC = false; 
 
  //std::cout <<" sel_MC 3"<< std::endl;
  	 
  //pt , eta of leptons	
  if (!(vTruth->at(0).Leptons_pt().at(0)> 10 && 
        vTruth->at(0).Leptons_pt().at(1)> 10 &&
	vTruth->at(0).Leptons_pt().at(2)> 10   )) sel_MC = false; 
  
  //std::cout <<" sel_MC 4"<< std::endl;

  if (!(fabs(vTruth->at(0).Leptons_eta().at(0)) <2.5 && 
        fabs(vTruth->at(0).Leptons_eta().at(1)) <2.5 &&
        fabs(vTruth->at(0).Leptons_eta().at(2)) <2.5   )) sel_MC = false; 
  
  //std::cout <<" sel_MC 5"<< std::endl;
	
  //lead. lepton
  if (!(vTruth->at(0).Leptons_pt().at(0) > 20 || 
        vTruth->at(0).Leptons_pt().at(1) > 20 ||
	vTruth->at(0).Leptons_pt().at(2) > 20  )) sel_MC = false; 
	
  TLorentzVector Lep1;
  Lep1.SetPtEtaPhiE(vTruth->at(0).Leptons_pt().at(0),  vTruth->at(0).Leptons_eta().at(0), vTruth->at(0).Leptons_phi().at(0), vTruth->at(0).Leptons_E().at(0));
  TLorentzVector Lep2;
  Lep2.SetPtEtaPhiE(vTruth->at(0).Leptons_pt().at(1),  vTruth->at(0).Leptons_eta().at(1), vTruth->at(0).Leptons_phi().at(1), vTruth->at(0).Leptons_E().at(1));
  TLorentzVector Lep3;
  Lep3.SetPtEtaPhiE(vTruth->at(0).Leptons_pt().at(2),  vTruth->at(0).Leptons_eta().at(2), vTruth->at(0).Leptons_phi().at(2), vTruth->at(0).Leptons_E().at(2));
 
  
  if ( !(( Lep1+Lep2 ).M()  > 12 && ( Lep1+Lep3 ).M()  > 12 && ( Lep2+Lep3 ).M()  > 12 )) sel_MC = false;
  
  return sel_MC;
   
}


void TTbarHiggsMultileptonAnalysis::PrintLHCOforMadweight_MC(int evt)
{  

  //std::cout <<" _MC 1"<< std::endl;

  /*if( proc<-1 || proc > 6 )f
  {
      std::cout << "proc can only take following values: -1,1,2,3,4,5,6" << std::endl;
      std::cout << "3l final state specific" << std::endl;    
      std::cout << "1,2,3,4: specific to the ttH final state, cf patches provided to madweight" << std::endl;
      std::cout << "5: specific to the ttZ with l+l-l+ final state" << std::endl;
      std::cout << "6: specific to the ttZ with l+l-l- final state" << std::endl;
      std::cout << "-1: no selection on the final state applied" << std::endl;

      return;
   }*/
      	
 
  fline0 = "0      " + std::string(Form("%d",evt)) + "     " + trig;
  
  int nobj = 1;	      
  
  int multilepton_Lepton1_Id_LHCO = -666;
  int multilepton_Lepton2_Id_LHCO = -666;
  int multilepton_Lepton3_Id_LHCO = -666;
 
  //
  // LHCO lepton ID convention
  //  
  if (abs(vTruth->at(0).Leptons_id().at(0))==11) multilepton_Lepton1_Id_LHCO = 1 ;
  else if (abs(vTruth->at(0).Leptons_id().at(0))==13) multilepton_Lepton1_Id_LHCO = 2 ;
  if (abs(vTruth->at(0).Leptons_id().at(1))==11) multilepton_Lepton2_Id_LHCO = 1 ;
  else if (abs(vTruth->at(0).Leptons_id().at(1))==13) multilepton_Lepton2_Id_LHCO = 2 ;
  if (abs(vTruth->at(0).Leptons_id().at(2))==11) multilepton_Lepton3_Id_LHCO = 1 ;
  else if (abs(vTruth->at(0).Leptons_id().at(2))==13) multilepton_Lepton3_Id_LHCO = 2 ;
 
  //
  // LHCO phi convention
  //
  float multilepton_Lepton1_phi       = Phi_0_2Pi(vTruth->at(0).Leptons_phi().at(0));
  float multilepton_Lepton2_phi       = Phi_0_2Pi(vTruth->at(0).Leptons_phi().at(1));
  float multilepton_Lepton3_phi       = Phi_0_2Pi(vTruth->at(0).Leptons_phi().at(2));
  float multilepton_Bjet1_phi         = Phi_0_2Pi(vTruth->at(0).Bjets_phi().at(0));
  float multilepton_Bjet2_phi         = Phi_0_2Pi(vTruth->at(0).Bjets_phi().at(1));	  
  float multilepton_JetHighestPt1_phi = Phi_0_2Pi(vTruth->at(0).JetsHighestPt_phi().at(0));
  float multilepton_JetHighestPt2_phi = Phi_0_2Pi(vTruth->at(0).JetsHighestPt_phi().at(1));
  float multilepton_JetClosestMw1_phi = Phi_0_2Pi(vTruth->at(0).JetsClosestMw_phi().at(0));
  float multilepton_JetClosestMw2_phi = Phi_0_2Pi(vTruth->at(0).JetsClosestMw_phi().at(1));
  float multilepton_JetLowestMjj1_phi = Phi_0_2Pi(vTruth->at(0).JetsLowestMjj_phi().at(0));
  float multilepton_JetLowestMjj2_phi = Phi_0_2Pi(vTruth->at(0).JetsLowestMjj_phi().at(1));
  
  //std::cout <<" _MC 22"<< std::endl;

  // l1
  std::string l1_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d", 
 					nobj,multilepton_Lepton1_Id_LHCO,vTruth->at(0).Leptons_eta().at(0),multilepton_Lepton1_phi,vTruth->at(0).Leptons_pt().at(0),0.0,vTruth->at(0).Leptons_id().at(0)/abs(vTruth->at(0).Leptons_id().at(0)),0,0,0,0));
  nobj++; 
  //std::cout <<" _MC 23"<< std::endl;
  // l2
  std::string l2_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d",
 					nobj,multilepton_Lepton2_Id_LHCO,vTruth->at(0).Leptons_eta().at(1),multilepton_Lepton2_phi,vTruth->at(0).Leptons_pt().at(1),0.0,vTruth->at(0).Leptons_id().at(1)/abs(vTruth->at(0).Leptons_id().at(1)),0,0,0,0));
  nobj++;	
  //std::cout <<" _MC 24"<< std::endl;
  // l3
  std::string l3_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d",
 					nobj,multilepton_Lepton3_Id_LHCO,vTruth->at(0).Leptons_eta().at(2),multilepton_Lepton3_phi,vTruth->at(0).Leptons_pt().at(2),0.0,vTruth->at(0).Leptons_id().at(2)/abs(vTruth->at(0).Leptons_id().at(2)),0,0,0,0));
  nobj++;
  //std::cout <<" _MC 3"<< std::endl;
			    
  //										    
  std::string j1_fline;
  std::string j2_fline;

  if ( _processLHCO_MC == 5 || _processLHCO_MC == 6 || _processLHCO_MC == 4 || _processLHCO_MC == 3 )
  {			     
    // j1
    j1_fline = std::string(Form("%d	  %d	 %.2f	    %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
 			      nobj,4,vTruth->at(0).JetsClosestMw_eta().at(0),multilepton_JetClosestMw1_phi,vTruth->at(0).JetsClosestMw_pt().at(0),0.0,1,0,0,0,0));
    nobj++;
 
    // j2
    j2_fline = std::string(Form("%d	  %d	  %.2f      %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
 			      nobj,4,vTruth->at(0).JetsClosestMw_eta().at(1),multilepton_JetClosestMw2_phi,vTruth->at(0).JetsClosestMw_pt().at(1),0.0,1,0,0,0,0));
    nobj++;
  }
  else if ( _processLHCO_MC == 1 || _processLHCO_MC == 2)
  {
    // j1
    j1_fline = std::string(Form("%d	  %d	  %.2f      %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
 			      nobj,4,vTruth->at(0).JetsLowestMjj_eta().at(0),multilepton_JetLowestMjj1_phi,vTruth->at(0).JetsLowestMjj_pt().at(0),0.0,1,0,0,0,0));
    nobj++;
 
    // j2
    j2_fline = std::string(Form("%d	  %d	   %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
 			      nobj,4,vTruth->at(0).JetsLowestMjj_eta().at(1),multilepton_JetLowestMjj2_phi,vTruth->at(0).JetsLowestMjj_pt().at(1),0.0,1,0,0,0,0));
    nobj++;
  }
  
  //
  TLorentzVector BJet1;
  BJet1.SetPtEtaPhiE(vTruth->at(0).Bjets_pt().at(0), vTruth->at(0).Bjets_eta().at(0), vTruth->at(0).Bjets_phi().at(0), vTruth->at(0).Bjets_E().at(0));
 
  TLorentzVector BJet2;
  BJet2.SetPtEtaPhiE(vTruth->at(0).Bjets_pt().at(1), vTruth->at(0).Bjets_eta().at(1), vTruth->at(0).Bjets_phi().at(1), vTruth->at(0).Bjets_E().at(1));
         

  // bj1
  std::string b1_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d",
 					nobj,4,vTruth->at(0).Bjets_eta().at(0),multilepton_Bjet1_phi,vTruth->at(0).Bjets_pt().at(0),BJet1.M(),1,2,0,0,0));
  nobj++;

 				     
  // bj2 
  std::string b2_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d",
 					nobj,4,vTruth->at(0).Bjets_eta().at(1),multilepton_Bjet2_phi,vTruth->at(0).Bjets_pt().at(1),BJet2.M(),1,2,0,0,0));
  nobj++;

  // met
  std::string met_fline = std::string(Form("%d     %d	  %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d	  %d",
 					 nobj,6,0.,Phi_0_2Pi(vTruth->at(0).metGen_phi()),vTruth->at(0).metGen_pt(),0.,0,0,0,0,0));
  nobj++;
  

 //    
 fout_MC << fline00   << std::endl;
 fout_MC << fline0    << std::endl;
 fout_MC << l1_fline  << std::endl;
 fout_MC << l2_fline  << std::endl;
 fout_MC << l3_fline  << std::endl;	  
 if (!(_processLHCO_MC == 7 || _processLHCO_MC == 8)) fout_MC << j1_fline  << std::endl;// don't print jets for ttW hypothesis
 if (!(_processLHCO_MC == 7 || _processLHCO_MC == 8)) fout_MC << j2_fline  << std::endl;// don't print jets for ttW hypothesis
 fout_MC << b1_fline  << std::endl;
 fout_MC << b2_fline  << std::endl;
 fout_MC << met_fline << std::endl;      
       
}

void TTbarHiggsMultileptonAnalysis::PrintLHCOforMadweight_RECO(int evt)
{
 
  if (vSelectedLeptons.size()!=3 || vSelectedBTagJets.size()<2 || vSelectedJets.size()<4 ) return; 

  fline0 = "0      " + std::string(Form("%d",evt)) + "     " + trig;
  
  int nobj = 1;	      
  
  int multilepton_Lepton1_Id_LHCO = -666;
  int multilepton_Lepton2_Id_LHCO = -666;
  int multilepton_Lepton3_Id_LHCO = -666;
 
  //
  // LHCO lepton ID convention
  //  
  if(abs(multilepton_Lepton1_Id)==11) multilepton_Lepton1_Id_LHCO = 1 ;
  else if(abs(multilepton_Lepton1_Id)==13) multilepton_Lepton1_Id_LHCO = 2 ;
  if(abs(multilepton_Lepton2_Id)==11) multilepton_Lepton2_Id_LHCO = 1 ;
  else if(abs(multilepton_Lepton2_Id)==13) multilepton_Lepton2_Id_LHCO = 2 ;
  if(abs(multilepton_Lepton3_Id)==11) multilepton_Lepton3_Id_LHCO = 1 ;
  else if(abs(multilepton_Lepton3_Id)==13) multilepton_Lepton3_Id_LHCO = 2 ;
 
  //
  // LHCO phi convention
  //
  
  float multilepton_Lepton1_phi       = Phi_0_2Pi(multilepton_Lepton1_P4.Phi());
  float multilepton_Lepton2_phi       = Phi_0_2Pi(multilepton_Lepton1_P4.Phi());
  float multilepton_Lepton3_phi       = Phi_0_2Pi(multilepton_Lepton1_P4.Phi());
  float multilepton_Bjet1_phi         = Phi_0_2Pi(multilepton_Bjet1_P4.Phi());
  float multilepton_Bjet2_phi         = Phi_0_2Pi(multilepton_Bjet2_P4.Phi());	 
  float multilepton_JetHighestPt1_phi = Phi_0_2Pi(multilepton_JetHighestPt1_P4.Phi());
  float multilepton_JetHighestPt2_phi = Phi_0_2Pi(multilepton_JetHighestPt1_P4.Phi());
  float multilepton_JetClosestMw1_phi = Phi_0_2Pi(multilepton_JetClosestMw1_P4.Phi());
  float multilepton_JetClosestMw2_phi = Phi_0_2Pi(multilepton_JetClosestMw2_P4.Phi());
  float multilepton_JetLowestMjj1_phi = Phi_0_2Pi(multilepton_JetLowestMjj1_P4.Phi());
  float multilepton_JetLowestMjj2_phi = Phi_0_2Pi(multilepton_JetLowestMjj2_P4.Phi());
  
	
  // l1
  std::string l1_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d", nobj,multilepton_Lepton1_Id_LHCO,multilepton_Lepton1_P4.Eta(),multilepton_Lepton1_phi,multilepton_Lepton1_P4.Pt(),0.0,multilepton_Lepton1_Id/abs(multilepton_Lepton1_Id),0,0,0,0));
  nobj++; 

  // l2
  std::string l2_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d", nobj,multilepton_Lepton2_Id_LHCO,multilepton_Lepton2_P4.Eta(),multilepton_Lepton2_phi,multilepton_Lepton2_P4.Pt(),0.0,multilepton_Lepton2_Id/abs(multilepton_Lepton2_Id),0,0,0,0));
  nobj++;	

  // l3
  std::string l3_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d", nobj,multilepton_Lepton3_Id_LHCO,multilepton_Lepton3_P4.Eta(),multilepton_Lepton3_phi,multilepton_Lepton3_P4.Pt(),0.0,multilepton_Lepton3_Id/abs(multilepton_Lepton3_Id),0,0,0,0));
  nobj++;
 				    
  //										    
  std::string j1_fline;
  std::string j2_fline;
  
  if ( _processLHCO_RECO == 5 || _processLHCO_RECO == 6 || _processLHCO_RECO == 4 || _processLHCO_RECO == 3 )
  {			     
    // j1
    j1_fline = std::string(Form("%d	  %d	 %.2f	    %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
     			      nobj,4,multilepton_JetClosestMw1_P4.Eta(),multilepton_JetClosestMw1_phi,multilepton_JetClosestMw1_P4.Pt(),0.0,1,0,0,0,0));
    nobj++;
  
    // j2
    j2_fline = std::string(Form("%d	  %d	  %.2f      %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
     			      nobj,4,multilepton_JetClosestMw2_P4.Eta(),multilepton_JetClosestMw2_phi,multilepton_JetClosestMw2_P4.Pt(),0.0,1,0,0,0,0));
    nobj++;
  }
  else if ( _processLHCO_RECO == 1 || _processLHCO_RECO == 2 ) 
  {
    // j1
    j1_fline = std::string(Form("%d	  %d	  %.2f      %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
     			      nobj,4,multilepton_JetLowestMjj1_P4.Eta(),multilepton_JetLowestMjj1_phi,multilepton_JetLowestMjj1_P4.Pt(),0.0,1,0,0,0,0));
    nobj++;
  
    // j2
    j2_fline = std::string(Form("%d	  %d	   %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d", 
     			      nobj,4,multilepton_JetLowestMjj2_P4.Eta(),multilepton_JetLowestMjj2_phi,multilepton_JetLowestMjj2_P4.Pt(),0.0,1,0,0,0,0));
    nobj++;
  }
  
  // for B-jet mass
  TLorentzVector BJet1;
  BJet1.SetPtEtaPhiE(multilepton_Bjet1_P4.Pt(), multilepton_Bjet1_P4.Eta(), multilepton_Bjet1_P4.Phi(), multilepton_Bjet1_P4.E());
 
  TLorentzVector BJet2;
  BJet2.SetPtEtaPhiE(multilepton_Bjet2_P4.Pt(), multilepton_Bjet2_P4.Eta(), multilepton_Bjet2_P4.Phi(), multilepton_Bjet2_P4.E());
         

  // bj1
  std::string b1_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d",
 					nobj,4,multilepton_Bjet1_P4.Eta(),multilepton_Bjet1_phi,multilepton_Bjet1_P4.Pt(),BJet1.M(),1,2,0,0,0));
  nobj++;

 				     
  // bj2 
  std::string b2_fline = std::string(Form("%d	  %d	 %.2f	  %.2f     %.2f     %.2f     %d     %d     %d	  %d	 %d",
 					nobj,4,multilepton_Bjet2_P4.Eta(),multilepton_Bjet2_phi,multilepton_Bjet2_P4.Pt(),BJet2.M(),1,2,0,0,0));
  nobj++;

  // met
  std::string met_fline = std::string(Form("%d     %d	  %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d	  %d",
 					 nobj,6,0.,Phi_0_2Pi(multilepton_mET.Phi()),multilepton_mET.Pt(),0.,0,0,0,0,0));
  nobj++;
  

 //    
 fout_RECO << fline00	<< std::endl;
 fout_RECO << fline0	<< std::endl;
 fout_RECO << l1_fline  << std::endl;
 fout_RECO << l2_fline  << std::endl;
 fout_RECO << l3_fline  << std::endl;	
 if (!(_processLHCO_RECO == 7 || _processLHCO_RECO == 8)) fout_RECO << j1_fline  << std::endl;// don't print jets for ttW hypothesis
 if (!(_processLHCO_RECO == 7 || _processLHCO_RECO == 8)) fout_RECO << j2_fline  << std::endl;// don't print jets for ttW hypothesis
 fout_RECO << b1_fline  << std::endl;
 fout_RECO << b2_fline  << std::endl;
 fout_RECO << met_fline << std::endl;      
       
}


float TTbarHiggsMultileptonAnalysis::Phi_0_2Pi(float phi)
{
 float phi_0_2pi = phi;
 if (phi>= TMath::Pi()) phi_0_2pi  -= 2.*TMath::Pi();
 if (phi<0.)            phi_0_2pi  += 2.*TMath::Pi();
 return phi_0_2pi;
}
