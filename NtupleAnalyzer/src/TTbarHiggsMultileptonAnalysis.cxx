#include "../include/TTbarHiggsMultileptonAnalysis.h"

TTbarHiggsMultileptonAnalysis::TTbarHiggsMultileptonAnalysis() 
{
   _printLHCO = false;
   _processLHCO = -1;
}

void TTbarHiggsMultileptonAnalysis::InitLHCO(int process)
{
   _printLHCO = true;
   _processLHCO = process;

   fout.open("LHCO.txt");

   nLHCOevts = 0 ;
   fline00 = "#   typ     eta    phi       pt  jmass  ntrk  btag   had/em  dummy dummy";
   del = "    ";
   trig = "8";
}

TTbarHiggsMultileptonAnalysis::TTbarHiggsMultileptonAnalysis(TString inputfilename, TTree *tree, TString theSampleName) 
{    
    if (tree == 0) 
    {
        TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(inputfilename.Data());
        if (!f || !f->IsOpen()) 
        {
            f = new TFile(inputfilename.Data());
        }
        f->GetObject("Nt",tree);
    }
    Init(tree);

    theHistoManager = new HistoManager();

    sampleName = theSampleName;

    std::string foutlog = "output.txt";
    fevc = fopen(foutlog.c_str(),"w");
}

void TTbarHiggsMultileptonAnalysis::createHistograms()
{    
    outputfile = new TFile("output.root", "recreate");

    theHistoManager->addHisto("CutFlow",      "noSel", "emu", sampleName.Data(), 10, 0, 10);

    theHistoManager->addHisto("MuonPt",      "noSel", "emu", sampleName.Data(), 100, 0, 200);
    theHistoManager->addHisto("MuonEta",     "noSel", "emu", sampleName.Data(), 100, -3, 3);

    theHistoManager->addHisto("ElectronPt",  "noSel", "emu", sampleName.Data(), 100, 0, 200);
    theHistoManager->addHisto("ElectronEta", "noSel", "emu", sampleName.Data(), 100, -3, 3);

    theHistoManager->addHisto("JetPt",       "noSel", "emu", sampleName.Data(), 100, 0, 200);
    theHistoManager->addHisto("JetEta",      "noSel", "emu", sampleName.Data(), 100, -3, 3);
}

void TTbarHiggsMultileptonAnalysis::writeHistograms()
{  
    outputfile->cd();

    std::vector<TH1F*> the1DHisto =  theHistoManager->getHisto1D_list();
    std::vector<TH2F*> the2DHisto =  theHistoManager->getHisto2D_list();

    for(unsigned int i=0; i<the1DHisto.size(); i++)  the1DHisto[i]->Write();
    for(unsigned int i=0; i<the2DHisto.size(); i++)  the2DHisto[i]->Write();         
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

        if(jentry > 1000000) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        //int nPart = vTruth->at(0).mc_truth_n();
        //std::vector<int> labTruth = vTruth->at(0).mc_truth_label();

        //std::cout << "Jusqu'ici tout va bien 1" << std::endl;

        int mc_truth_h0W1_id = 0;
        int mc_truth_h0W2_id = 0;
        int mc_truth_h0Z1_id = 0;
        int mc_truth_h0Z2_id = 0;
        int mc_truth_h0tau1_id = 0;
        int mc_truth_h0tau2_id = 0;

        /*
        
        for(int i=0;i<nPart;i++)
        {	     
            if( labTruth[i] == 12 ) // h0W1
            {
                mc_truth_h0W1_id = vTruth->at(0).mc_truth_id()[i];
            }	     
            else if( labTruth[i] == 13 ) // h0W2
            {
                mc_truth_h0W2_id = vTruth->at(0).mc_truth_id()[i];
            }	     
            else if( labTruth[i] == 14 ) // h0Z1
            {
                mc_truth_h0Z1_id = vTruth->at(0).mc_truth_id()[i];
            }	     
            else if( labTruth[i] == 15 ) // h0Z2
            {
                mc_truth_h0Z2_id = vTruth->at(0).mc_truth_id()[i];
            }	     
            else if( labTruth[i] == 16 ) // h0tau1
            {
                mc_truth_h0tau1_id = vTruth->at(0).mc_truth_id()[i];
            }	     
            else if( labTruth[i] == 17 ) // h0tau2
            {
                mc_truth_h0tau2_id = vTruth->at(0).mc_truth_id()[i];
            }	     
        }

        std::cout << "Jusqu'ici tout va bien 1" << std::endl;

        // select only ttH multilepton events at truth level
        bool isHtoWW = (abs(mc_truth_h0W1_id) == 24 &&
                abs(mc_truth_h0W2_id) == 24);
        bool isHtoZZ = (abs(mc_truth_h0Z1_id) == 23 &&
                abs(mc_truth_h0Z2_id) == 23);
        bool isHtoTT = (abs(mc_truth_h0tau1_id) == 15 &&
                abs(mc_truth_h0tau2_id) == 15);
        if( !(isHtoWW || isHtoZZ || isHtoTT) ) continue;

        */

        float theweight = vEvent->at(0).mc_weight();

        std::vector<Muon>     vSelectedMuons;
        std::vector<Electron> vSelectedElectrons;
        std::vector<Lepton>   vSelectedLeptons;
        std::vector<Jet>      vSelectedNonBTagJets;
        std::vector<Jet>      vSelectedBTagJets;
        std::vector<Jet>      vSelectedJets;

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

        ThreeLeptonSelection(vSelectedLeptons, vSelectedJets, vSelectedBTagJets, vSelectedNonBTagJets, nLooseBJets, nMediumBJets);

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

      if (_printLHCO) PrintLHCOforMadweight();

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

void TTbarHiggsMultileptonAnalysis::PrintLHCOforMadweight()
{

  /*if( proc<-1 || proc > 6 )
  {
      std::cout << "proc can only take following values: -1,1,2,3,4,5,6" << std::endl;
      std::cout << "3l final state specific" << std::endl;    
      std::cout << "1,2,3,4: specific to the ttH final state, cf patches provided to madweight" << std::endl;
      std::cout << "5: specific to the ttZ with l+l-l+ final state" << std::endl;
      std::cout << "6: specific to the ttZ with l+l-l- final state" << std::endl;
      std::cout << "-1: no selection on the final state applied" << std::endl;

      return;
   }*/


  //
  //std::string fline00 = "#   typ     eta    phi       pt  jmass  ntrk  btag   had/em  dummy dummy";
  //std::string del = "   ";
  //std::string trig = "8";
  //std::cout <<"proc " << proc << std::endl;

  //std::cout << "printLHCO 0" << std::endl;

  //if (vTruth->at(0).ttbar_decay()==-1 || vTruth->at(0).boson_decay()==-1) return;

/*  if (!((vTruth->at(0).ttbar_decay()==1 && vTruth->at(0).boson_decay()==0) || (vTruth->at(0).ttbar_decay()==2 && vTruth->at(0).boson_decay()==1))) return;
  if (vTruth->at(0).JetsHighestPt_phi().size()<2 || vTruth->at(0).JetsClosestMw_phi().size()<2 || vTruth->at(0).JetsLowestMjj_phi().size()<2) return;

  fline0 = "0      " + std::string(Form("%d",nLHCOevts)) + "     " + trig;

  int nobj = 1;

  int multilepton_Lepton1_Id_LHCO = -666;
  int multilepton_Lepton2_Id_LHCO = -666;
  int multilepton_Lepton3_Id_LHCO = -666;

  //
  // LHCO lepton ID convention
  //

  //std::cout << vTruth->at(0).Leptons_id().size()<<std::endl;

  //if (vTruth->at(0).ttbar_decay()==-1 || vTruth->at(0).boson_decay()==-1) return;

  if (!((vTruth->at(0).ttbar_decay()==1 && vTruth->at(0).boson_decay()==0) || (vTruth->at(0).ttbar_decay()==2 && vTruth->at(0).boson_decay()==1))) return;
  if (vTruth->at(0).JetsHighestPt_phi().size()<2 || vTruth->at(0).JetsClosestMw_phi().size()<2 || vTruth->at(0).JetsLowestMjj_phi().size()<2) return;

  fline0 = "0      " + std::string(Form("%d",nLHCOevts)) + "     " + trig;

  int nobj = 1;

  int multilepton_Lepton1_Id_LHCO = -666;
  int multilepton_Lepton2_Id_LHCO = -666;
  int multilepton_Lepton3_Id_LHCO = -666;

  //
  // LHCO lepton ID convention
  //

  //std::cout << vTruth->at(0).Leptons_id().size()<<std::endl;

  if(abs(vTruth->at(0).Leptons_id().at(0))==11) multilepton_Lepton1_Id_LHCO = 1 ;
  else if(abs(vTruth->at(0).Leptons_id().at(0))==13) multilepton_Lepton1_Id_LHCO = 2 ;
  if(abs(vTruth->at(0).Leptons_id().at(1))==11) multilepton_Lepton2_Id_LHCO = 1 ;
  else if(abs(vTruth->at(0).Leptons_id().at(1))==13) multilepton_Lepton2_Id_LHCO = 2 ;
  if(abs(vTruth->at(0).Leptons_id().at(2))==11) multilepton_Lepton3_Id_LHCO = 1 ;
  else if(abs(vTruth->at(0).Leptons_id().at(2))==13) multilepton_Lepton3_Id_LHCO = 2 ;


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

  // l1
  std::string l1_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
                    nobj,multilepton_Lepton1_Id_LHCO,vTruth->at(0).Leptons_eta().at(0),multilepton_Lepton1_phi,vTruth->at(0).Leptons_pt().at(0),0.0,vTruth->at(0).Leptons_id().at(0)/         abs(vTruth->at(0).Leptons_id().at(0)),0,0,0,0));
  nobj++;

  // l2
  std::string l2_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
                    nobj,multilepton_Lepton2_Id_LHCO,vTruth->at(0).Leptons_eta().at(1),multilepton_Lepton2_phi,vTruth->at(0).Leptons_pt().at(1),0.0,vTruth->at(0).Leptons_id().at(1)/         abs(vTruth->at(0).Leptons_id().at(1)),0,0,0,0));
  nobj++;

  // l3
  std::string l3_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
                    nobj,multilepton_Lepton3_Id_LHCO,vTruth->at(0).Leptons_eta().at(2),multilepton_Lepton3_phi,vTruth->at(0).Leptons_pt().at(2),0.0,vTruth->at(0).Leptons_id().at(2)/         abs(vTruth->at(0).Leptons_id().at(2)),0,0,0,0));
  nobj++;

  //                                            
  std::string j1_fline;
  std::string j2_fline;

  if ( _processLHCO == 5 || _processLHCO == 6 || _processLHCO == 4 || _processLHCO == 3 )
  {
    // j1
    j1_fline = std::string(Form("%d   %d     %.2f       %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
                  nobj,4,vTruth->at(0).JetsClosestMw_eta().at(0),multilepton_JetClosestMw1_phi,vTruth->at(0).JetsClosestMw_pt().at(0),0.0,1,0,0,0,0));
    nobj++;

    // j2
    j2_fline = std::string(Form("%d   %d      %.2f      %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
                  nobj,4,vTruth->at(0).JetsClosestMw_eta().at(1),multilepton_JetClosestMw2_phi,vTruth->at(0).JetsClosestMw_pt().at(1),0.0,1,0,0,0,0));
    nobj++;
  }
  else
  {
    // j1
    j1_fline = std::string(Form("%d   %d      %.2f      %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
                  nobj,4,vTruth->at(0).JetsLowestMjj_eta().at(0),multilepton_JetLowestMjj1_phi,vTruth->at(0).JetsLowestMjj_pt().at(0),0.0,1,0,0,0,0));
    nobj++;

    // j2
    j2_fline = std::string(Form("%d   %d       %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
                  nobj,4,vTruth->at(0).JetsLowestMjj_eta().at(1),multilepton_JetLowestMjj2_phi,vTruth->at(0).JetsLowestMjj_pt().at(1),0.0,1,0,0,0,0));
    nobj++;
  }
  //
  TLorentzVector BJet1;
  BJet1.SetPtEtaPhiE(vTruth->at(0).Bjets_pt().at(0), vTruth->at(0).Bjets_eta().at(0), vTruth->at(0).Bjets_phi().at(0), vTruth->at(0).Bjets_E().at(0));

  TLorentzVector BJet2;
  BJet2.SetPtEtaPhiE(vTruth->at(0).Bjets_pt().at(0), vTruth->at(0).Bjets_eta().at(0), vTruth->at(0).Bjets_phi().at(0), vTruth->at(0).Bjets_E().at(0));


  // bj1
  std::string b1_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
                    nobj,4,vTruth->at(0).Bjets_eta().at(0),multilepton_Bjet1_phi,vTruth->at(0).Bjets_pt().at(0),BJet1.M(),1,2,0,0,0));
  nobj++;


  // bj2 
  std::string b2_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
                    nobj,4,vTruth->at(0).Bjets_eta().at(1),multilepton_Bjet2_phi,vTruth->at(0).Bjets_pt().at(1),BJet2.M(),1,2,0,0,0));
  nobj++;

  // met
  std::string met_fline = std::string(Form("%d     %d     %.2f     %.2f     %.2f     %.2f     %d     %d     %d     %d     %d",
                     nobj,6,0.,Phi_0_2Pi(vTruth->at(0).metGen_phi()),vTruth->at(0).metGen_pt(),0.,0,0,0,0,0));
  nobj++;


 //    
 fout << fline00   << std::endl;
 fout << fline0    << std::endl;
 fout << l1_fline  << std::endl;
 fout << l2_fline  << std::endl;
 fout << l3_fline  << std::endl;
 fout << j1_fline  << std::endl;
 fout << j2_fline  << std::endl;
 fout << b1_fline  << std::endl;
 fout << b2_fline  << std::endl;
 fout << met_fline << std::endl;

 nLHCOevts++; */

}

float TTbarHiggsMultileptonAnalysis::Phi_0_2Pi(float phi)
{
 float phi_0_2pi = phi;
 if (phi>= TMath::Pi()) phi_0_2pi  -= 2.*TMath::Pi();
 if (phi<0.)            phi_0_2pi  += 2.*TMath::Pi();
 return phi_0_2pi;
}

void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection(std::vector<Lepton> vSelectedLeptons, 
                                                         std::vector<Jet>    vSelectedJets, 
                                                         std::vector<Jet>    vSelectedBTagJets, 
                                                         std::vector<Jet>    vSelectedNonBTagJets,
                                                         int                 nLooseBJets, 
                                                         int                 nMediumBJets)
{
    std::cout << "Number of selected leptons:    " << vSelectedLeptons.size()       << std::endl;
    std::cout << "Number of selected jets   :    " << (vSelectedBTagJets.size() 
                                                     + vSelectedNonBTagJets.size()) << std::endl;
    std::cout << "Number of loose b jets    :    " << nLooseBJets                   << std::endl;
    std::cout << "Number of medium b jets   :    " << nMediumBJets                  << std::endl;

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

    std::cout << "nLep: "              << nLep 
              << " nJets: "            << nJets 
              << " pass_OSSF: "        << pass_OSSF 
              << " leading_lep_pt: "   << leading_lep_pt 
              << " following_lep_pt: " << following_lep_pt
              << " passMll12Gt12: "    << passMll12Gt12 
              << " nLooseBtag: "       << nLooseBtag 
              << " nMediumBtag: "      << nMediumBtag << std::endl; 

    if ( nLep && nJets && pass_OSSF && leading_lep_pt && following_lep_pt && passMll12Gt12 && (nLooseBtag || nMediumBtag) )
    {

        // #############################
        // # FILLING MULTILEPTON CLASS #
        // #############################

        std::cout << "******** SELECTED EVENT OF THE THREE LEPTON CATEGORY! ********" << std::endl;
        std::cout << "Number of selected leptons:    " << vSelectedLeptons.size()     << std::endl;
        std::cout << "Number of selected jets   :    " << (vSelectedBTagJets.size() + vSelectedNonBTagJets.size()) << std::endl;
        std::cout << "Number of loose b jets    :    " << nLooseBJets                 << std::endl;
        std::cout << "Number of medium b jets   :    " << nMediumBJets                << std::endl;

    }
}
