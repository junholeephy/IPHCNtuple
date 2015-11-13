#include "../include/TTbarHiggsMultileptonAnalysis.h"

TTbarHiggsMultileptonAnalysis::TTbarHiggsMultileptonAnalysis() {}

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

        if( nLep )
        {	 
            // PrintEventList(vSelectedLeptons,vSelectedJets);
            // theHistoManager->fillHisto("CutFlow", "noSel", "emu", sampleName.Data(),  1, theweight);

            if( (vSelectedBTagJets.size() + vSelectedNonBTagJets.size()) >= 4)
            {
                theHistoManager->fillHisto("CutFlow", "noSel", "emu", sampleName.Data(),  2, theweight);

                if(vSelectedBTagJets.size() >=1)
                {		     
                    theHistoManager->fillHisto("CutFlow", "noSel", "emu", sampleName.Data(),  3, theweight);
                }//NBjet selection
            }//end NJet selection
        }//end nlepton selection

        ThreeLeptonSelection(vSelectedLeptons, vSelectedJets, vSelectedBTagJets, vSelectedNonBTagJets, nLooseBtag, nMediumBtag);

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

    bool nLep        = ( vSelectedLeptons.size()                                  == 3 );
    bool nJets       = ( (vSelectedBTagJets.size() + vSelectedNonBTagJets.size()) >= 4 );
    bool nLooseBtag  = ( nLooseBJets                                              >= 2 );
    bool nMediumBtag = ( nMediumBJets                                             >= 1 );

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


    }
}
