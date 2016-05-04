#include "../include/TTbarHiggsChargeFlip.h"
#include "TSystem.h"
#include "Helper.cxx"


TTbarHiggsChargeFlip::TTbarHiggsChargeFlip() 
{
}

TTbarHiggsChargeFlip::TTbarHiggsChargeFlip(TString inputFileName, TChain *tree, TString sampleName, TString treeName, TString outputFileName, bool isdata, float xsec, float lumi, int nowe, int nmax)
{    

    //
    _isdata = isdata;
    _xsec = xsec;
    _lumi = lumi;
    _nowe = nowe;
    _nmax = nmax;
    _outputFileName = outputFileName;
    _sampleName = sampleName;

    //_file_PVreweighting = TFile::Open("/home-pbs/lebihan/someone/ttH_070116/ttH/NtupleAnalyzer/test/PUweight.root");
    //_h_PV = (TH1F*)_file_PVreweighting->Get("PU_reweighting");

    //
    tree = new TChain(treeName.Data());

    std::ifstream infile;
    infile.open(inputFileName.Data());
    std::string ifile = "";
    while( getline(infile, ifile) )
    {
        std::string fnameStr = std::string(ifile);
        tree->Add(fnameStr.c_str());
    }   
    infile.close();
    Init(tree);

    theHistoManager = new HistoManager();

    TString outputfileNameRoot = _outputFileName+".root";
    outputfile = new TFile(outputfileNameRoot.Data(), "recreate");  

}

void TTbarHiggsChargeFlip::createHistograms()
{    
    outputfile->cd();
    initializeOutputTree();

    // General
    theHistoManager->addHisto("CutFlow",                                     "noSel",        "",   "",  10,   0,     10);

    // Electron charge misassignement probabilities

    theHistoManager->addHisto("DiElectronInvariantMassSameSignTCA",          "noSel",        "",   "",  60,   60,   120);
    theHistoManager->addHisto("DiElectronInvariantMassSameSign",             "noSel",        "",   "",  60,   60,   120);

    //theHistoManager->addHisto("ProbabilityVspTBarrel",                       "noSel",        "",   "",  60,    0,  1000);
    //theHistoManager->addHisto("ProbabilityVspTEndcap",                       "noSel",        "",   "",  60,    0,  1000);

    // Background estimation CR 1 DY

    theHistoManager->addHisto("InvariantMassDilepton",                      "CR1_DY",        "",   "",  60,   60,   120);
    theHistoManager->addHisto("NJets",                                      "CR1_DY",        "",   "",  10, -0.5,   9.5);
    theHistoManager->addHisto("MET",                                        "CR1_DY",        "",   "",  30,    0,   200);
    theHistoManager->addHisto("LeadingLeptonpT",                            "CR1_DY",        "",   "",  45,   20,   200);
    theHistoManager->addHisto("SubLeadingLeptonpT",                         "CR1_DY",        "",   "",  45,   10,   100);

    // Background estimation CR 2 tt

    theHistoManager->addHisto("InvariantMassDilepton",                      "CR1_TT",        "",   "",  60,   60,   120);
    theHistoManager->addHisto("NJets",                                      "CR1_TT",        "",   "",  10, -0.5,   9.5);
    theHistoManager->addHisto("MET",                                        "CR1_TT",        "",   "",  30,    0,   200);
    theHistoManager->addHisto("LeadingLeptonpT",                            "CR1_TT",        "",   "",  45,   20,   200);
    theHistoManager->addHisto("SubLeadingLeptonpT",                         "CR1_TT",        "",   "",  45,   10,   100);
}


void TTbarHiggsMultileptonAnalysis::writeHistograms()
{  
    outputfile->cd();

    std::vector<TH1F*> the1DHisto =  theHistoManager->getHisto1D_list();
    std::vector<TH2F*> the2DHisto =  theHistoManager->getHisto2D_list();

    for(unsigned int i=0; i<the1DHisto.size(); i++)  the1DHisto[i]->Write();
    for(unsigned int i=0; i<the2DHisto.size(); i++)  the2DHisto[i]->Write();

    tOutput->Write();
    outputfile->Close();
}


void TTbarHiggsMultileptonAnalysis::Init(TChain *tree)
{
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;

    vEvent    = new std::vector<Event>();
    vElectron = new std::vector<Electron>();
    vMuon     = new std::vector<Muon>();
    vTau      = new std::vector<Tau>();
    vJet      = new std::vector<Jet>();
    vTruth    = new std::vector<Truth>();

    fChain->SetBranchAddress("Event",    &vEvent   );
    fChain->SetBranchAddress("Electron", &vElectron);
    fChain->SetBranchAddress("Muon",     &vMuon    );
    fChain->SetBranchAddress("Tau",      &vTau     );
    fChain->SetBranchAddress("Jet",      &vJet     );
    fChain->SetBranchAddress("Truth",    &vTruth   );

    Load_MVA();
}


void TTbarHiggsMultileptonAnalysis::Loop()
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntries();
    int nentries_max = nentries;
    if ( _nmax != -1 && _nmax < nentries ) nentries_max = _nmax;

    std::cout << "Number of input events = " << nentries << std::endl;
    std::cout << "Number of processed events = " << nentries_max << std::endl;

    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<nentries_max;jentry++) 
    {

        Long64_t ientry = fChain->LoadTree(jentry);
        if (ientry < 0) break;

        if(jentry%10000 == 0) std::cout << "number of processed events " << jentry << std::endl;
        //std::cout << "number of processed events " << jentry << std::endl;

        //if(jentry > 100000) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        //
        int pvn = vEvent->at(0).pv_n();
        theHistoManager->fillHisto("NumberOfPrimaryVertex", "noSel", "", "",  pvn, 1);


        if ( !_isdata )
        {
            weight = _lumi*_xsec/_nowe;
            mc_weight = vEvent->at(0).mc_weight();
            //weight_PV = _h_PV->GetBinContent(pvn);
            weight = weight * mc_weight; //*weight_PV;
        }
        else 
        {
            weight    = 1.;
            mc_weight = 1.;
            weight_PV = 1.; 

            /////////////////////////////////////
            // remove double counting in data

            // MuonEG
            // HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v
            // HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v
            // HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v
            // HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v
            // 
            // DoubleMuon
            // HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v
            // HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v
            // HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v
            // HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v
            // HLT_TripleMu_12_10_5_v
            // 
            // DoubleEG
            // HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
            // HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v
            // HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v
            // 
            // SingleMu
            // HLT_IsoMu20_v
            // HLT_IsoTkMu20_v
            // 
            // SingleEG
            // HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v
            // HLT_Ele23_WPLoose_Gsf_v

            // détricotage
            bool TRIGm = false, TRIGe = false, TRIGeData = false, TRIGmTk = false; 
            bool TRIGee = false, TRIGmm = false, TRIGme = false, TRIGem = false, TRIGmmTk = false;
            bool TRIGeee = false, TRIGmme = false, TRIGeem = false, TRIGmmm = false;

            int a = ( vEvent->at(0).ev_trigger_pass_byname_1() )%10;
            int b = ((vEvent->at(0).ev_trigger_pass_byname_1() -a)/10)%10;
            int c = ((vEvent->at(0).ev_trigger_pass_byname_1() -a-10*b)/100)%10;
            int d = ((vEvent->at(0).ev_trigger_pass_byname_1() -a-10*b-100*c)/1000)%10;
            int e =   vEvent->at(0).ev_trigger_pass_byname_1();

            if (a==1 || a==3 || a==6 || a==8) TRIGeee  = true; // HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*
            if (a==2 || a==3 || a==7 || a==8) TRIGme   = true; // HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v* 
            if (a>=5 )                        TRIGem   = true; // HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v* 
            if (b==1 || b==3 || b==6 || b==8) TRIGee   = true; // HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v* 
            if (b==2 || b==3 || b==7 || b==8) TRIGmm   = true; // HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v* 
            if (b>=5 )                        TRIGmmTk = true; // HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v* 
            if (c==1 || c==3 || c==6 || c==8) TRIGm    = true; // HLT_IsoMu20_v* 
            if (c==2 || c==3 || c==7 || c==8) TRIGmTk  = true; // HLT_IsoTkMu20_v* 
            if (c>=5 )                        TRIGe    = true; // HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v* 
            if (d==1 || d==3 || d==6 || d==8) TRIGeData= true; // HLT_Ele23_WPLoose_Gsf_v* 
            if (d==2 || d==3 || d==7 || d==8) TRIGmme  = true; // HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v*  but doesn't exist ?
            if (d>=5 )                        TRIGeem  = true; // HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v* but doesn't exist ?
            if (e>=10000)                     TRIGmmm  = true; // HLT_TripleMu_12_10_5_v*             but doesn't exist ?

            bool E = false, M = false, EE = false, MM = false, EM = false;
            if ( TRIGme || TRIGem || TRIGeem || TRIGmme ) EM = true;
            if ( TRIGmm || TRIGmmTk || TRIGmmm )          MM = true;
            if ( TRIGee || TRIGeee )	              EE = true;
            if ( TRIGm  || TRIGmTk )                      M  = true;
            if ( TRIGe  || TRIGeData )                    E  = true;

            // new code from Xavier (with next Ntuple production)
            //         EM = ( vEvent->at(0).ev_pass_eem() || vEvent->at(0).ev_pass_em()     || vEvent->at(0).ev_pass_mme()    || vEvent->at(0).ev_pass_me() );
            //         MM = ( vEvent->at(0).ev_pass_mm()  || vEvent->at(0).ev_pass_mmTk()   || vEvent->at(0).ev_pass_mmnoDz() || vEvent->at(0).ev_pass_mmTknoDz() || vEvent->at(0).ev_pass_mmm() ) ;
            //         EE = ( vEvent->at(0).ev_pass_ee()  || vEvent->at(0).ev_pass_eenoDz() || vEvent->at(0).ev_pass_eee() );
            //         M  = ( vEvent->at(0).ev_pass_m()   || vEvent->at(0).ev_pass_mTk() );
            //         E  = ( vEvent->at(0).ev_pass_e()   || vEvent->at(0).ev_pass_eData() );

            // old code from Xavier
            //         int tab[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            //         int a = vEvent->at(0).ev_trigger_pass_byname_1();
            //         int n = 0, size = 0;
            // 
            //         do{
            //             n = a%10;
            //             a = a/10;
            //             tab[size] = n;
            //             size = size + 1;
            //         }while(a!=0);
            // 
            //         bool E = false, M = false, EE = false, MM = false, EM = false;
            //         int result_trigger = 0;
            //         if (size > 2) {if ( tab[3] == 1                ) E  = true;}
            //         if (size > 1) {if ( tab[2] == 1 || tab[2] == 2 ) M  = true;}
            //         if (size > 0) {if ( tab[1] == 1                ) EE = true;}
            //         if (size > 1) {if ( tab[1] == 2 || tab[1] == 5 ) MM = true;}
            //         if ( tab[0] == 2 || tab[0] == 5 )                EM = true;

            bool emdataset = _sampleName.Contains("MuonEG");
            bool mmdataset = _sampleName.Contains("DoubleM");
            bool eedataset = _sampleName.Contains("DoubleE");
            bool mdataset  = _sampleName.Contains("SingleM");
            bool edataset  = _sampleName.Contains("SingleE");

            int result_trigger = 0;
            if ( EM  &&                               (emdataset) ) result_trigger = 1;
            if ( !EM && MM  &&                        (mmdataset) ) result_trigger = 1;
            if ( !EM && !MM && EE  &&                 (eedataset) ) result_trigger = 1;
            if ( !EM && !MM && !EE && M  &&           (mdataset ) ) result_trigger = 1;
            if ( !EM && !MM && !EE && !M && E &&      (edataset ) ) result_trigger = 1;

            if(result_trigger == 1)
            {
                weight = 1;
                //std::cout << " EM: " << EM 
                //          << " MM: " << MM 
                //          << " EE: " << EE
                //          << " M:  " << M
                //          << " E:  " << E                    << std::endl;
                //std::cout << " sampleName: " << _sampleName   << std::endl;
                //std::cout << " weight:     " << weight       << std::endl; 
            }
            else
            {
                weight = 0;
                //std::cout << " EM: " << EM 
                //          << " MM: " << MM  
                //          << " EE: " << EE
                //          << " M:  " << M
                //          << " E:  " << E                    << std::endl;
                //std::cout << " sampleName: " << _sampleName   << std::endl;
                //std::cout << " weight:     " << weight       << std::endl;  
            }
        }


        //---------------------------
        // initialisation
        //---------------------------
        vLeptons.clear();
        vSelectedMuons.clear();
        vSelectedElectrons.clear();
        vSelectedLeptons.clear();
        vFakeMuons.clear();
        vFakeElectrons.clear();
        vFakeLeptons.clear();	
        vSelectedNonBTagJets.clear();
        vSelectedBTagJets.clear();
        vSelectedMediumBTagJets.clear();
        vSelectedJets.clear();

        is_2lss_TTH_SR    = false;
        is_2lss_JM_SB     = false;
        is_2lss_LepMVA_SB = false;
        is_emu_TT_CR      = false;

        is_3l_TTH_SR      = false;
        is_3l_WZ_CR       = false; 
        is_3l_WZrel_CR    = false;
        is_3l_TTZ_CR      = false;
        is_Zl_CR          = false;

        //---------------------------
        //trigger
        //---------------------------
        is_trigger = false;
        if ( vEvent->at(0).ev_trigger_pass_byname_1() >= 1 ) is_trigger = true;

        //---------------------------
        //muons
        //---------------------------
        for(unsigned int imuon=0; imuon < vMuon->size() ; imuon++)
        {   
            Lepton l; l.setLepton(&vMuon->at(imuon),imuon,0);

            if ( vMuon->at(imuon).isTightTTH() )
            {
                vSelectedMuons.push_back(vMuon->at(imuon));
                vSelectedLeptons.push_back(l);
            }
            else if ( vMuon->at(imuon).isFakeableTTH() )
            {
                vFakeMuons.push_back(vMuon->at(imuon));
                vFakeLeptons.push_back(l);
            }

            vLeptons.push_back(l);

            theHistoManager->fillHisto("MuonPt",                            "noSel",        "",   "",  vMuon->at(imuon).pt(),             weight);
            theHistoManager->fillHisto("MuonEta",                           "noSel",        "",   "",  vMuon->at(imuon).eta(),            weight);
            theHistoManager->fillHisto("MuonMVA",                           "noSel",        "",   "",  vMuon->at(imuon).lepMVA(),         weight);
        }     

        //---------------------------
        // electrons
        //---------------------------
        for(unsigned int ielectron=0; ielectron < vElectron->size() ; ielectron++)
        {   
            Lepton l; l.setLepton(&vElectron->at(ielectron),ielectron,1);

            if ( vElectron->at(ielectron).isTightTTH() )
            {
                vSelectedElectrons.push_back(vElectron->at(ielectron));	     
                vSelectedLeptons.push_back(l);
            }
            else if ( vElectron->at(ielectron).isFakeableTTH() )
            {
                vFakeElectrons.push_back(vElectron->at(ielectron));	     
                vFakeLeptons.push_back(l);
            }

            vLeptons.push_back(l);

            theHistoManager->fillHisto("ElectronPt",                        "noSel",        "",   "",  vElectron->at(ielectron).pt(),     weight);
            theHistoManager->fillHisto("ElectronEta",                       "noSel",        "",   "",  vElectron->at(ielectron).eta(),    weight);
            theHistoManager->fillHisto("ElectronMVA",                       "noSel",        "",   "",  vElectron->at(ielectron).lepMVA(), weight);
        }  

        //---------------------------
        // taus
        //---------------------------
        for(unsigned int itau=0; itau < vTau->size() ; itau++)
        {
            Lepton l; l.setLepton(&vTau->at(itau),itau,1);

            vSelectedTaus.push_back(vTau->at(itau));
            //vSelectedLeptons.push_back(l);

            //vLeptons.push_back(l);

            theHistoManager->fillHisto("TauPt",                             "noSel",        "",   "",  vTau->at(itau).pt(),               weight);
            theHistoManager->fillHisto("TauEta",                            "noSel",        "",   "",  vTau->at(itau).eta(),              weight);
        }

        std::sort(vLeptons.begin(), vLeptons.end(), SortingLeptonPt);
        std::sort(vSelectedLeptons.begin(), vSelectedLeptons.end(), SortingLeptonPt);
        std::sort(vFakeLeptons.begin(), vFakeLeptons.end(), SortingLeptonPt);

        //---------------------------
        // b-jets
        //---------------------------
        nLooseBJets  = 0;
        nMediumBJets = 0;

        for(unsigned int ijet=0; ijet < vJet->size() ; ijet++)
        {

            // version for 74x as in the AN...:
            if( vJet->at(ijet).CSVv2() > 0.423 ) nLooseBJets++;
            if( vJet->at(ijet).CSVv2() > 0.814 ) nMediumBJets++;

            if(vJet->at(ijet).CSVv2() >= 0.423 ) vSelectedBTagJets.push_back(vJet->at(ijet));
            else                                 vSelectedNonBTagJets.push_back(vJet->at(ijet));
            if(vJet->at(ijet).CSVv2() >= 0.814 ) vSelectedMediumBTagJets.push_back(vJet->at(ijet));

            // to be updated for 76x ???
            //             if( vJet->at(ijet).CSVv2() > 0.460 ) nLooseBJets++;
            //             if( vJet->at(ijet).CSVv2() > 0.800 ) nMediumBJets++;
            // 
            //             if(vJet->at(ijet).CSVv2() >= 0.460 ) vSelectedBTagJets.push_back(vJet->at(ijet));
            //             else                                 vSelectedNonBTagJets.push_back(vJet->at(ijet));
            //             if(vJet->at(ijet).CSVv2() >= 0.800 ) vSelectedMediumBTagJets.push_back(vJet->at(ijet));

            vSelectedJets.push_back(vJet->at(ijet));

            theHistoManager->fillHisto("JetPt",                             "noSel",        "",   "",  vJet->at(ijet).pt(),               weight);
            theHistoManager->fillHisto("JetEta",                            "noSel",        "",   "",  vJet->at(ijet).eta(),              weight);
            theHistoManager->fillHisto("JetCSVv2",                          "noSel",        "",   "",  vJet->at(ijet).CSVv2(),            weight);
        }

        theHistoManager->fillHisto("MET",                               "noSel",        "",   "",  vEvent->at(0).metpt(),            weight );

        theHistoManager->fillHisto("CutFlow",                        "noSel", "", "", 1, 1);

        if( vMuon->size()+vElectron->size() == 3 )     
        {
            theHistoManager->fillHisto("CutFlow",             "ThreePreselected", "", "", 1, 1);
        }


        //---------------------------
        //Selection for signal and control regions
        //---------------------------

        TwoLeptonsSameSignSelection_TTH2l(jentry);
        //TwoLeptonsSameSignSelection_LepMVA_sideband(jentry);
        //TwoLeptonsSameSignSelection_JetMultiplicity_sideband(jentry);
        DiLeptonSelection_TT_CR(jentry);

        ThreeLeptonSelection_TTH3l(jentry);
        ThreeLeptonSelection_CR_WZ(jentry);
        ThreeLeptonSelection_CR_WZrelaxed(jentry);
        ThreeLeptonSelection_CR_Zl(jentry);
        ThreeLeptonSelection_TTZ(jentry);

        //std::cout <<is_CR_TTl<<" "<< is_Zl_CR <<" " << is_CR_WZ<<" " << is_TTH3l<< std::endl;
        //if (is_TTH3l==true ) std::cout <<"is_TTH3l" << std::endl;
        if ( is_2lss_TTH_SR || is_3l_TTH_SR ) fillOutputTree();

        //---------------------------
        //Madweight LHCO stuff
        //---------------------------
        if ( !_isdata && _printLHCO_MC && ThreeLeptonSelection_TTH3l_MC()) PrintLHCOforMadweight_MC(jentry);

        // Common Selection:
        if ( !(vLeptons.size() >= 2
                    && vSelectedJets.size() >= 2) ) continue;

        float MET = vEvent->at(0).metpt();
        float METphi = vEvent->at(0).metphi();
        float METx = MET * TMath::Cos(METphi);
        float METy = MET * TMath::Sin(METphi);
        float METsum = vEvent->at(0).metsumet();
        int nlepsel = 0;
        float jet_px = 0, jet_py = 0, lep_px = 0, lep_py = 0, MHT = 0, met_ld = 0, jetht = 0;
        for (int i=0; i<vSelectedLeptons.size(); i++) {
            lep_px += vSelectedLeptons.at(i).p4().Px();
            lep_py += vSelectedLeptons.at(i).p4().Py();
            if ( vSelectedLeptons.at(i).pt() > 10. ) nlepsel++;
        }
        if ( nlepsel >= 2 && vSelectedLeptons.at(0).pt() < 20. ) nlepsel = -1.;
        TLorentzVector jetp4;
        for (int ijet=0; ijet < vSelectedJets.size(); ijet++) {
            jetp4.SetPtEtaPhiE(vSelectedJets.at(ijet).pt(), vSelectedJets.at(ijet).eta(), vSelectedJets.at(ijet).phi(), vSelectedJets.at(ijet).E());
            jet_px += jetp4.Px();
            jet_py += jetp4.Py();
            jetht += vSelectedJets.at(ijet).pt();
            theHistoManager->fillHisto("JetPt",  "Trig", "", "", vJet->at(ijet).pt(), weight);
        }
        MHT = sqrt( (jet_px+lep_px)*(jet_px+lep_px) + (jet_py+lep_py)*(jet_py+lep_py) );
        met_ld = 0.00397 * MET + 0.00265 * MHT;

        float Mllmin = 1000., Mllbest = 1000., Deltabest = 1000.;
        theHistoManager->fillHisto("nLep",   "PreSel", "", "", nlepsel, weight);
        theHistoManager->fillHisto("nLep loose", "PreSel", "", "", vLeptons.size(), weight);

        if ( nlepsel >= 2 ) {
            theHistoManager->fillHisto("lep1Pt", "PreSel", "", "", vSelectedLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep1Eta","PreSel", "", "", vSelectedLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lep2Pt", "PreSel", "", "", vSelectedLeptons.at(1).pt(), weight);
            theHistoManager->fillHisto("lep2Eta","PreSel", "", "", vSelectedLeptons.at(1).eta(), weight);
            if ( nlepsel >= 3 ) theHistoManager->fillHisto("lep3Pt", "PreSel", "", "", vSelectedLeptons.at(2).pt(), weight);
            if ( nlepsel >= 3 ) theHistoManager->fillHisto("lep3Eta","PreSel", "", "", vSelectedLeptons.at(2).eta(), weight);
            float lepq = vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge();
            if ( nlepsel >= 3 ) lepq += vSelectedLeptons.at(2).charge();
            theHistoManager->fillHisto("lepQ", "PreSel", "", "", lepq, weight);
            for (int i=0; i<vSelectedLeptons.size()-1; i++) {
                for (int j=i+1; j<vSelectedLeptons.size(); j++) {
                    float mll = ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M();
                    if ( mll < Mllmin ) Mllmin = mll;
                    if ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id()) {
                        if ( fabs(mll - 91.188) < Deltabest ) {
                            Mllbest = mll;
                            Deltabest = fabs(mll - 91.188);
                        }
                        theHistoManager->fillHisto("Mll", "PreSel", "", "", mll, weight);
                        if ( nlepsel >= 3 ) theHistoManager->fillHisto("Mll", "PreSel 3l", "", "", mll, weight);
                        if ( nLooseBJets>=2 || nMediumBJets>=1 ) theHistoManager->fillHisto("Mll", "PreSel btag", "", "", mll, weight);
                        if ( (nLooseBJets>=2 || nMediumBJets>=1 ) && nlepsel >= 3 ) theHistoManager->fillHisto("Mll", "PreSel btag 3l", "", "", mll, weight);
                    }
                }
            }
        }

        theHistoManager->fillHisto("nJets",    "PreSel", "", "", vSelectedJets.size(), weight);
        theHistoManager->fillHisto("nLooseB",  "PreSel", "", "", nLooseBJets, weight);
        theHistoManager->fillHisto("nMediumB", "PreSel", "", "", nMediumBJets, weight);
        for (int ijet=0; ijet < vSelectedJets.size(); ijet++) {
            theHistoManager->fillHisto("JetPt",  "PreSel", "", "", vJet->at(ijet).pt(), weight);
            theHistoManager->fillHisto("JetEta", "PreSel", "", "", vJet->at(ijet).eta(), weight);
            theHistoManager->fillHisto("CSVv2",  "PreSel", "", "", vJet->at(ijet).CSVv2(), weight);
        }
        theHistoManager->fillHisto("METpx",  "PreSel", "", "", METx, weight);
        theHistoManager->fillHisto("METpy",  "PreSel", "", "", METy, weight);
        theHistoManager->fillHisto("MET"  ,  "PreSel", "", "", MET, weight);
        theHistoManager->fillHisto("METphi", "PreSel", "", "", METphi, weight);
        theHistoManager->fillHisto("METsum", "PreSel", "", "", METsum, weight);
        theHistoManager->fillHisto("MHT",    "PreSel", "", "", MHT, weight);
        theHistoManager->fillHisto("MET LD", "PreSel", "", "", met_ld, weight);

        int lepid2 = -1, lepid3 = -1;
        if ( nlepsel >= 2 ) {
            if ( abs(vSelectedLeptons.at(0).id()) == 11 
                    && abs(vSelectedLeptons.at(1).id()) == 11 ) lepid2 = 0; // ee
            if ( abs(vSelectedLeptons.at(0).id()) == 11 
                    && abs(vSelectedLeptons.at(1).id()) == 13 ) lepid2 = 1; // emu
            if ( abs(vSelectedLeptons.at(0).id()) == 13 
                    && abs(vSelectedLeptons.at(1).id()) == 11 ) lepid2 = 1;
            if ( abs(vSelectedLeptons.at(0).id()) == 13 
                    && abs(vSelectedLeptons.at(1).id()) == 13 ) lepid2 = 2; // mumu
        }
        if ( nlepsel >= 3 ) {
            if ( lepid2 == 0 && abs(vSelectedLeptons.at(2).id()) == 11 ) lepid3 = 0; // eee
            if ( lepid2 == 0 && abs(vSelectedLeptons.at(2).id()) == 13 ) lepid3 = 1; // eemu
            if ( lepid2 == 1 && abs(vSelectedLeptons.at(2).id()) == 11 ) lepid3 = 1;
            if ( lepid2 == 1 && abs(vSelectedLeptons.at(2).id()) == 13 ) lepid3 = 2; // emumu
            if ( lepid2 == 2 && abs(vSelectedLeptons.at(2).id()) == 11 ) lepid3 = 2;
            if ( lepid2 == 2 && abs(vSelectedLeptons.at(2).id()) == 13 ) lepid3 = 3; // mumumu
        }

        // BDT variables
        float dr;
        float lep1_mtw = -1., lep1_dr_min = 100., lep2_dr_min = 100., lep_eta_max = -1.;

        if ( vSelectedLeptons.size() >= 2
                && vSelectedLeptons.at(0).pt() > 20 && vSelectedLeptons.at(1).pt() > 10 ) {
            lep_eta_max = fabs(vSelectedLeptons.at(0).eta());
            if ( fabs(vSelectedLeptons.at(1).eta()) > lep_eta_max ) lep_eta_max = fabs(vSelectedLeptons.at(1).eta());
            for (int ijet=0; ijet < vJet->size() ; ijet++) {
                dr = GetDeltaR( vSelectedLeptons.at(0).eta(), vSelectedLeptons.at(0).phi(), vJet->at(ijet).eta(), vJet->at(ijet).phi() );
                if ( dr < lep1_dr_min ) lep1_dr_min = dr;
                dr = GetDeltaR( vSelectedLeptons.at(1).eta(), vSelectedLeptons.at(1).phi(), vJet->at(ijet).eta(), vJet->at(ijet).phi() );
                if ( dr < lep2_dr_min ) lep2_dr_min = dr;
            }
            lep1_mtw = sqrt( 2 * vSelectedLeptons.at(0).p4().Pt() * MET * (1 - cos( vSelectedLeptons.at(0).phi() - METphi )));
        }	   

        int njj = 0; 
        float jet_dr_av = 0.;
        for (int ijet=0; ijet < vJet->size()-1 ; ijet++) {
            for (int kjet=ijet+1; kjet < vJet->size() ; kjet++) {
                jet_dr_av += GetDeltaR( vJet->at(ijet).eta(), vJet->at(ijet).phi(), vJet->at(kjet).eta(), vJet->at(kjet).phi() );
                njj++;
            }
        }
        if ( njj > 0 ) jet_dr_av = jet_dr_av / njj;

        // TTH3l
        if ( is_3l_TTH_SR ) {

            if ( _isdata ) {
                std::cout << "Run Event " << vEvent->at(0).run() <<  " " << vEvent->at(0).id() 
                    << " nlep/tight/fake " << vLeptons.size() <<  "/" << vSelectedLeptons.size() <<  "/" <<  vFakeLeptons.size() 
                    << " njet/loose/medium " << vSelectedJets.size() 
                    <<  "/" << nLooseBJets <<  "/" << nMediumBJets << std::endl;
                std::cout << "lep id " << lepid3 << " pT " << vSelectedLeptons.at(0).pt() <<  " " 
                    << vSelectedLeptons.at(1).pt() <<  " "
                    << vSelectedLeptons.at(2).pt() <<  " " << std::endl;
                std::cout << " " << std::endl;
            }

            theHistoManager->fillHisto("nLep",   "TTH3l", "", "", vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("nFake",  "TTH3l", "", "", vFakeLeptons.size()+vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("lep1Pt", "TTH3l", "", "", vSelectedLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep2Pt", "TTH3l", "", "", vSelectedLeptons.at(1).pt(), weight);
            theHistoManager->fillHisto("lep3Pt", "TTH3l", "", "", vSelectedLeptons.at(2).pt(), weight);
            theHistoManager->fillHisto("lep1Eta", "TTH3l", "", "", vSelectedLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lep2Eta", "TTH3l", "", "", vSelectedLeptons.at(1).eta(), weight);
            theHistoManager->fillHisto("lep3Eta", "TTH3l", "", "", vSelectedLeptons.at(2).eta(), weight);
            if ( vSelectedLeptons.size() >= 4 ) {
                theHistoManager->fillHisto("lep4Pt",  "TTH3l", "", "", vSelectedLeptons.at(3).pt(), weight);
                theHistoManager->fillHisto("lep4Eta", "TTH3l", "", "", vSelectedLeptons.at(3).eta(), weight);
                theHistoManager->fillHisto("lepQ",    "TTH3l", "", "", vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge()+vSelectedLeptons.at(2).charge()+vSelectedLeptons.at(3).charge(), weight);
            }
            else theHistoManager->fillHisto("lepQ", "TTH3l", "", "", vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge()+vSelectedLeptons.at(2).charge(), weight);
            theHistoManager->fillHisto("nJets",     "TTH3l", "", "", vSelectedJets.size(), weight);
            theHistoManager->fillHisto("nLooseB",   "TTH3l", "", "", nLooseBJets, weight);
            theHistoManager->fillHisto("nMediumB",  "TTH3l", "", "", nMediumBJets, weight);
            for (unsigned int ijet=0; ijet < vJet->size() ; ijet++) {
                theHistoManager->fillHisto("JetPt",  "TTH3l", "", "", vJet->at(ijet).pt(), weight);
                theHistoManager->fillHisto("JetEta", "TTH3l", "", "", vJet->at(ijet).eta(), weight);
            }
            theHistoManager->fillHisto("METpx",  "TTH3l", "", "", METx, weight);
            theHistoManager->fillHisto("METpy",  "TTH3l", "", "", METy, weight);
            theHistoManager->fillHisto("MET"  ,  "TTH3l", "", "", MET, weight);
            theHistoManager->fillHisto("METphi", "TTH3l", "", "", METphi, weight);
            theHistoManager->fillHisto("METsum", "TTH3l", "", "", METsum, weight);
            theHistoManager->fillHisto("MHT",    "TTH3l", "", "", MHT, weight);
            if ( vSelectedJets.size() <= 3 ) theHistoManager->fillHisto("MET LD", "TTH3l", "", "", met_ld, weight);
            theHistoManager->fillHisto("lep1 Mt",     "TTH3l", "", "", lep1_mtw, weight);
            theHistoManager->fillHisto("lep1 dRmin",  "TTH3l", "", "", lep1_dr_min, weight);
            theHistoManager->fillHisto("lep2 dRmin",  "TTH3l", "", "", lep2_dr_min, weight);
            theHistoManager->fillHisto("lep Eta max", "TTH3l", "", "", lep_eta_max, weight);
            theHistoManager->fillHisto("jet dR av",   "TTH3l", "", "", jet_dr_av, weight);
            theHistoManager->fillHisto("LepId", "TTH3l", "", "", lepid3, weight);
            theHistoManager->fillHisto("HT",     "TTH3l", "", "", jetht, weight);
            theHistoManager->fillHisto("Mll min","TTH3l", "", "", Mllmin, weight);
            theHistoManager->fillHisto("best MZ","TTH3l", "", "", Mllbest, weight);
        }

        // WZ
        if ( is_3l_WZ_CR ) {
            theHistoManager->fillHisto("nLep",   "WZ", "", "", vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("nFake",  "WZ", "", "", vFakeLeptons.size()+vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("lep1Pt", "WZ", "", "", vSelectedLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep2Pt", "WZ", "", "", vSelectedLeptons.at(1).pt(), weight);
            theHistoManager->fillHisto("lep3Pt", "WZ", "", "", vSelectedLeptons.at(2).pt(), weight);
            theHistoManager->fillHisto("lep1Eta", "WZ", "", "", vSelectedLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lep2Eta", "WZ", "", "", vSelectedLeptons.at(1).eta(), weight);
            theHistoManager->fillHisto("lep3Eta", "WZ", "", "", vSelectedLeptons.at(2).eta(), weight);
            theHistoManager->fillHisto("lepQ", "WZ", "", "", vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge()+vSelectedLeptons.at(2).charge(), weight);
            int lepw = -1;
            float zpt = -1.;
            for (int i=0; i<vSelectedLeptons.size()-1; i++) {
                for (int j=i+1; j<vSelectedLeptons.size(); j++) {
                    if ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id()) {
                        float mz = ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M();
                        if ( fabs(mz - 91.188) < 10. ) {
                            if ( i == 0 && j == 1 ) lepw = 2;
                            if ( i == 0 && j == 2 ) lepw = 1;    
                            if ( i == 1 && j == 2 ) lepw = 0;
                            zpt = ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).Pt();
                        }
                    }
                }
            }
            theHistoManager->fillHisto("nJets",     "WZ", "", "", vSelectedJets.size(), weight);
            theHistoManager->fillHisto("nLooseB",   "WZ", "", "", nLooseBJets, weight);
            theHistoManager->fillHisto("nMediumB",  "WZ", "", "", nMediumBJets, weight);
            for (unsigned int ijet=0; ijet < vJet->size() ; ijet++) {
                theHistoManager->fillHisto("JetPt",  "WZ", "", "", vJet->at(ijet).pt(), weight);
                theHistoManager->fillHisto("JetEta", "WZ", "", "", vJet->at(ijet).eta(), weight);
            }
            theHistoManager->fillHisto("METpx",  "WZ", "", "", METx, weight);
            theHistoManager->fillHisto("METpy",  "WZ", "", "", METy, weight);
            theHistoManager->fillHisto("MET"  ,  "WZ", "", "", MET, weight);
            theHistoManager->fillHisto("METphi", "WZ", "", "", METphi, weight);
            theHistoManager->fillHisto("METsum", "WZ", "", "", METsum, weight);
            theHistoManager->fillHisto("MHT",    "WZ", "", "", MHT, weight);
            if ( vSelectedJets.size() <= 3 ) theHistoManager->fillHisto("MET LD", "WZ", "", "", met_ld, weight);
            theHistoManager->fillHisto("lep1 Mt",     "WZ", "", "", lep1_mtw, weight);
            theHistoManager->fillHisto("lep1 dRmin",  "WZ", "", "", lep1_dr_min, weight);
            theHistoManager->fillHisto("lep2 dRmin",  "WZ", "", "", lep2_dr_min, weight);
            theHistoManager->fillHisto("lep Eta max", "WZ", "", "", lep_eta_max, weight);
            theHistoManager->fillHisto("jet dR av",   "WZ", "", "", jet_dr_av, weight);
            if ( lepw >= 0 ) {
                float mtw = sqrt( 2 * vSelectedLeptons.at(lepw).p4().Pt() * MET * (1 - cos( vSelectedLeptons.at(lepw).phi() - METphi )));
                theHistoManager->fillHisto("W Mt", "WZ", "", "", mtw, weight);
                theHistoManager->fillHisto("Z Pt", "WZ", "", "", zpt, weight);
            }
            TLorentzVector all_lep_invmass_p4;
            for (int i=0; i<vSelectedLeptons.size(); i++) {
                all_lep_invmass_p4 += vSelectedLeptons.at(i).p4();
            }
            float all_lep_invmass = all_lep_invmass_p4.M();
            float all_lep_sumofpt = sqrt( (lep_px*lep_px) + (lep_py*lep_py) );
            theHistoManager->fillHisto("inv.mass(l)", "WZ", "", "", all_lep_invmass, weight);
            theHistoManager->fillHisto("Pt Sum(l)",   "WZ", "", "", all_lep_sumofpt, weight);
            theHistoManager->fillHisto("LepId", "WZ", "", "", lepid3, weight);
            theHistoManager->fillHisto("HT",     "WZ", "", "", jetht, weight);
            theHistoManager->fillHisto("Mll min","WZ", "", "", Mllmin, weight);
            theHistoManager->fillHisto("best MZ","WZ", "", "", Mllbest, weight);
        }

        // WZrelaxed
        if ( is_3l_WZrel_CR ) {
            theHistoManager->fillHisto("nLep",   "WZrelaxed", "", "", vLeptons.size(), weight);
            theHistoManager->fillHisto("nFake",  "WZrelaxed", "", "", vFakeLeptons.size()+vLeptons.size(), weight);
            theHistoManager->fillHisto("lep1Pt", "WZrelaxed", "", "", vLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep2Pt", "WZrelaxed", "", "", vLeptons.at(1).pt(), weight);
            theHistoManager->fillHisto("lep3Pt", "WZrelaxed", "", "", vLeptons.at(2).pt(), weight);
            theHistoManager->fillHisto("lep1Eta", "WZrelaxed", "", "", vLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lep2Eta", "WZrelaxed", "", "", vLeptons.at(1).eta(), weight);
            theHistoManager->fillHisto("lep3Eta", "WZrelaxed", "", "", vLeptons.at(2).eta(), weight);
            theHistoManager->fillHisto("lepQ", "WZrelaxed", "", "", vLeptons.at(0).charge()+vLeptons.at(1).charge()+vLeptons.at(2).charge(), weight);
            int lepw = -1;
            float zpt = -1.;
            for (int i=0; i<vLeptons.size()-1; i++) {
                for (int j=i+1; j<vLeptons.size(); j++) {
                    if ( vLeptons.at(i).id() == -vLeptons.at(j).id()) {
                        float mz = ( vLeptons.at(i).p4() + vLeptons.at(j).p4() ).M();
                        if ( fabs(mz - 91.188) < 10. ) {
                            if ( i == 0 && j == 1 ) lepw = 2;
                            if ( i == 0 && j == 2 ) lepw = 1;    
                            if ( i == 1 && j == 2 ) lepw = 0;
                            zpt = ( vLeptons.at(i).p4() + vLeptons.at(j).p4() ).Pt();
                        }
                    }
                }
            }
            theHistoManager->fillHisto("nJets",     "WZrelaxed", "", "", vSelectedJets.size(), weight);
            theHistoManager->fillHisto("nLooseB",   "WZrelaxed", "", "", nLooseBJets, weight);
            theHistoManager->fillHisto("nMediumB",  "WZrelaxed", "", "", nMediumBJets, weight);
            for (unsigned int ijet=0; ijet < vJet->size() ; ijet++) {
                theHistoManager->fillHisto("JetPt",  "WZrelaxed", "", "", vJet->at(ijet).pt(), weight);
                theHistoManager->fillHisto("JetEta", "WZrelaxed", "", "", vJet->at(ijet).eta(), weight);
            }
            theHistoManager->fillHisto("METpx",  "WZrelaxed", "", "", METx, weight);
            theHistoManager->fillHisto("METpy",  "WZrelaxed", "", "", METy, weight);
            theHistoManager->fillHisto("MET"  ,  "WZrelaxed", "", "", MET, weight);
            theHistoManager->fillHisto("METphi", "WZrelaxed", "", "", METphi, weight);
            theHistoManager->fillHisto("METsum", "WZrelaxed", "", "", METsum, weight);
            theHistoManager->fillHisto("MHT",    "WZrelaxed", "", "", MHT, weight);
            if ( vSelectedJets.size() <= 3 ) theHistoManager->fillHisto("MET LD", "WZrelaxed", "", "", met_ld, weight);
            theHistoManager->fillHisto("lep1 Mt",     "WZrelaxed", "", "", lep1_mtw, weight);
            theHistoManager->fillHisto("lep1 dRmin",  "WZrelaxed", "", "", lep1_dr_min, weight);
            theHistoManager->fillHisto("lep2 dRmin",  "WZrelaxed", "", "", lep2_dr_min, weight);
            theHistoManager->fillHisto("lep Eta max", "WZrelaxed", "", "", lep_eta_max, weight);
            theHistoManager->fillHisto("jet dR av",   "WZrelaxed", "", "", jet_dr_av, weight);
            if ( lepw >= 0 ) {
                float mtw = sqrt( 2 * vLeptons.at(lepw).p4().Pt() * MET * (1 - cos( vLeptons.at(lepw).phi() - METphi )));
                theHistoManager->fillHisto("W Mt", "WZrelaxed", "", "", mtw, weight);
                theHistoManager->fillHisto("Z Pt", "WZrelaxed", "", "", zpt, weight);
            }
            TLorentzVector all_lep_invmass_p4;
            for (int i=0; i<vLeptons.size(); i++) {
                all_lep_invmass_p4 += vLeptons.at(i).p4();
            }
            float all_lep_invmass = all_lep_invmass_p4.M();
            float all_lep_sumofpt = sqrt( (lep_px*lep_px) + (lep_py*lep_py) );
            theHistoManager->fillHisto("inv.mass(l)", "WZrelaxed", "", "", all_lep_invmass, weight);
            theHistoManager->fillHisto("Pt Sum(l)",   "WZrelaxed", "", "", all_lep_sumofpt, weight);
            theHistoManager->fillHisto("LepId", "WZrelaxed", "", "", lepid3, weight);
            theHistoManager->fillHisto("HT",     "WZrelaxed", "", "", jetht, weight);
            theHistoManager->fillHisto("Mll min","WZrelaxed", "", "", Mllmin, weight);
            theHistoManager->fillHisto("best MZ","WZrelaxed", "", "", Mllbest, weight);
        }

        // TTZ
        if ( is_3l_TTZ_CR ) {
            theHistoManager->fillHisto("nLep",    "TTZ", "", "", vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("nFake",   "TTZ", "", "", vFakeLeptons.size()+vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("lep1Pt",  "TTZ", "", "", vSelectedLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep2Pt",  "TTZ", "", "", vSelectedLeptons.at(1).pt(), weight);
            theHistoManager->fillHisto("lep3Pt",  "TTZ", "", "", vSelectedLeptons.at(2).pt(), weight);
            theHistoManager->fillHisto("lep1Eta", "TTZ", "", "", vSelectedLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lep2Eta", "TTZ", "", "", vSelectedLeptons.at(1).eta(), weight);
            theHistoManager->fillHisto("lep3Eta", "TTZ", "", "", vSelectedLeptons.at(2).eta(), weight);
            if ( vSelectedLeptons.size() >= 4 ) {
                theHistoManager->fillHisto("lep4Pt",  "TTZ", "", "", vSelectedLeptons.at(3).pt(), weight);
                theHistoManager->fillHisto("lep4Eta", "TTZ", "", "", vSelectedLeptons.at(3).eta(), weight);
                theHistoManager->fillHisto("lepQ",    "TTZ", "", "", vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge()+vSelectedLeptons.at(2).charge()+vSelectedLeptons.at(3).charge(), weight);
            }
            else theHistoManager->fillHisto("lepQ", "TTZ", "", "", vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge()+vSelectedLeptons.at(2).charge(), weight);
            theHistoManager->fillHisto("nJets",     "TTZ", "", "", vSelectedJets.size(), weight);
            theHistoManager->fillHisto("nLooseB",   "TTZ", "", "", nLooseBJets, weight);
            theHistoManager->fillHisto("nMediumB",  "TTZ", "", "", nMediumBJets, weight);
            for (unsigned int ijet=0; ijet < vJet->size() ; ijet++) {
                theHistoManager->fillHisto("JetPt",  "TTZ", "", "", vJet->at(ijet).pt(), weight);
                theHistoManager->fillHisto("JetEta", "TTZ", "", "", vJet->at(ijet).eta(), weight);
            }
            theHistoManager->fillHisto("METpx",  "TTZ", "", "", METx, weight);
            theHistoManager->fillHisto("METpy",  "TTZ", "", "", METy, weight);
            theHistoManager->fillHisto("MET"  ,  "TTZ", "", "", MET, weight);
            theHistoManager->fillHisto("METphi", "TTZ", "", "", METphi, weight);
            theHistoManager->fillHisto("METsum", "TTZ", "", "", METsum, weight);
            theHistoManager->fillHisto("MHT",    "TTZ", "", "", MHT, weight);
            if ( vSelectedJets.size() <= 3 ) theHistoManager->fillHisto("MET LD", "TTZ", "", "", met_ld, weight);
            theHistoManager->fillHisto("lep1 Mt",     "TTZ", "", "", lep1_mtw, weight);
            theHistoManager->fillHisto("lep1 dRmin",  "TTZ", "", "", lep1_dr_min, weight);
            theHistoManager->fillHisto("lep2 dRmin",  "TTZ", "", "", lep2_dr_min, weight);
            theHistoManager->fillHisto("lep Eta max", "TTZ", "", "", lep_eta_max, weight);
            theHistoManager->fillHisto("jet dR av",   "TTZ", "", "", jet_dr_av, weight);
            theHistoManager->fillHisto("LepId", "TTZ", "", "", lepid3, weight);
            theHistoManager->fillHisto("HT",     "TTZ", "", "", jetht, weight);
            theHistoManager->fillHisto("Mll min","TTZ", "", "", Mllmin, weight);
            theHistoManager->fillHisto("best MZ","TTZ", "", "", Mllbest, weight);
        }

        // Zl
        if ( is_Zl_CR ) {
            theHistoManager->fillHisto("nLep",   "Zl", "", "", vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("nFake",  "Zl","", "",  vFakeLeptons.size()+vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("lep1Pt", "Zl", "", "", vSelectedLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep2Pt", "Zl", "", "", vSelectedLeptons.at(1).pt(), weight);
            theHistoManager->fillHisto("lep3Pt", "Zl", "", "", vFakeLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep1Eta", "Zl", "", "", vSelectedLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lep2Eta", "Zl", "", "", vSelectedLeptons.at(1).eta(), weight);
            theHistoManager->fillHisto("lep3Eta", "Zl", "", "", vFakeLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lepQ",    "Zl", "", "", vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge(), weight);
            theHistoManager->fillHisto("nJets",     "Zl", "", "", vSelectedJets.size(), weight);
            theHistoManager->fillHisto("nLooseB",   "Zl", "", "", nLooseBJets, weight);
            theHistoManager->fillHisto("nMediumB",  "Zl", "", "", nMediumBJets, weight);
            for (unsigned int ijet=0; ijet < vJet->size() ; ijet++) {
                theHistoManager->fillHisto("JetPt",  "Zl", "", "", vJet->at(ijet).pt(), weight);
                theHistoManager->fillHisto("JetEta", "Zl", "", "", vJet->at(ijet).eta(), weight);
            }
            theHistoManager->fillHisto("METpx",  "Zl", "", "", METx, weight);
            theHistoManager->fillHisto("METpy",  "Zl", "", "", METy, weight);
            theHistoManager->fillHisto("MET"  ,  "Zl", "", "", MET, weight);
            theHistoManager->fillHisto("METphi", "Zl", "", "", METphi, weight);
            theHistoManager->fillHisto("METsum", "Zl", "", "", METsum, weight);
            theHistoManager->fillHisto("MHT",    "Zl", "", "", MHT, weight);
            if ( vSelectedJets.size() <= 3 ) theHistoManager->fillHisto("MET LD", "Zl", "", "", met_ld, weight);
            theHistoManager->fillHisto("lep1 Mt",     "Zl", "", "", lep1_mtw, weight);
            theHistoManager->fillHisto("lep1 dRmin",  "Zl", "", "", lep1_dr_min, weight);
            theHistoManager->fillHisto("lep2 dRmin",  "Zl", "", "", lep2_dr_min, weight);
            theHistoManager->fillHisto("lep Eta max", "Zl", "", "", lep_eta_max, weight);
            theHistoManager->fillHisto("jet dR av",   "Zl", "", "", jet_dr_av, weight);
            theHistoManager->fillHisto("LepId", "Zl", "", "", lepid2, weight);
            theHistoManager->fillHisto("HT",     "Zl", "", "", jetht, weight);
            theHistoManager->fillHisto("Mll min","Zl", "", "", Mllmin, weight);
            theHistoManager->fillHisto("best MZ","Zl", "", "", Mllbest, weight);
        }

        // TTH2l
        if ( is_2lss_TTH_SR ) {
            theHistoManager->fillHisto("nLep",   "TTH2l", "", "", vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("nFake",  "TTH2l","", "",  vFakeLeptons.size()+vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("lep1Pt", "TTH2l", "", "", vSelectedLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep2Pt", "TTH2l", "", "", vSelectedLeptons.at(1).pt(), weight);
            theHistoManager->fillHisto("lep1Eta", "TTH2l", "", "", vSelectedLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lep2Eta", "TTH2l", "", "", vSelectedLeptons.at(1).eta(), weight);
            theHistoManager->fillHisto("lepQ", "TTH2l", "", "", vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge(), weight);
            theHistoManager->fillHisto("nJets",     "TTH2l", "", "", vSelectedJets.size(), weight);
            theHistoManager->fillHisto("nLooseB",   "TTH2l", "", "", nLooseBJets, weight);
            theHistoManager->fillHisto("nMediumB",  "TTH2l", "", "", nMediumBJets, weight);
            for (unsigned int ijet=0; ijet < vJet->size() ; ijet++) {
                theHistoManager->fillHisto("JetPt",  "TTH2l", "", "", vJet->at(ijet).pt(), weight);
                theHistoManager->fillHisto("JetEta", "TTH2l", "", "", vJet->at(ijet).eta(), weight);
            }
            theHistoManager->fillHisto("METpx",  "TTH2l", "", "", METx, weight);
            theHistoManager->fillHisto("METpy",  "TTH2l", "", "", METy, weight);
            theHistoManager->fillHisto("MET"  ,  "TTH2l", "", "", MET, weight);
            theHistoManager->fillHisto("METphi", "TTH2l", "", "", METphi, weight);
            theHistoManager->fillHisto("METsum", "TTH2l", "", "", METsum, weight);
            theHistoManager->fillHisto("MHT",    "TTH2l", "", "", MHT, weight);
            if ( vSelectedJets.size() <= 3 ) theHistoManager->fillHisto("MET LD", "TTH2l", "", "", met_ld, weight);
            theHistoManager->fillHisto("lep1 Mt",     "TTH2l", "", "", lep1_mtw, weight);
            theHistoManager->fillHisto("lep1 dRmin",  "TTH2l", "", "", lep1_dr_min, weight);
            theHistoManager->fillHisto("lep2 dRmin",  "TTH2l", "", "", lep2_dr_min, weight);
            theHistoManager->fillHisto("lep Eta max", "TTH2l", "", "", lep_eta_max, weight);
            theHistoManager->fillHisto("jet dR av",   "TTH2l", "", "", jet_dr_av, weight);
            theHistoManager->fillHisto("LepId", "TTH2l", "", "", lepid2, weight);
            theHistoManager->fillHisto("HT",     "TTH2l", "", "", jetht, weight);
            theHistoManager->fillHisto("Mll min","TTH2l", "", "", Mllmin, weight);
            theHistoManager->fillHisto("best MZ","TTH2l", "", "", Mllbest, weight);
        }

        // TTdilep
        if ( is_emu_TT_CR ) {
            theHistoManager->fillHisto("nLep",   "TTemu", "", "", vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("nFake",  "TTemu","", "",  vFakeLeptons.size()+vSelectedLeptons.size(), weight);
            theHistoManager->fillHisto("lep1Pt", "TTemu", "", "", vSelectedLeptons.at(0).pt(), weight);
            theHistoManager->fillHisto("lep2Pt", "TTemu", "", "", vSelectedLeptons.at(1).pt(), weight);
            theHistoManager->fillHisto("lep1Eta", "TTemu", "", "", vSelectedLeptons.at(0).eta(), weight);
            theHistoManager->fillHisto("lep2Eta", "TTemu", "", "", vSelectedLeptons.at(1).eta(), weight);
            theHistoManager->fillHisto("lepQ", "TTemu", "", "", vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge(), weight);
            theHistoManager->fillHisto("nJets",     "TTemu", "", "", vSelectedJets.size(), weight);
            theHistoManager->fillHisto("nLooseB",   "TTemu", "", "", nLooseBJets, weight);
            theHistoManager->fillHisto("nMediumB",  "TTemu", "", "", nMediumBJets, weight);
            for (unsigned int ijet=0; ijet < vJet->size() ; ijet++) {
                theHistoManager->fillHisto("JetPt",  "TTemu", "", "", vJet->at(ijet).pt(), weight);
                theHistoManager->fillHisto("JetEta", "TTemu", "", "", vJet->at(ijet).eta(), weight);
            }
            theHistoManager->fillHisto("METpx",  "TTemu", "", "", METx, weight);
            theHistoManager->fillHisto("METpy",  "TTemu", "", "", METy, weight);
            theHistoManager->fillHisto("MET"  ,  "TTemu", "", "", MET, weight);
            theHistoManager->fillHisto("METphi", "TTemu", "", "", METphi, weight);
            theHistoManager->fillHisto("METsum", "TTemu", "", "", METsum, weight);
            theHistoManager->fillHisto("MHT",    "TTemu", "", "", MHT, weight);
            if ( vSelectedJets.size() <= 3 ) theHistoManager->fillHisto("MET LD", "TTemu", "", "", met_ld, weight);
            theHistoManager->fillHisto("lep1 Mt",     "TTemu", "", "", lep1_mtw, weight);
            theHistoManager->fillHisto("lep1 dRmin",  "TTemu", "", "", lep1_dr_min, weight);
            theHistoManager->fillHisto("lep2 dRmin",  "TTemu", "", "", lep2_dr_min, weight);
            theHistoManager->fillHisto("lep Eta max", "TTemu", "", "", lep_eta_max, weight);
            theHistoManager->fillHisto("jet dR av",   "TTemu", "", "", jet_dr_av, weight);
            theHistoManager->fillHisto("LepId", "TTemu", "", "", lepid2, weight);
            theHistoManager->fillHisto("HT",     "TTemu", "", "", jetht, weight);
            theHistoManager->fillHisto("Mll min","TTemu", "", "", Mllmin, weight);
            theHistoManager->fillHisto("best MZ","TTemu", "", "", Mllbest, weight);
        }

    }

}

void TTbarHiggsMultileptonAnalysis::TwoLeptonsSameSignSelection_TTH2l(int evt)
{

    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel",   "ttH2lss",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel",   "ttH2lss",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection #
    // ####################

    bool nLep               = ( vSelectedLeptons.size()     == 2 );
    if(!nLep)               return;

    bool leading_lep_pt     = ( vSelectedLeptons.at(0).pt() > 20 );
    if(!leading_lep_pt)     return;

    //bool following_lep_pt   = ( vSelectedLeptons.at(1).pt()               > 10 );
    bool following_lep_pt   = (  ( (abs(vSelectedLeptons.at(1).id()) == 11) && (vSelectedLeptons.at(1).pt() > 15) )
            || ( (abs(vSelectedLeptons.at(1).id()) == 13) && (vSelectedLeptons.at(1).pt() > 10) ) );
    if(!following_lep_pt)   return;

    bool passMll12Gt12      = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M()  > 12);
    if(!passMll12Gt12)      return;

    bool nJets              = ( vSelectedJets.size()        >= 4 );
    if(!nJets)              return;

    bool nLooseBtag         = ( nLooseBJets                 >= 2 );
    bool nMediumBtag        = ( nMediumBJets                >= 1 );
    if(!nLooseBtag && !nMediumBtag)      return;


    // ###############################
    // # Two leptons event selection #
    // ###############################

    bool same_sign      = ( vSelectedLeptons.at(0).charge() == vSelectedLeptons.at(1).charge() );
    if(!same_sign)      return;

    // ##########
    // # Z veto # here for leptons of same charge only !
    // ##########

    bool pass_Zveto = true;
    for(int i=0; i<vSelectedLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vSelectedLeptons.size(); j++)
        {
            if (  ( vSelectedLeptons.at(i).id() == vSelectedLeptons.at(j).id()                            )
                    && ( fabs( ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M() - 91.188) < 10 ) )
            { pass_Zveto = false ;}
        }
    }

    theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "ttH2lss",   "",   2, weight);

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "ttH2lss", _sampleName.Data(),   3, weight);

    if(  (abs(vSelectedLeptons.at(0).id()) == 11)
            && (abs(vSelectedLeptons.at(1).id()) == 11)
            && (pass_Zveto                            )
            && (met_ld                           > 0.2) // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "FullThreeLeptons",   "ttH2lee",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "ttH2lee",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "ttH2lee",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "FullThreeLeptons",   "ttH2lee",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "FullThreeLeptons",   "ttH2lee",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "FullThreeLeptons",   "ttH2lee",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "FullThreeLeptons",   "ttH2lee",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "FullThreeLeptons",   "ttH2lee",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "ttH2lee",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "ttH2lee",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "ttH2lee",   "", vSelectedLeptons.size()        , weight);
    }

    if(  ( (abs(vSelectedLeptons.at(0).id()) == 11) && (abs(vSelectedLeptons.at(1).id()) == 13) )
            || ( (abs(vSelectedLeptons.at(0).id()) == 13) && (abs(vSelectedLeptons.at(1).id()) == 11) ) // smarter way to do is probably exists
            && true                                      // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "FullThreeLeptons",   "ttH2lem",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "ttH2lem",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "ttH2lem",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "FullThreeLeptons",   "ttH2lem",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "FullThreeLeptons",   "ttH2lem",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "FullThreeLeptons",   "ttH2lem",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "FullThreeLeptons",   "ttH2lem",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "FullThreeLeptons",   "ttH2lem",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "ttH2lem",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "ttH2lem",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "ttH2lem",   "", vSelectedLeptons.size()        , weight);
    }

    if(  (abs(vSelectedLeptons.at(0).id()) == 13)
            && (abs(vSelectedLeptons.at(1).id()) == 13)
            && true                                      // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "FullThreeLeptons",   "ttH2lmm",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "ttH2lmm",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "ttH2lmm",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "FullThreeLeptons",   "ttH2lmm",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "FullThreeLeptons",   "ttH2lmm",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "FullThreeLeptons",   "ttH2lmm",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "FullThreeLeptons",   "ttH2lmm",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "FullThreeLeptons",   "ttH2lmm",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "ttH2lmm",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "ttH2lmm",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "ttH2lmm",   "", vSelectedLeptons.size()        , weight);
    }

    if (   (abs(vSelectedLeptons.at(0).id()) == 11)
            && (abs(vSelectedLeptons.at(1).id()) == 11)
            && ( !pass_Zveto || met_ld < 0.2)         ) return;

    is_2lss_TTH_SR = true;   

    // Calcul of input variables of the 2D BDT
    max_Lep_eta     = std::max( abs(vSelectedLeptons.at(0).eta()), abs(vSelectedLeptons.at(1).eta()) ) ;

    numJets_float   = vSelectedJets.size() ;

    mindr_lep1_jet = 1000.;

    for (int i=0; i<vSelectedJets.size(); i++)
    {
        if( DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ) < mindr_lep1_jet )
        { mindr_lep1_jet = DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ); }
    }

    mindr_lep2_jet  = 1000. ;

    for (int i=0; i<vSelectedJets.size(); i++)
    {
        if( DeltaRLeptonJet( vSelectedLeptons.at(1), vSelectedJets.at(i) ) < mindr_lep1_jet )
        { mindr_lep1_jet = DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ); }
    }

    met             = vEvent->at(0).metpt() ;

    //avg_dr_jet      = DeltaRJets( vSelectedNonBTagJets.at(0), vSelectedNonBTagJets.at(1) ) ;
    avg_dr_jet      = DeltaRJets( vSelectedJets.at(0), vSelectedJets.at(1) ) ;

    MT_met_lep1     = sqrt( 2 * vSelectedLeptons.at(0).p4().Pt() * vEvent->at(0).metpt() * (1 - cos( vSelectedLeptons.at(0).phi() - vEvent->at(0).metphi() )));

    LepGood_conePt0 = vSelectedLeptons.at(0).pt() ;

    LepGood_conePt1 = vSelectedLeptons.at(1).pt() ;

    mhtJet25_Recl   = sqrt( (jet_px*jet_px) + (jet_py*jet_py) );

    signal_2lss_TT_MVA  = mva_2lss_tt->EvaluateMVA("BDTG method");
    signal_2lss_TTV_MVA = mva_2lss_ttV->EvaluateMVA("BDTG method");

    //std::cout << " signal 2lss TT MVA: "  << signal_2lss_TT_MVA
    //          << " signal 2lss TTV MVA: " << signal_2lss_TTV_MVA << std::endl;

    theHistoManager->fillHisto("Signal_2lss_TT_MVA",                       "FinalCut", "ttH2lss",   "",  signal_2lss_TT_MVA,   weight);
    theHistoManager->fillHisto("Signal_2lss_TTV_MVA",                      "FinalCut", "ttH2lss",   "",  signal_2lss_TTV_MVA,  weight);


    //fillOutputTree();

    if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);
}

void TTbarHiggsMultileptonAnalysis::TwoLeptonsSameSignSelection_LepMVA_sideband(int evt)
{

    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel",   "LepMVA_2l_SB",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel",   "LepMVA_2l_SB",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection #
    // ####################

    bool nLep               = ( vSelectedLeptons.size()                   == 2 );
    if(!nLep)               return;

    bool leading_lep_pt     = ( vSelectedLeptons.at(0).pt()               > 20 );
    if(!leading_lep_pt)     return;

    //bool following_lep_pt   = ( vSelectedLeptons.at(1).pt()               > 10 );
    bool following_lep_pt   = (  ( (abs(vSelectedLeptons.at(1).id()) == 11) && (vSelectedLeptons.at(1).pt() > 15) )
            || ( (abs(vSelectedLeptons.at(1).id()) == 13) && (vSelectedLeptons.at(1).pt() > 10) ) );
    if(!following_lep_pt)   return;

    bool passMll12Gt12      = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M()  > 12);
    if(!passMll12Gt12)      return;

    bool nJets              = ( vSelectedJets.size()                      >= 4 );
    if(!nJets)              return;

    bool nLooseBtag         = ( nLooseBJets                               >= 2 );
    bool nMediumBtag        = ( nMediumBJets                              >= 1 );
    if(!nLooseBtag && !nMediumBtag)      return;


    // ###############################
    // # Two leptons event selection #
    // ###############################

    // ##########
    // # Z veto #
    // ##########

    bool pass_Zveto = true;
    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        for(int j=0; j<vSelectedLeptons.size(); j++)
        {
            if (  ( i                           != j                                                       )
                    && ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id()                            )
                    && ( fabs( ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M() - 91.188) < 10 ) )
            { pass_Zveto = false ;}
        }
    }

    theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "LepMVA_2l_SB",   "",   2, weight);

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "LepMVA_2l_SB", _sampleName.Data(),   3, weight);

    if(  (abs(vSelectedLeptons.at(0).id()) == 11)
            && (abs(vSelectedLeptons.at(1).id()) == 11)
            && (pass_Zveto                            )
            && (met_ld                           > 0.2) // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "FullThreeLeptons",   "LepMVA_2l_SB",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "FullThreeLeptons",   "LepMVA_2l_SB",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "FullThreeLeptons",   "LepMVA_2l_SB",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "FullThreeLeptons",   "LepMVA_2l_SB",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedLeptons.size()        , weight);
    }

    if(  ( (abs(vSelectedLeptons.at(0).id()) == 11) && (abs(vSelectedLeptons.at(1).id()) == 13) )
            || ( (abs(vSelectedLeptons.at(0).id()) == 13) && (abs(vSelectedLeptons.at(1).id()) == 11) ) // smarter way to do is probably exists
            && true                                      // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "FullThreeLeptons",   "LepMVA_2l_SB",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "FullThreeLeptons",   "LepMVA_2l_SB",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "FullThreeLeptons",   "LepMVA_2l_SB",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "FullThreeLeptons",   "LepMVA_2l_SB",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedLeptons.size()        , weight);
    }

    if(  (abs(vSelectedLeptons.at(0).id()) == 13)
            && (abs(vSelectedLeptons.at(1).id()) == 13)
            && true                                      // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "FullThreeLeptons",   "LepMVA_2l_SB",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "FullThreeLeptons",   "LepMVA_2l_SB",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "FullThreeLeptons",   "LepMVA_2l_SB",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "FullThreeLeptons",   "LepMVA_2l_SB",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "LepMVA_2l_SB",   "", vSelectedLeptons.size()        , weight);
    }

    is_2lss_LepMVA_SB = true;   

    //fillOutputTree();

    if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);
}

void TTbarHiggsMultileptonAnalysis::TwoLeptonsSameSignSelection_JetMultiplicity_sideband(int evt)
{

    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel",   "JM_2l_SB",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel",   "JM_2l_SB",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection #
    // ####################

    bool nLep               = ( vSelectedLeptons.size()                   == 2 );
    if(!nLep)               return;

    bool leading_lep_pt     = ( vSelectedLeptons.at(0).pt()               > 20 );
    if(!leading_lep_pt)     return;

    bool following_lep_pt   = (  ( (abs(vSelectedLeptons.at(1).id()) == 11) && (vSelectedLeptons.at(1).pt() > 15) )
            || ( (abs(vSelectedLeptons.at(1).id()) == 13) && (vSelectedLeptons.at(1).pt() > 10) ) );
    if(!following_lep_pt)   return;

    bool passMll12Gt12      = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M()  > 12);
    if(!passMll12Gt12)      return;

    bool nJets              = ( vSelectedJets.size()                      == 3 );
    if(!nJets)              return;

    bool nLooseBtag         = ( nLooseBJets                               >= 2 );
    bool nMediumBtag        = ( nMediumBJets                              >= 1 );
    if(!nLooseBtag && !nMediumBtag)      return;


    // ###############################
    // # Two leptons event selection #
    // ###############################

    // ##########
    // # Z veto #
    // ##########

    bool pass_Zveto = true;
    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        for(int j=0; j<vSelectedLeptons.size(); j++)
        {
            if (  ( i                           != j                                                       )
                    && ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id()                            )
                    && ( fabs( ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M() - 91.188) < 10 ) )
            { pass_Zveto = false ;}
        }
    }

    theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "JM_2l_SB",   "",   2, weight);

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "JM_2l_SB", _sampleName.Data(),   3, weight);

    if(  (abs(vSelectedLeptons.at(0).id()) == 11)
            && (abs(vSelectedLeptons.at(1).id()) == 11)
            && (pass_Zveto                            )
            && (met_ld                           > 0.2) // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "FullThreeLeptons",   "JM_2l_SB",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "FullThreeLeptons",   "JM_2l_SB",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "FullThreeLeptons",   "JM_2l_SB",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "FullThreeLeptons",   "JM_2l_SB",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedLeptons.size()        , weight);
    }

    if(  ( (abs(vSelectedLeptons.at(0).id()) == 11) && (abs(vSelectedLeptons.at(1).id()) == 13) )
            || ( (abs(vSelectedLeptons.at(0).id()) == 13) && (abs(vSelectedLeptons.at(1).id()) == 11) ) // smarter way to do is probably exists
            && true                                      // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "FullThreeLeptons",   "JM_2l_SB",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "FullThreeLeptons",   "JM_2l_SB",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "FullThreeLeptons",   "JM_2l_SB",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "FullThreeLeptons",   "JM_2l_SB",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedLeptons.size()        , weight);
    }

    if(  (abs(vSelectedLeptons.at(0).id()) == 13)
            && (abs(vSelectedLeptons.at(1).id()) == 13)
            && true                                      // remaining conditions for ee       
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "FullThreeLeptons",   "JM_2l_SB",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "FullThreeLeptons",   "JM_2l_SB",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MHT",                              "FullThreeLeptons",   "JM_2l_SB",   "", MHT                            , weight);
        theHistoManager->fillHisto("MetLD",                            "FullThreeLeptons",   "JM_2l_SB",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "JM_2l_SB",   "", vSelectedLeptons.size()        , weight);
    }

    is_2lss_JM_SB = true;   

    //fillOutputTree();

    if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);
}

void TTbarHiggsMultileptonAnalysis::DiLeptonSelection_TT_CR(int evt)
{

    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel",   "TT_2l_CR",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel",   "TT_2l_CR",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection #
    // ####################

    bool nLep               = ( vSelectedLeptons.size()                   == 2 );
    if(!nLep)               return;

    bool leading_lep_pt     = ( vSelectedLeptons.at(0).pt()               > 20 );
    if(!leading_lep_pt)     return;

    bool following_lep_pt   = (  ( (abs(vSelectedLeptons.at(1).id()) == 11) && (vSelectedLeptons.at(1).pt() > 15) )
            || ( (abs(vSelectedLeptons.at(1).id()) == 13) && (vSelectedLeptons.at(1).pt() > 10) ) );
    if(!following_lep_pt)   return;

    bool passMll12Gt12      = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M()  > 12);
    if(!passMll12Gt12)      return;

    bool nJets              = ( vSelectedJets.size()                      >= 2 );
    if(!nJets)              return;

    bool nMediumBtag        = ( nMediumBJets                              >= 1 );
    if(!nMediumBtag)      return;

    // ###############################
    // #        e+mu- selection      #
    // ###############################

    if ( vSelectedLeptons.at(0).charge()+vSelectedLeptons.at(1).charge() != 0 ) return;
    if ( vSelectedLeptons.at(0).id() == -vSelectedLeptons.at(1).id() ) return;

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if(  ( (abs(vSelectedLeptons.at(0).id()) == 11) && (abs(vSelectedLeptons.at(1).id()) == 13) )
            || ( (abs(vSelectedLeptons.at(0).id()) == 13) && (abs(vSelectedLeptons.at(1).id()) == 11) ) // smarter way to do this probably exists
            && (  vSelectedLeptons.at(0).charge()         == -vSelectedLeptons.at(1).charge()         )
      )
    {
        theHistoManager->fillHisto("CutFlow",                          "FullThreeLeptons",   "TT_2l_CR",   "", 8                              , weight);

        theHistoManager->fillHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "TT_2l_CR",   "", vSelectedLeptons.at(0).pt()    , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "TT_2l_CR",   "", vSelectedLeptons.at(1).pt()    , weight);
        theHistoManager->fillHisto("MET",                              "FullThreeLeptons",   "TT_2l_CR",   "", vEvent->at(0).metpt()          , weight);
        theHistoManager->fillHisto("MetLD",                            "FullThreeLeptons",   "TT_2l_SB",   "", met_ld                         , weight);
        theHistoManager->fillHisto("TauMultiplicity",                  "FullThreeLeptons",   "TT_2l_CR",   "", vSelectedTaus.size()           , weight);
        theHistoManager->fillHisto("JetMultiplicity",                  "FullThreeLeptons",   "TT_2l_CR",   "", vSelectedJets.size()           , weight);
        theHistoManager->fillHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "TT_2l_CR",   "", vSelectedBTagJets.size()       , weight);
        theHistoManager->fillHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "TT_2l_CR",   "", vSelectedMediumBTagJets.size() , weight);
        theHistoManager->fillHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "TT_2l_CR",   "", vSelectedLeptons.size()        , weight);
    }

    is_emu_TT_CR = true;   

    //fillOutputTree();

    if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);
}



void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_TTH3l(int evt)
{

    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel",   "ttH3l",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel",   "ttH3l",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection #
    // ####################

    bool nLep             = ( vSelectedLeptons.size()                   >= 2 );
    if(!nLep)             return;

    bool leading_lep_pt   = ( vSelectedLeptons.at(0).pt()               > 20 );
    if(!leading_lep_pt)   return;

    bool following_lep_pt = ( vSelectedLeptons.at(1).pt()               > 10 );
    if(!following_lep_pt) return;

    bool passMll12Gt12    = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M()  > 12);
    if(!passMll12Gt12)    return;

    bool nJets            = ( vSelectedJets.size()                      >= 2 );
    if(!nJets)            return;

    bool nLooseBtag       = ( nLooseBJets                               >= 2 );
    bool nMediumBtag      = ( nMediumBJets                              >= 1 );
    if(!nLooseBtag && !nMediumBtag)      return;


    // #################################
    // # Three leptons event selection #
    // #################################

    nLep                  = ( vSelectedLeptons.size()				    >= 3 );
    //nLep                  = ( vSelectedLeptons.size()                 == 3 );
    if(!nLep)             return;

    bool third_lep_pt     = ( vSelectedLeptons.at(2).pt()               > 10 );
    if(!third_lep_pt)     return;

    theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "ttH3l",   "",    1, weight);

    passMll12Gt12  = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(2).p4() ).M()  > 12
            && ( vSelectedLeptons.at(1).p4() + vSelectedLeptons.at(2).p4() ).M()  > 12 );
    if(!passMll12Gt12)    return;

    // ##########
    // # Z veto # with tight or loose leptons ???
    // ##########

    bool pass_Zveto = true;
    for(int i=0; i<vSelectedLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vSelectedLeptons.size(); j++)
        {
            if (  ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id()                            )
                    && ( fabs( ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M() - 91.188) < 10 ) )
            { pass_Zveto = false ;}
        }
    }
    if(!pass_Zveto)       return;

    theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "ttH3l",   "",   2, weight);

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "ttH3l", _sampleName.Data(),   3, weight);
    if( vSelectedJets.size() >= 4  )  theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "ttH3l", _sampleName.Data(),   4, weight);

    // ########
    // # SFOS #
    // ########

    bool isSFOS = false;
    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        for(int j=0; j<vSelectedLeptons.size(); j++)
        {
            if (  ( i                           != j                            )
                    && ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id() ) )
            { isSFOS = true ;}
        }
    }

    if( met_ld               > 0.3 ) theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "ttH3l",   "",   6, weight);
    if( isSFOS                     ) theHistoManager->fillHisto("CutFlow",                "FullThreeLeptons",   "ttH3l",   "",   7, weight);

    if(vSelectedJets.size() < 4 && (met_ld < (0.2 + 0.1 * isSFOS)) ) return;

    int sum_charges = 0;
    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        sum_charges = sum_charges + vSelectedLeptons.at(i).charge();
    }

    int sum_charges_3l = 0;
    for(int i=0; i<3; i++)
    {
        sum_charges_3l = sum_charges_3l + vSelectedLeptons.at(i).charge();
    }
    if( fabs(sum_charges_3l) != 1 ) return;

    theHistoManager->fillHisto("CutFlow",                          "FullThreeLeptons",   "ttH3l",   "", 8                              , weight);

    theHistoManager->fillHisto("LeadingLeptonPt",                  "FullThreeLeptons",   "ttH3l",   "", vSelectedLeptons.at(0).pt()    , weight);
    theHistoManager->fillHisto("SubLeadingLeptonPt",               "FullThreeLeptons",   "ttH3l",   "", vSelectedLeptons.at(1).pt()    , weight);
    theHistoManager->fillHisto("ThirdLeptonPt",                    "FullThreeLeptons",   "ttH3l",   "", vSelectedLeptons.at(2).pt()    , weight);
    theHistoManager->fillHisto("MET",                              "FullThreeLeptons",   "ttH3l",   "", vEvent->at(0).metpt()          , weight);
    theHistoManager->fillHisto("MHT",                              "FullThreeLeptons",   "ttH3l",   "", MHT                            , weight);
    theHistoManager->fillHisto("MetLD",                            "FullThreeLeptons",   "ttH3l",   "", met_ld                         , weight);
    theHistoManager->fillHisto("TauMultiplicity",                  "FullThreeLeptons",   "ttH3l",   "", vSelectedTaus.size()           , weight);
    theHistoManager->fillHisto("JetMultiplicity",                  "FullThreeLeptons",   "ttH3l",   "", vSelectedJets.size()           , weight);
    theHistoManager->fillHisto("LooseBJetMultiplicity",            "FullThreeLeptons",   "ttH3l",   "", vSelectedBTagJets.size()       , weight);
    theHistoManager->fillHisto("MediumBJetMultiplicity",           "FullThreeLeptons",   "ttH3l",   "", vSelectedMediumBTagJets.size() , weight);
    theHistoManager->fillHisto("SumOfLeptonsCharges",              "FullThreeLeptons",   "ttH3l",   "", sum_charges                    , weight);
    theHistoManager->fillHisto("SumOfThreeLeptonsCharges",         "FullThreeLeptons",   "ttH3l",   "", sum_charges_3l                 , weight);
    theHistoManager->fillHisto("NumberOfSelectedLeptons",          "FullThreeLeptons",   "ttH3l",   "", vSelectedLeptons.size()        , weight);

    is_3l_TTH_SR = true;   

    // Calcul of input variables of the 2D BDT
    max_Lep_eta     = std::max( abs(vSelectedLeptons.at(0).eta()), abs(vSelectedLeptons.at(1).eta()) ) ;

    numJets_float   = vSelectedJets.size() ;

    mindr_lep1_jet = 1000.;

    for (int i=0; i<vSelectedJets.size(); i++)
    {
        if( DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ) < mindr_lep1_jet )
        { mindr_lep1_jet = DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ); }
    }

    mindr_lep2_jet  = 1000. ;

    for (int i=0; i<vSelectedJets.size(); i++)
    {
        if( DeltaRLeptonJet( vSelectedLeptons.at(1), vSelectedJets.at(i) ) < mindr_lep1_jet )
        { mindr_lep1_jet = DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ); }
    }

    met             = vEvent->at(0).metpt() ;

    //avg_dr_jet      = DeltaRJets( vSelectedNonBTagJets.at(0), vSelectedNonBTagJets.at(1) ) ;
    avg_dr_jet      = DeltaRJets( vSelectedJets.at(0), vSelectedJets.at(1) ) ;

    MT_met_lep1     = sqrt( 2 * vSelectedLeptons.at(0).p4().Pt() * vEvent->at(0).metpt() * (1 - cos( vSelectedLeptons.at(0).phi() - vEvent->at(0).metphi() )));

    LepGood_conePt0 = vSelectedLeptons.at(0).pt() ;

    LepGood_conePt1 = vSelectedLeptons.at(1).pt() ;

    mhtJet25_Recl   = sqrt( (jet_px*jet_px) + (jet_py*jet_py) );

    signal_3l_TT_MVA    = mva_3l_tt->EvaluateMVA("BDTG method");
    signal_3l_TTV_MVA   = mva_3l_ttV->EvaluateMVA("BDTG method");

    //std::cout << " signal 3l   TT MVA: "  << signal_3l_TT_MVA
    //          << " signal 3l   TTV MVA: " << signal_3l_TTV_MVA << std::endl;

    theHistoManager->fillHisto("Signal_3l_TT_MVA",                         "FinalCut",   "ttH3l",   "",  signal_3l_TT_MVA,   weight);
    theHistoManager->fillHisto("Signal_3l_TTV_MVA",                        "FinalCut",   "ttH3l",   "",  signal_3l_TTV_MVA,  weight);



    //fillOutputTree();

    if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);
}

void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_CR_WZ(int evt)
{

    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel", "WZZZ_CR",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel", "WZZZ_CR",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection #
    // ####################

    bool nLep             = ( vSelectedLeptons.size()                   >= 2 );
    if(!nLep)             return;

    bool leading_lep_pt   = ( vSelectedLeptons.at(0).pt()               > 20 );
    if(!leading_lep_pt)   return;

    bool following_lep_pt = ( vSelectedLeptons.at(1).pt()               > 10 );
    if(!following_lep_pt) return;

    bool passMll12Gt12    = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M()  > 12);
    if(!passMll12Gt12)    return;

    bool nJets            = ( vSelectedJets.size()                      >= 2 );
    if(!nJets)            return;

    // #################################
    // # Three leptons event selection #
    // #################################

    nLep                  = ( vSelectedLeptons.size()		        >= 3 );
    if(!nLep)             return;

    bool third_lep_pt     = ( vSelectedLeptons.at(2).pt()               > 10 );
    if(!third_lep_pt)     return;

    passMll12Gt12  = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(2).p4() ).M()  > 12
            && ( vSelectedLeptons.at(1).p4() + vSelectedLeptons.at(2).p4() ).M()  > 12 );
    if(!passMll12Gt12)    return;

    theHistoManager->fillHisto("CutFlow",                           "noSel", "WZZZ_CR",   "",   2, weight);
    theHistoManager->fillHisto("CutFlow",                 "PassingTightMVA", "WZZZ_CR",   "",   1, weight);

    // ##########
    // # Z peak #
    // ##########

    bool pass_Zpeak = false;
    float Delta = 10, ZM = -1., Zpt = -1., Lep1Z = -1, Lep2Z = -1, LepW = -1;

    for(int i=0; i<vSelectedLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vSelectedLeptons.size(); j++)
        {
            if (  ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id()                               )
                    && ( fabs( ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M() - 91.188) < Delta ) )
            {
                Lep1Z      = i;
                Lep2Z      = j;
                ZM         = ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M();
                Zpt        = ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).Pt();
                Delta      = fabs(ZM - 91.188);
                pass_Zpeak = true;
            }
        }
    }
    if(!pass_Zpeak)       return;

    if( ( (Lep1Z == 0) && (Lep2Z == 1) ) || ( (Lep1Z == 1) && (Lep2Z == 0) ) ) LepW = 2;
    if( ( (Lep1Z == 0) && (Lep2Z == 2) ) || ( (Lep1Z == 2) && (Lep2Z == 0) ) ) LepW = 1;    
    if( ( (Lep1Z == 1) && (Lep2Z == 2) ) || ( (Lep1Z == 2) && (Lep2Z == 1) ) ) LepW = 0;

    float MTW = 0. ;
    if (LepW >=0) MTW = sqrt( 2 * vSelectedLeptons.at(LepW).p4().Pt() * vEvent->at(0).metpt() * (1 - cos( vSelectedLeptons.at(LepW).phi() - vEvent->at(0).metphi() )));

    theHistoManager->fillHisto("CutFlow",                           "noSel", "WZZZ_CR",   "",   3, weight);
    theHistoManager->fillHisto("CutFlow",                        "PassingZ", "WZZZ_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass",        "PassingZ", "WZZZ_CR",   "",  ZM, weight);

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                    "PassingMETLD", "WZZZ_CR",   "",   1, weight);
    if( vSelectedJets.size() >= 4  )  theHistoManager->fillHisto("CutFlow",                     "PassingJets", "WZZZ_CR",   "",   1, weight); 

    if( met_ld < 0.2 && vSelectedJets.size() < 4 ) return;
    theHistoManager->fillHisto("CutFlow",                 "PassingMETLDorJets", "WZZZ_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass", "PassingMETLDorJets", "WZZZ_CR",   "",  ZM, weight);

    // ########
    // # SFOS #
    // ########

    bool isSFOS = false;
    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        for(int j=0; j<vSelectedLeptons.size(); j++)
        {
            if (  ( i                           != j                            )
                    && ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id() ) )
            { isSFOS = true ;}
        }
    }

    theHistoManager->fillHisto("CutFlow",                 "PassingMETLDorSFOS", "WZZZ_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass", "PassingMETLDorSFOS", "WZZZ_CR",   "",  ZM, weight);

    if(vSelectedJets.size() < 4 && (met_ld < (0.2 + 0.1 * isSFOS)) ) return;

    // ##############
    // # b-jet veto #
    // ##############

    bool nLooseBtag       = ( nLooseBJets                               == 0 );
    bool nMediumBtag      = ( nMediumBJets                              == 0 );
    if(!nLooseBtag || !nMediumBtag)      return;

    theHistoManager->fillHisto("CutFlow",                 "PassingbJetsVeto", "WZZZ_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass", "PassingbJetsVeto", "WZZZ_CR",   "",  ZM, weight);

    int sum_charges_3l = 0;
    for(int i=0; i<3; i++)
    {
        sum_charges_3l = sum_charges_3l + vSelectedLeptons.at(i).charge();
    }
    if( fabs(sum_charges_3l) != 1 ) return;

    // Building arbitrary variables with all leptons...

    TLorentzVector all_lep_invmass_p4;

    int   all_lep_sumofcharges = 0;
    float all_lep_invmass      = 0;
    float all_lep_sumofpt      = 0;

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        all_lep_sumofcharges = all_lep_sumofcharges + vSelectedLeptons.at(i).charge();
        all_lep_invmass_p4   = all_lep_invmass_p4   + vSelectedLeptons.at(i).p4();
    }

    all_lep_invmass = all_lep_invmass_p4.M();
    all_lep_sumofpt = sqrt( (lepton_px*lepton_px) + (lepton_py*lepton_py) );

    // Building arbitrary variables with remaining (after Z peak determination) leptons...

    TLorentzVector rem_lep_invmass_p4;

    int   rem_lep_sumofcharges = 0;
    float rem_lep_invmass      = 0;
    float rem_lep_sumofpt      = 0;
    float rem_lep_px           = 0;
    float rem_lep_py           = 0;


    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        if( i!=Lep1Z && i!=Lep2Z)
        {
            rem_lep_sumofcharges = rem_lep_sumofcharges + vSelectedLeptons.at(i).charge();
            rem_lep_px           = rem_lep_px           + vSelectedLeptons.at(i).p4().Px();
            rem_lep_py           = rem_lep_py           + vSelectedLeptons.at(i).p4().Py();
            rem_lep_invmass_p4   = rem_lep_invmass_p4   + vSelectedLeptons.at(i).p4();
        }
    }

    rem_lep_invmass = rem_lep_invmass_p4.M();
    rem_lep_sumofpt = sqrt( (rem_lep_px*rem_lep_px) + (rem_lep_py*rem_lep_py) );

    theHistoManager->fillHisto("ZCandidateInvariantMass",        "FinalCut", "WZZZ_CR",   "",   ZM                    , weight);
    theHistoManager->fillHisto("ZCandidateTransverseMomentum",   "FinalCut", "WZZZ_CR",   "",   Zpt                   , weight);
    theHistoManager->fillHisto("MET",                            "FinalCut", "WZZZ_CR",   "",   vEvent->at(0).metpt() , weight);
    theHistoManager->fillHisto("MHT",                            "FinalCut", "WZZZ_CR",   "",   MHT                   , weight);
    theHistoManager->fillHisto("MetLD",                          "FinalCut", "WZZZ_CR",   "",   met_ld                , weight);
    theHistoManager->fillHisto("TauMultiplicity",                "FinalCut", "WZZZ_CR",   "",   vSelectedTaus.size()  , weight);
    theHistoManager->fillHisto("JetMultiplicity",                "FinalCut", "WZZZ_CR",   "",   vSelectedJets.size()  , weight);

    theHistoManager->fillHisto("MTW",                            "FinalCut", "WZZZ_CR",   "",   MTW                   , weight);
    theHistoManager->fillHisto("InvariantMassOfSelectedLeptons", "FinalCut", "WZZZ_CR",   "",   all_lep_invmass       , weight);
    theHistoManager->fillHisto("SumOfSelectedLeptonsCharges",    "FinalCut", "WZZZ_CR",   "",   all_lep_sumofcharges  , weight);
    theHistoManager->fillHisto("SumVecPtSelectedLeptons",        "FinalCut", "WZZZ_CR",   "",   all_lep_sumofpt       , weight);

    theHistoManager->fillHisto("InvMassRemainingLepton",         "FinalCut", "WZZZ_CR",   "",   rem_lep_invmass       , weight);
    theHistoManager->fillHisto("SumOfRemainingLeptonsCharges",   "FinalCut", "WZZZ_CR",   "",   rem_lep_sumofcharges  , weight);
    theHistoManager->fillHisto("SumVecPtRemainingLeptons",       "FinalCut", "WZZZ_CR",   "",   rem_lep_sumofpt       , weight);

    theHistoManager->fillHisto2D("InvMassLastLeptonVSZMass",     "FinalCut", "WZZZ_CR",   "",   rem_lep_invmass,    ZM, weight);
    theHistoManager->fillHisto2D("SumPtLepVSZMass",              "FinalCut", "WZZZ_CR",   "",   rem_lep_sumofpt,    ZM, weight);
    theHistoManager->fillHisto2D("METLDVSZMass",                 "FinalCut", "WZZZ_CR",   "",   met_ld,             ZM, weight);

    is_3l_WZ_CR = true;
}

void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_CR_WZrelaxed(int evt)
{

    theHistoManager->fillHisto2D("SelectedLeptonsVsJets",           "noSel", "WZrel_CR",   "",    vLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("SelectedLeptonsVsBJets",          "noSel", "WZrel_CR",   "",    vLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection #
    // ####################

    bool nLep             = ( vLeptons.size()                   >= 2 );
    if(!nLep)             return;

    bool leading_lep_pt   = ( vLeptons.at(0).pt()               > 20 );
    if(!leading_lep_pt)   return;

    bool following_lep_pt = ( vLeptons.at(1).pt()               > 10 );
    if(!following_lep_pt) return;

    bool passMll12Gt12    = ( ( vLeptons.at(0).p4() + vLeptons.at(1).p4() ).M()  > 12);
    if(!passMll12Gt12)    return;

    bool nJets            = ( vSelectedJets.size()                      >= 2 );
    if(!nJets)            return;

    // #################################
    // # Three leptons event selection #
    // #################################

    nLep                  = ( vLeptons.size()				    >= 3 );
    if(!nLep)             return;

    bool third_lep_pt     = ( vLeptons.at(2).pt()               > 10 );
    if(!third_lep_pt)     return;

    passMll12Gt12  = ( ( vLeptons.at(0).p4() + vLeptons.at(2).p4() ).M()  > 12
            && ( vLeptons.at(1).p4() + vLeptons.at(2).p4() ).M()  > 12 );
    if(!passMll12Gt12)    return;

    theHistoManager->fillHisto("CutFlow",                           "noSel", "WZrel_CR",   "",   2, weight);
    theHistoManager->fillHisto("CutFlow",                 "PassingTightMVA", "WZrel_CR",   "",   1, weight);

    // ##########
    // # Z peak #
    // ##########

    bool pass_Zpeak = false;
    float Delta = 10, ZM = -1., Zpt = -1., Lep1Z = -1, Lep2Z = -1;

    for(int i=0; i<vLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vLeptons.size(); j++)
        {
            if (  ( vLeptons.at(i).id() == -vLeptons.at(j).id()                               )
                    && ( fabs( ( vLeptons.at(i).p4() + vLeptons.at(j).p4() ).M() - 91.188) < Delta ) )
            {
                Lep1Z      = i;
                Lep2Z      = j;
                ZM         = ( vLeptons.at(i).p4() + vLeptons.at(j).p4() ).M();
                Zpt        = ( vLeptons.at(i).p4() + vLeptons.at(j).p4() ).Pt();
                Delta      = fabs(ZM - 91.188);
                pass_Zpeak = true;
            }
        }
    }
    if(!pass_Zpeak)       return;

    theHistoManager->fillHisto("CutFlow",                           "noSel", "WZrel_CR",   "",   3, weight);
    theHistoManager->fillHisto("CutFlow",                        "PassingZ", "WZrel_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass",        "PassingZ", "WZrel_CR",   "",  ZM, weight);

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vLeptons.size(); i++)
    {
        lepton_px = lepton_px + vLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                    "PassingMETLD", "WZrel_CR",   "",   1, weight);
    if( vSelectedJets.size() >= 4  )  theHistoManager->fillHisto("CutFlow",                     "PassingJets", "WZrel_CR",   "",   1, weight); 

    if( met_ld < 0.2 && vSelectedJets.size() < 4 ) return;
    theHistoManager->fillHisto("CutFlow",                 "PassingMETLDorJets", "WZrel_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass", "PassingMETLDorJets", "WZrel_CR",   "",  ZM, weight);

    // ########
    // # SFOS #
    // ########

    bool isSFOS = false;
    for(int i=0; i<vLeptons.size(); i++)
    {
        for(int j=0; j<vLeptons.size(); j++)
        {
            if (  ( i                   != j                            )
                    && ( vLeptons.at(i).id() == -vLeptons.at(j).id() ) )
            { isSFOS = true ;}
        }
    }

    theHistoManager->fillHisto("CutFlow",                 "PassingMETLDorSFOS", "WZrel_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass", "PassingMETLDorSFOS", "WZrel_CR",   "",  ZM, weight);

    if(vSelectedJets.size() < 4 && (met_ld < (0.2 + 0.1 * isSFOS)) ) return;

    // ##############
    // # b-jet veto #
    // ##############

    bool nMediumBtag      = ( nMediumBJets                              == 0 );
    if(!nMediumBtag)      return;

    theHistoManager->fillHisto("CutFlow",                 "PassingbJetsVeto", "WZrel_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass", "PassingbJetsVeto", "WZrel_CR",   "",  ZM, weight);

    int sum_charges_3l = 0;
    for(int i=0; i<3; i++)
    {
        sum_charges_3l = sum_charges_3l + vLeptons.at(i).charge();
    }
    if( fabs(sum_charges_3l) != 1 ) return;

    // Building arbitrary variables with all leptons...

    TLorentzVector all_lep_invmass_p4;

    int   all_lep_sumofcharges = 0;
    float all_lep_invmass      = 0;
    float all_lep_sumofpt      = 0;

    for(int i=0; i<vLeptons.size(); i++)
    {
        all_lep_sumofcharges = all_lep_sumofcharges + vLeptons.at(i).charge();
        all_lep_invmass_p4   = all_lep_invmass_p4   + vLeptons.at(i).p4();
    }

    all_lep_invmass = all_lep_invmass_p4.M();
    all_lep_sumofpt = sqrt( (lepton_px*lepton_px) + (lepton_py*lepton_py) );

    // Building arbitrary variables with remaining (after Z peak determination) leptons...

    TLorentzVector rem_lep_invmass_p4;

    int   rem_lep_sumofcharges = 0;
    float rem_lep_invmass      = 0;
    float rem_lep_sumofpt      = 0;
    float rem_lep_px           = 0;
    float rem_lep_py           = 0;


    for(int i=0; i<vLeptons.size(); i++)
    {
        if( i!=Lep1Z && i!=Lep2Z)
        {
            rem_lep_sumofcharges = rem_lep_sumofcharges + vLeptons.at(i).charge();
            rem_lep_px           = rem_lep_px           + vLeptons.at(i).p4().Px();
            rem_lep_py           = rem_lep_py           + vLeptons.at(i).p4().Py();
            rem_lep_invmass_p4   = rem_lep_invmass_p4   + vLeptons.at(i).p4();
        }
    }

    rem_lep_invmass = rem_lep_invmass_p4.M();
    rem_lep_sumofpt = sqrt( (rem_lep_px*rem_lep_px) + (rem_lep_py*rem_lep_py) );

    theHistoManager->fillHisto("ZCandidateInvariantMass",        "FinalCut", "WZrel_CR",   "",   ZM                    , weight);
    theHistoManager->fillHisto("ZCandidateTransverseMomentum",   "FinalCut", "WZrel_CR",   "",   Zpt                   , weight);
    theHistoManager->fillHisto("MET",                            "FinalCut", "WZrel_CR",   "",   vEvent->at(0).metpt() , weight);
    theHistoManager->fillHisto("MHT",                            "FinalCut", "WZrel_CR",   "",   MHT                   , weight);
    theHistoManager->fillHisto("MetLD",                          "FinalCut", "WZrel_CR",   "",   met_ld                , weight);
    theHistoManager->fillHisto("TauMultiplicity",                "FinalCut", "WZrel_CR",   "",   vSelectedTaus.size()  , weight);
    theHistoManager->fillHisto("JetMultiplicity",                "FinalCut", "WZrel_CR",   "",   vSelectedJets.size()  , weight);

    theHistoManager->fillHisto("InvariantMassOfSelectedLeptons", "FinalCut", "WZrel_CR",   "",   all_lep_invmass       , weight);
    theHistoManager->fillHisto("SumOfSelectedLeptonsCharges",    "FinalCut", "WZrel_CR",   "",   all_lep_sumofcharges  , weight);
    theHistoManager->fillHisto("SumVecPtSelectedLeptons",        "FinalCut", "WZrel_CR",   "",   all_lep_sumofpt       , weight);

    theHistoManager->fillHisto("InvMassRemainingLepton",         "FinalCut", "WZrel_CR",   "",   rem_lep_invmass       , weight);
    theHistoManager->fillHisto("SumOfRemainingLeptonsCharges",   "FinalCut", "WZrel_CR",   "",   rem_lep_sumofcharges  , weight);
    theHistoManager->fillHisto("SumVecPtRemainingLeptons",       "FinalCut", "WZrel_CR",   "",   rem_lep_sumofpt       , weight);

    theHistoManager->fillHisto2D("InvMassLastLeptonVSZMass",     "FinalCut", "WZrel_CR",   "",   rem_lep_invmass,    ZM, weight);
    theHistoManager->fillHisto2D("SumPtLepVSZMass",              "FinalCut", "WZrel_CR",   "",   rem_lep_sumofpt,    ZM, weight);
    theHistoManager->fillHisto2D("METLDVSZMass",                 "FinalCut", "WZrel_CR",   "",   met_ld,             ZM, weight);

    is_3l_WZrel_CR = true;
}

void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_TTZ(int evt)
{

    theHistoManager->fillHisto2D("LeptonsVsJets",           "noSel", "TTZ_CR",   "",    vSelectedLeptons.size(),   vSelectedJets.size(),      weight);
    theHistoManager->fillHisto2D("LeptonsVsBJets",          "noSel", "TTZ_CR",   "",    vSelectedLeptons.size(),   vSelectedBTagJets.size(),  weight);

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)

    // ####################
    // # Common selection # with tight leptons !
    // ####################

    bool nLep             = ( vSelectedLeptons.size()                   >= 2 );
    if(!nLep)             return;

    bool leading_lep_pt   = ( vSelectedLeptons.at(0).pt()               > 20 );
    if(!leading_lep_pt)   return;

    bool following_lep_pt = ( vSelectedLeptons.at(1).pt()               > 10 );
    if(!following_lep_pt) return;

    bool passMll12Gt12    = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M()  > 12);
    if(!passMll12Gt12)    return;

    bool nJets            = ( vSelectedJets.size()                      >= 2 );
    if(!nJets)            return;

    // #################################
    // # Three leptons event selection #
    // #################################

    nLep                  = ( vSelectedLeptons.size()			>= 3 );
    if(!nLep)             return;

    bool third_lep_pt     = ( vSelectedLeptons.at(2).pt()               > 10 );
    if(!third_lep_pt)     return;

    passMll12Gt12  = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(2).p4() ).M()  > 12
            && ( vSelectedLeptons.at(1).p4() + vSelectedLeptons.at(2).p4() ).M()  > 12 );
    if(!passMll12Gt12)    return;

    theHistoManager->fillHisto("CutFlow",                           "noSel", "TTZ_CR",   "",   2, weight);
    theHistoManager->fillHisto("CutFlow",                 "PassingTightMVA", "TTZ_CR",   "",   1, weight);

    // ##########
    // # Z peak #
    // ##########

    bool pass_Zpeak = false;
    float Delta = 10, ZM = -1., Zpt = -1., Lep1Z = -1, Lep2Z = -1;

    for(int i=0; i<vSelectedLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vSelectedLeptons.size(); j++)
        {
            if (  ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id()                               )
                    && ( fabs( ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M() - 91.188) < Delta ) )
            {
                Lep1Z      = i;
                Lep2Z      = j;
                ZM         = ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M();
                Zpt        = ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).Pt();
                Delta      = fabs(ZM - 91.188);
                pass_Zpeak = true;
            }
        }
    }
    if(!pass_Zpeak)       return;

    theHistoManager->fillHisto("CutFlow",                           "noSel", "TTZ_CR",   "",   3, weight);
    theHistoManager->fillHisto("CutFlow",                        "PassingZ", "TTZ_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass",        "PassingZ", "TTZ_CR",   "",  ZM, weight);

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if( met_ld               > 0.2 )  theHistoManager->fillHisto("CutFlow",                    "PassingMETLD", "TTZ_CR",   "",   1, weight);
    if( vSelectedJets.size() >= 4  )  theHistoManager->fillHisto("CutFlow",                     "PassingJets", "TTZ_CR",   "",   1, weight); 

    if( met_ld < 0.2 && vSelectedJets.size() < 4 ) return;
    theHistoManager->fillHisto("CutFlow",                 "PassingMETLDorJets", "TTZ_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass", "PassingMETLDorJets", "TTZ_CR",   "",  ZM, weight);

    // ########
    // # SFOS #
    // ########

    bool isSFOS = false;
    for(int i=0; i<vSelectedLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vSelectedLeptons.size(); j++)
        {
            if (  ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id() ) )
            { isSFOS = true ;}
        }
    }

    theHistoManager->fillHisto("CutFlow",                 "PassingMETLDorSFOS", "TTZ_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass", "PassingMETLDorSFOS", "TTZ_CR",   "",  ZM, weight);

    if(vSelectedJets.size() < 4 && (met_ld < (0.2 + 0.1 * isSFOS)) ) return;

    // #####################
    // # b-jet requirement #
    // #####################

    bool nLooseBtag       = ( nLooseBJets                               >= 2 );
    bool nMediumBtag      = ( nMediumBJets                              >= 1 );
    if(!nLooseBtag || !nMediumBtag)      return;

    theHistoManager->fillHisto("CutFlow",                 "PassingbJetsVeto", "TTZ_CR",   "",   1, weight);
    theHistoManager->fillHisto("ZCandidateInvariantMass", "PassingbJetsVeto", "TTZ_CR",   "",  ZM, weight);

    int sum_charges_3l = 0;
    for(int i=0; i<3; i++)
    {
        sum_charges_3l = sum_charges_3l + vSelectedLeptons.at(i).charge();
    }
    if( fabs(sum_charges_3l) != 1 ) return;

    // Building arbitrary variables with all leptons...

    TLorentzVector all_lep_invmass_p4;

    int   all_lep_sumofcharges = 0;
    float all_lep_invmass      = 0;
    float all_lep_sumofpt      = 0;

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        all_lep_sumofcharges = all_lep_sumofcharges + vSelectedLeptons.at(i).charge();
        all_lep_invmass_p4   = all_lep_invmass_p4   + vSelectedLeptons.at(i).p4();
    }

    all_lep_invmass = all_lep_invmass_p4.M();
    all_lep_sumofpt = sqrt( (lepton_px*lepton_px) + (lepton_py*lepton_py) );

    // Building arbitrary variables with remaining (after Z peak determination) leptons...

    TLorentzVector rem_lep_invmass_p4;

    int   rem_lep_sumofcharges = 0;
    float rem_lep_invmass      = 0;
    float rem_lep_sumofpt      = 0;
    float rem_lep_px           = 0;
    float rem_lep_py           = 0;


    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        if( i!=Lep1Z && i!=Lep2Z)
        {
            rem_lep_sumofcharges = rem_lep_sumofcharges + vSelectedLeptons.at(i).charge();
            rem_lep_px           = rem_lep_px           + vSelectedLeptons.at(i).p4().Px();
            rem_lep_py           = rem_lep_py           + vSelectedLeptons.at(i).p4().Py();
            rem_lep_invmass_p4   = rem_lep_invmass_p4   + vSelectedLeptons.at(i).p4();
        }
    }

    rem_lep_invmass = rem_lep_invmass_p4.M();
    rem_lep_sumofpt = sqrt( (rem_lep_px*rem_lep_px) + (rem_lep_py*rem_lep_py) );

    theHistoManager->fillHisto("LeadingLeptonPt",                "FinalCut", "TTZ_CR",   "",   vSelectedLeptons.at(0).pt()   , weight);
    theHistoManager->fillHisto("SubLeadingLeptonPt",             "FinalCut", "TTZ_CR",   "",   vSelectedLeptons.at(1).pt()   , weight);    

    theHistoManager->fillHisto("ZCandidateInvariantMass",        "FinalCut", "TTZ_CR",   "",   ZM                    , weight);
    theHistoManager->fillHisto("ZCandidateTransverseMomentum",   "FinalCut", "TTZ_CR",   "",   Zpt                   , weight);
    theHistoManager->fillHisto("MET",                            "FinalCut", "TTZ_CR",   "",   vEvent->at(0).metpt() , weight);
    theHistoManager->fillHisto("MHT",                            "FinalCut", "TTZ_CR",   "",   MHT                   , weight);
    theHistoManager->fillHisto("MetLD",                          "FinalCut", "TTZ_CR",   "",   met_ld                , weight);
    theHistoManager->fillHisto("TauMultiplicity",                "FinalCut", "TTZ_CR",   "",   vSelectedTaus.size()  , weight);
    theHistoManager->fillHisto("JetMultiplicity",                "FinalCut", "TTZ_CR",   "",   vSelectedJets.size()  , weight);

    theHistoManager->fillHisto("InvariantMassOfSelectedLeptons", "FinalCut", "TTZ_CR",   "",   all_lep_invmass       , weight);
    theHistoManager->fillHisto("SumOfSelectedLeptonsCharges",    "FinalCut", "TTZ_CR",   "",   all_lep_sumofcharges  , weight);
    theHistoManager->fillHisto("SumVecPtSelectedLeptons",        "FinalCut", "TTZ_CR",   "",   all_lep_sumofpt       , weight);

    theHistoManager->fillHisto("InvMassRemainingLepton",         "FinalCut", "TTZ_CR",   "",   rem_lep_invmass       , weight);
    theHistoManager->fillHisto("SumOfRemainingLeptonsCharges",   "FinalCut", "TTZ_CR",   "",   rem_lep_sumofcharges  , weight);
    theHistoManager->fillHisto("SumVecPtRemainingLeptons",       "FinalCut", "TTZ_CR",   "",   rem_lep_sumofpt       , weight);

    theHistoManager->fillHisto2D("InvMassLastLeptonVSZMass",     "FinalCut", "TTZ_CR",   "",   rem_lep_invmass,    ZM, weight);
    theHistoManager->fillHisto2D("SumPtLepVSZMass",              "FinalCut", "TTZ_CR",   "",   rem_lep_sumofpt,    ZM, weight);
    theHistoManager->fillHisto2D("METLDVSZMass",                 "FinalCut", "TTZ_CR",   "",   met_ld,             ZM, weight);

    if(vSelectedJets.size() >= 4)
    {
        theHistoManager->fillHisto("LeadingLeptonPt",                "FinalCut4j", "TTZ_CR",   "",   vSelectedLeptons.at(0).pt()  , weight);
        theHistoManager->fillHisto("SubLeadingLeptonPt",             "FinalCut4j", "TTZ_CR",   "",   vSelectedLeptons.at(1).pt()  , weight);

        theHistoManager->fillHisto("MET",                            "FinalCut4j", "TTZ_CR",   "",   vEvent->at(0).metpt() , weight);
        theHistoManager->fillHisto("JetMultiplicity",                "FinalCut4j", "TTZ_CR",   "",   vSelectedJets.size()  , weight);
        theHistoManager->fillHisto("ZCandidateInvariantMass",        "FinalCut4j", "TTZ_CR",   "",   ZM                    , weight);
    }

    is_3l_TTZ_CR = true;
}

void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_CR_Zl(int evt)
{
    bool nLep     = ( vSelectedLeptons.size() == 2 ); 
    bool nLepFake = ( vFakeLeptons.size() == 1 );
    bool nJets    = ( vSelectedNonBTagJets.size() + vSelectedBTagJets.size() >= 2 );

    float MZ = -1.;

    bool pass_OSSF = false;
    if (nLep && vSelectedLeptons.at(0).id() == -vSelectedLeptons.at(1).id() &&  fabs( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M() - 91.188 ) < 10  )
    {
        MZ = ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M();
        pass_OSSF = true;
    }

    // common selection
    bool leading_lep_pt = 0;
    if (nLep) leading_lep_pt = ( vSelectedLeptons.at(0).pt() > 20 );
    bool following_lep_pt = 0;
    if (nLep) following_lep_pt = ( vSelectedLeptons.at(1).pt() > 10 );

    bool passMll12Gt12 = 0;
    if (nLep) passMll12Gt12  = ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M()  > 12 );

    bool nLooseBtag  = ( nLooseBJets  >= 2 ); 

    /* std::cout << "nLep: "    << nLep
       << " nLepFake: "	     << nLepFake 
       << " nJets: "	     << nJets 
       << " pass_OSSF: "	     << pass_OSSF 
       << " leading_lep_pt: "   << leading_lep_pt 
       << " following_lep_pt: " << following_lep_pt
       << " passMll12Gt12: "    << passMll12Gt12 
       << " nLooseBtag: "       << nLooseBtag << std::endl; */


    if ( nLep && nLepFake && nJets && pass_OSSF && leading_lep_pt && following_lep_pt && passMll12Gt12 && nLooseBtag )
    {        
        is_Zl_CR = true;
        theHistoManager->fillHisto("ZCandidateInvariantMass", "CR_Zl", "", "", MZ, weight*weight_PV*mc_weight);
        //fillOutputTree();
    }

}

bool TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_TTH3l_MC() 
{ 
    bool sel_MC = true;

    //Check decays 
    if (!((vTruth->at(0).ttbar_decay()==1 && vTruth->at(0).boson_decay()==0) || //ttH
                (vTruth->at(0).ttbar_decay()==2 && vTruth->at(0).boson_decay()==1) || //ttH
                (vTruth->at(0).ttbar_decay()==1 && vTruth->at(0).boson_decay()==2) || //tt semi-lep, ttZ
                (vTruth->at(0).ttbar_decay()==2 && vTruth->at(0).boson_decay()==3)    //tt di-lep, ttW
         )) 
    { 
        sel_MC = false; 
        return sel_MC;}

        // 
        if (vTruth->at(0).JetsHighestPt_phi().size()<2 || vTruth->at(0).JetsClosestMw_phi().size()<2 || vTruth->at(0).JetsLowestMjj_phi().size()<2) 
        { 
            sel_MC = false; 
            return sel_MC;}

            //pt, eta of leptons	
            if (!(vTruth->at(0).Leptons_pt().at(0)> 10 && 
                        vTruth->at(0).Leptons_pt().at(1)> 10 &&
                        vTruth->at(0).Leptons_pt().at(2)> 10   )) sel_MC = false; 


            if (!(fabs(vTruth->at(0).Leptons_eta().at(0)) <2.5 && 
                        fabs(vTruth->at(0).Leptons_eta().at(1)) <2.5 &&
                        fabs(vTruth->at(0).Leptons_eta().at(2)) <2.5   )) sel_MC = false; 


            //lead. lepton
            if (!(vTruth->at(0).Leptons_pt().at(0) > 20 || 
                        vTruth->at(0).Leptons_pt().at(1) > 20 ||
                        vTruth->at(0).Leptons_pt().at(2) > 20  )) sel_MC = false; 


            //SFOS && M(ll) not in 81-101 ??? 
            int SFOSpair = -1;
            if ((vTruth->at(0).Leptons_id().at(0)==-vTruth->at(0).Leptons_id().at(1))) SFOSpair = 0;
            if ((vTruth->at(0).Leptons_id().at(1)==-vTruth->at(0).Leptons_id().at(2))) SFOSpair = 1;
            if ((vTruth->at(0).Leptons_id().at(0)==-vTruth->at(0).Leptons_id().at(2))) SFOSpair = 2;


            TLorentzVector Lep1;
            Lep1.SetPtEtaPhiE(vTruth->at(0).Leptons_pt().at(0),  vTruth->at(0).Leptons_eta().at(0), vTruth->at(0).Leptons_phi().at(0), vTruth->at(0).Leptons_E().at(0));
            TLorentzVector Lep2;
            Lep2.SetPtEtaPhiE(vTruth->at(0).Leptons_pt().at(1),  vTruth->at(0).Leptons_eta().at(1), vTruth->at(0).Leptons_phi().at(1), vTruth->at(0).Leptons_E().at(1));
            TLorentzVector Lep3;
            Lep3.SetPtEtaPhiE(vTruth->at(0).Leptons_pt().at(2),  vTruth->at(0).Leptons_eta().at(2), vTruth->at(0).Leptons_phi().at(2), vTruth->at(0).Leptons_E().at(2));


            if ( !(( Lep1+Lep2 ).M()  > 12 && ( Lep1+Lep3 ).M()  > 12 && ( Lep2+Lep3 ).M()  > 12 )) sel_MC = false;

            if ( (SFOSpair == 0 && fabs( (Lep1+Lep2 ).M()-91.188 ) < 10. ) || 
                    (SFOSpair == 1 && fabs( (Lep2+Lep3 ).M()-91.188 ) < 10. ) ||  
                    (SFOSpair == 2 && fabs( (Lep1+Lep3 ).M()-91.188 ) < 10. )    ) sel_MC = false;

            return sel_MC;

}

float TTbarHiggsMultileptonAnalysis::Phi_0_2Pi(float phi)
{
    float phi_0_2pi = phi;
    if (phi>= TMath::Pi()) phi_0_2pi  -= 2.*TMath::Pi();
    if (phi<0.)            phi_0_2pi  += 2.*TMath::Pi();
    return phi_0_2pi;
}

float TTbarHiggsMultileptonAnalysis::GetDeltaR(float eta1,float phi1,float eta2,float phi2)
{
    float DeltaPhi = TMath::Abs(phi2 - phi1);
    if (DeltaPhi > 3.141593 ) DeltaPhi -= 2.*3.141593;
    return TMath::Sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}
