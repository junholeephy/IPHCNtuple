#include "../include/TTbarHiggsTriggerAnalysis.h"
#include "TSystem.h"

TTbarHiggsTriggerAnalysis::TTbarHiggsTriggerAnalysis() 
{
 
}


TTbarHiggsTriggerAnalysis::TTbarHiggsTriggerAnalysis(TString inputFileName, TChain *tree, TString sampleName, TString treeName, TString outputFileName, bool isdata, float xsec, float lumi, int nowe, int nmax)
{    
   
   //
   _isdata = isdata;
   _xsec = xsec;
   _lumi = lumi;
   _nowe = nowe;
   _nmax = nmax;
   _outputFileName = outputFileName;
   _sampleName = sampleName;
   
   
   tree = new TChain(treeName.Data());

   std::ifstream infile;
   infile.open(inputFileName.Data());
   std::string ifile = "";
   while( getline(infile, ifile) )
     {
	std::string fnameStr = std::string(ifile);
	
	tree->Add(fnameStr.c_str());
	
	std::cout << "file: " << fnameStr << std::endl;
     }   
   infile.close();
   
   Init(tree);

   theHistoManager = new HistoManager();
   
   
   TString outputFileNameRoot = _outputFileName+".root";
   _outputFile = new TFile(outputFileNameRoot.Data(), "recreate");  
  
}

void TTbarHiggsTriggerAnalysis::createHistograms()
{    
 
    _outputFile->cd();

    theHistoManager->addHisto("ProbeMuonPt", "", "", "", 30, 0, 40);
    theHistoManager->addHisto("ProbeMuonPt_IsoMu20_Matched", "", "", "", 30, 0, 40);
    
    theHistoManager->addHisto("ProbeElPt", "", "", "", 30, 0, 40);
    theHistoManager->addHisto("ProbeElPt_Ele23_Matched", "", "", "", 30, 0, 40);
    
    theHistoManager->addHisto("ProbeElPtLeg1", "", "", "", 30, 0, 40);
    theHistoManager->addHisto("ProbeElPtLeg1_Matched", "", "", "", 30, 0, 40);

    theHistoManager->addHisto("ProbeElPtLeg2", "", "", "", 30, 0, 40);
    theHistoManager->addHisto("ProbeElPtLeg2_Matched", "", "", "", 30, 0, 40);

    theHistoManager->addHisto("ProbeMuonPt_Mu17Mu8Leg1", "", "", "", 30, 0, 40);
    theHistoManager->addHisto("ProbeMuonPt_Mu17Mu8Leg1_Matched", "", "", "", 30, 0, 40);
    
    theHistoManager->addHisto("ProbeMuonPt_Mu17Mu8Leg2", "", "", "", 30, 0, 40);
    theHistoManager->addHisto("ProbeMuonPt_Mu17Mu8Leg2_Matched", "", "", "", 30, 0, 40);
    
    theHistoManager->addHisto("ProbeMuonPt_Mu17TkMu8Leg1", "", "", "", 30, 0, 40);
    theHistoManager->addHisto("ProbeMuonPt_Mu17TkMu8Leg1_Matched", "", "", "", 30, 0, 40);
    
    theHistoManager->addHisto("ProbeMuonPt_Mu17TkMu8Leg2", "", "", "", 30, 0, 40);
    theHistoManager->addHisto("ProbeMuonPt_Mu17TkMu8Leg2_Matched", "", "", "", 30, 0, 40);

}  

void TTbarHiggsTriggerAnalysis::writeHistograms()
{  
    _outputFile->cd();
        
    TH1F* eff_IsoMu20 = (TH1F*)theHistoManager->getHisto1D("ProbeMuonPt_IsoMu20_Matched", "", "", "")->Clone("Eff_IsoMu20");
    eff_IsoMu20->SetTitle("Eff_IsoMu20");
    eff_IsoMu20->Divide(theHistoManager->getHisto1D("ProbeMuonPt", "", "", ""));
    theHistoManager->addHisto1D(eff_IsoMu20);

    TH1F* eff_Ele23 = (TH1F*)theHistoManager->getHisto1D("ProbeElPt_Ele23_Matched", "", "", "")->Clone("Eff_Ele23");
    eff_Ele23->SetTitle("Eff_Ele23");
    eff_Ele23->Divide(theHistoManager->getHisto1D("ProbeElPt", "", "", ""));
    theHistoManager->addHisto1D(eff_Ele23);

    TH1F* eff_EleLeg2 = (TH1F*)theHistoManager->getHisto1D("ProbeElPtLeg2_Matched", "", "", "")->Clone("Eff_EleLeg2");
    eff_EleLeg2->SetTitle("Eff_EleLeg2");
    eff_EleLeg2->Divide(theHistoManager->getHisto1D("ProbeElPtLeg2", "", "", ""));
    theHistoManager->addHisto1D(eff_EleLeg2);
    
    TH1F* eff_EleLeg1 = (TH1F*)theHistoManager->getHisto1D("ProbeElPtLeg1_Matched", "", "", "")->Clone("Eff_EleLeg1");
    eff_EleLeg1->SetTitle("Eff_EleLeg1");
    eff_EleLeg1->Divide(theHistoManager->getHisto1D("ProbeElPtLeg1", "", "", ""));
    theHistoManager->addHisto1D(eff_EleLeg1);
   
    TH1F* eff_MuLeg1 = (TH1F*)theHistoManager->getHisto1D("ProbeMuonPt_Mu17Mu8Leg1_Matched", "", "", "")->Clone("Eff_MuLeg1");
    eff_MuLeg1->SetTitle("Eff_MuLeg1");
    eff_MuLeg1->Divide(theHistoManager->getHisto1D("ProbeMuonPt_Mu17Mu8Leg1", "", "", ""));
    theHistoManager->addHisto1D(eff_MuLeg1);
    
    TH1F* eff_MuLeg2 = (TH1F*)theHistoManager->getHisto1D("ProbeMuonPt_Mu17Mu8Leg2_Matched", "", "", "")->Clone("Eff_MuLeg2");
    eff_MuLeg2->SetTitle("Eff_MuLeg2");
    eff_MuLeg2->Divide(theHistoManager->getHisto1D("ProbeMuonPt_Mu17Mu8Leg2", "", "", ""));
    theHistoManager->addHisto1D(eff_MuLeg2);
    
    TH1F* eff_TkMuLeg1 = (TH1F*)theHistoManager->getHisto1D("ProbeMuonPt_Mu17TkMu8Leg1_Matched", "", "", "")->Clone("Eff_TkMuLeg1");
    eff_TkMuLeg1->SetTitle("Eff_TkMuLeg1");
    eff_TkMuLeg1->Divide(theHistoManager->getHisto1D("ProbeMuonPt_Mu17TkMu8Leg1", "", "", ""));
    theHistoManager->addHisto1D(eff_TkMuLeg1);
    
    TH1F* eff_TkMuLeg2 = (TH1F*)theHistoManager->getHisto1D("ProbeMuonPt_Mu17TkMu8Leg2_Matched", "", "", "")->Clone("Eff_TkMuLeg2");
    eff_TkMuLeg2->SetTitle("Eff_TkMuLeg2");
    eff_TkMuLeg2->Divide(theHistoManager->getHisto1D("ProbeMuonPt_Mu17TkMu8Leg2", "", "", ""));
    theHistoManager->addHisto1D(eff_TkMuLeg2);
    
    //
  
    std::vector<TH1F*> the1DHisto =  theHistoManager->getHisto1D_list();
    std::vector<TH2F*> the2DHisto =  theHistoManager->getHisto2D_list();

    for(unsigned int i=0; i<the1DHisto.size(); i++)  the1DHisto[i]->Write();
    for(unsigned int i=0; i<the2DHisto.size(); i++)  the2DHisto[i]->Write();

    _outputFile->Close();
    
    std::cout <<"Dz filter efficiency di_mu :    " <<  _n_Mu17_Mu8/_n_Mu17_Mu8_noDz        <<"+- " << _n_Mu17_Mu8   <<" "<<_n_Mu17_Mu8_noDz<<" "<< sqrt((_n_Mu17_Mu8/_n_Mu17_Mu8_noDz)*(1.-_n_Mu17_Mu8/_n_Mu17_Mu8_noDz))/ sqrt(_n_Mu17_Mu8_noDz)<<  std::endl;
    std::cout <<"Dz filter efficiency di_mu_tk : " <<  _n_Mu17_TkMu8/_n_Mu17_TkMu8_noDz    <<"+- " << _n_Mu17_TkMu8 <<" "<<_n_Mu17_TkMu8_noDz<<" "<< sqrt((_n_Mu17_TkMu8/_n_Mu17_TkMu8_noDz)*(1.-_n_Mu17_TkMu8/_n_Mu17_TkMu8_noDz))/sqrt(_n_Mu17_TkMu8_noDz)<<std::endl;
    std::cout <<"Dz filter efficiency di_el    : " <<  _n_Ele17_Ele12/_n_Ele17_Ele12_noDz  <<"+- " << _n_Ele17_Ele12<<" "<< _n_Ele17_Ele12_noDz<<" "<< sqrt((_n_Mu17_TkMu8/_n_Mu17_TkMu8_noDz)*(1.-_n_Mu17_TkMu8/_n_Mu17_TkMu8_noDz))/sqrt(_n_Ele17_Ele12_noDz)<<std::endl;
    
}


void TTbarHiggsTriggerAnalysis::Init(TChain *tree)
{
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;

    vElectron   = new std::vector<Electron>();
    vMuon       = new std::vector<Muon>();
    vEvent      = new std::vector<Event>();
    vJet        = new std::vector<Jet>();
    vTruth      = new std::vector<Truth>();
    vTriggerObj = new std::vector<TriggerObj>();

    fChain->SetBranchAddress("Electron",   &vElectron);
    fChain->SetBranchAddress("Muon",       &vMuon);
    fChain->SetBranchAddress("Jet",        &vJet);
    fChain->SetBranchAddress("Event",      &vEvent);
    fChain->SetBranchAddress("Truth",      &vTruth);
    fChain->SetBranchAddress("TriggerObj", &vTriggerObj); 
    
    _n_Mu17_Mu8_noDz= 0;
    _n_Mu17_Mu8= 0;
    _n_Mu17_TkMu8_noDz= 0;
    _n_Mu17_TkMu8 = 0;
    _n_Ele17_Ele12_noDz = 0; 
    _n_Ele17_Ele12 = 0;
  
}

void TTbarHiggsTriggerAnalysis::Loop()
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
        nb = fChain->GetEntry(jentry);   nbytes += nb;
	
	Base base;

	////////////////////////////////////////////////
	//
	////////////////////////////////////////////////
	
        vSelectedMuons.clear();
	vSelectedElectrons.clear();
	vSelectedLeptons.clear();
	
	vSelectedTriggerObj_IsoMu20.clear();
	vSelectedTriggerObj_IsoMu20_recoMatched.clear();
	
	vSelectedTriggerObj_Ele23.clear();
	vSelectedTriggerObj_Ele23_recoMatched.clear();
        
	vSelectedTriggerObj_Ele17Ele12Leg1.clear();
        vSelectedTriggerObj_Ele17Ele12Leg2.clear();
	
	vSelectedTriggerObj_Mu17Mu8Leg1.clear();
	vSelectedTriggerObj_Mu17Mu8Leg2.clear();
	vSelectedTriggerObj_Mu17TkMu8Leg1.clear();
	vSelectedTriggerObj_Mu17TkMu8Leg2.clear();
	
	////////////////////////////////////////////////
	//
	////////////////////////////////////////////////

	  
	for(unsigned int imuon=0; imuon < vMuon->size() ; imuon++)
        {             
	    if ( vMuon->at(imuon).lepMVA() > 0.65 && vMuon->at(imuon).isMedium() == true )
	    { 
	      Lepton l; l.setLepton(&vMuon->at(imuon),imuon,0,1);
              vSelectedMuons.push_back(vMuon->at(imuon));    
	      vSelectedLeptons.push_back(l);	    
	     }
	}
	
	
        for(unsigned int ielectron=0; ielectron < vElectron->size() ; ielectron++)
        {   	 
	    if ( vElectron->at(ielectron).lepMVA() > 0.65 )
            { 
	      Lepton l; l.setLepton(&vElectron->at(ielectron),ielectron,1,0);
              vSelectedElectrons.push_back(vElectron->at(ielectron));
	      vSelectedLeptons.push_back(l);	     
              }
	}
	
		
	////////////////////////////////////////////////
	//
	////////////////////////////////////////////////
	
	bool isOSSmu = false;
        bool isOSSel = false;
      
	bool nMu = ( vSelectedMuons.size() == 2 );
	bool nEl = ( vSelectedElectrons.size() == 2 );
	
	if (nMu ) 
        { 
         if ( ( ( vSelectedLeptons.at(0).id() == -vSelectedLeptons.at(1).id() ) && ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M() - 91.188 < 10 ) ) )
         isOSSmu = true ; }
	
	if (nEl ) 
        { 
         if ( ( ( vSelectedLeptons.at(0).id() == -vSelectedLeptons.at(1).id() ) && ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M() - 91.188 < 10 ) ) )
         isOSSel = true ; }
	
	
	////////////////////////////////////////////////
	//
	////////////////////////////////////////////////
	
	int  trigCode = 0;//vEvent->at(0).ev_trigger_pass_byname_1();
	int  trigCode_d = trigCode /10 % 10;
	int  trigCode_c = trigCode /100 % 10;
	int  trigCode_m = trigCode /1000 % 10;
	
	int  trigCode_noDz = -1;//vEvent->at(0).ev_trigger_pass_byname_1_noDz();
	
	////////////////////////////////////////////////
	// MUONS
	////////////////////////////////////////////////
	     
	if (isOSSmu)
	{  
	   //std::cout <<"isOSSmu " << isOSSmu <<" "<< vTriggerObj->size() <<" " << trigCode << std::endl;
	   
	   for (unsigned int i=0 ; i<vTriggerObj->size() ; i++)
           {  
	   
              std::vector<std::string> pathNamesAll = vTriggerObj->at(i).pathNamesAll();
	      std::vector<std::string> filterLabels = vTriggerObj->at(i).filterLabels();
	      
	      //--------------------------------
	      
	      for(int j=0 ; j<vTriggerObj->at(i).pathNamesAll_n() ; j++)
              { 
	        //
	        //Single muon triggers
		//	   	
           	std::size_t ok1 = pathNamesAll.at(j).find("HLT_IsoMu20_v");
		std::size_t ok2 = pathNamesAll.at(j).find("HLT_IsoTkMu20_v");
		
		//std::cout << pathNamesAll.at(j) << std::endl;
		
	   	if ( (ok1!=std::string::npos || ok2!=std::string::npos) && (trigCode_c == 1 || trigCode_c == 2 || trigCode_c == 3 ||
		                                                            trigCode_c == 6 || trigCode_c == 7 || trigCode_c == 8 ) )
		{  		
		  vSelectedTriggerObj_IsoMu20.push_back(vTriggerObj->at(i));
		  
		 /* 
		  if ( base.GetDeltaR(vTriggerObj->at(i).eta(),vTriggerObj->at(i).phi(),vSelectedMuons.at(0).eta(),vSelectedMuons.at(0).phi()) < 0.3 )
		  vSelectedTriggerObj_IsoMu20_recoMatched.push_back(0);
		  else if ( base.GetDeltaR(vTriggerObj->at(i).eta(),vTriggerObj->at(i).phi(),vSelectedMuons.at(1).eta(),vSelectedMuons.at(1).phi()) < 0.3 )
		  vSelectedTriggerObj_IsoMu20_recoMatched.push_back(1);
		  else vSelectedTriggerObj_IsoMu20_recoMatched.push_back(-1);
		  */
		  break;
		}
	      }   
   	   
           //--------------------------------	   
	     bool ok3 = false;
	     bool ok4 = false;
	     bool ok5 = false;
	     bool ok6 = false;
	     bool ok7 = false;
	     bool ok8 = false;
	   
	     for(int k=0 ; k<vTriggerObj->at(i).filterLabels().size() ; k++)
             {
	        //std::cout << filterLabels.at(k) << std::endl;
		
		if ( filterLabels.at(k).find("hltL3fL1sDoubleMu103p5L1f0L2f10L3Filtered17") != std::string::npos ) 
		{ ok3 = true;
		  //std::cout <<"ok 3 true  "<< std::endl;	
		  }						 
		if ( filterLabels.at(k).find("hltDiMuonGlbFiltered17TrkFiltered8") != std::string::npos )
		{  ok4 = true;
		   //std::cout <<"ok 4 true "<< std::endl;	
		   }
		if ( filterLabels.at(k).find("hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4") != std::string::npos ) 
		{  ok5 = true;
		   //std::cout <<"ok 5 true "<< std::endl;	
		   }
		
		if ( filterLabels.at(k).find("hltL3fL1sDoubleMu103p5L1f0L2f10OneMuL3Filtered17") != std::string::npos ) ok6 = true;
		if ( filterLabels.at(k).find("hltL3pfL1sDoubleMu103p5L1f0L2pf0L3PreFiltered8") != std::string::npos ) ok7 = true;
		if ( filterLabels.at(k).find("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4") != std::string::npos ) ok8 = true;
		                                        
	     }
	     //std::cout <<"OK : "<< ok3<<" "<<ok4<<" "<<ok5<<" "<<ok6<<" " <<ok7<<" "<< ok8 << std::endl;
	     if ( ok3 )
	     //{
	       vSelectedTriggerObj_Mu17TkMu8Leg1.push_back(vTriggerObj->at(i));
	       //std::cout <<"yoouppi "<< std::endl;}
	     if ( ok4 && ok5 )
	       vSelectedTriggerObj_Mu17TkMu8Leg2.push_back(vTriggerObj->at(i)); 	               
	     if ( ok6 && ok8 )
	       vSelectedTriggerObj_Mu17Mu8Leg1.push_back(vTriggerObj->at(i));
	     if ( ok7 && ok8 )
	       vSelectedTriggerObj_Mu17Mu8Leg2.push_back(vTriggerObj->at(i));		     
           
	   }

	  
    	   //--------------------------------	   
     	
	   int mu0_IsoMu20Matched = -1;
	   int mu1_IsoMu20Matched = -1; 
	   	    
	   for (unsigned int i=0; i < vSelectedTriggerObj_IsoMu20.size() ; i++)
	   {
	     if (base.GetDeltaR(vSelectedTriggerObj_IsoMu20.at(i).eta(),vSelectedTriggerObj_IsoMu20.at(i).phi(),vSelectedMuons.at(0).eta(),vSelectedMuons.at(0).phi()) < 0.3) 
	       mu0_IsoMu20Matched = i;
	     if (base.GetDeltaR(vSelectedTriggerObj_IsoMu20.at(i).eta(),vSelectedTriggerObj_IsoMu20.at(i).phi(),vSelectedMuons.at(1).eta(),vSelectedMuons.at(1).phi()) < 0.3)
	       mu1_IsoMu20Matched = i;
	   }
	   
	   
	   //--------------------------------	   

	   int mu0_Mu17TkMu8Leg1Matched = -1;
	   int mu1_Mu17TkMu8Leg1Matched = -1;
	   int mu0_Mu17TkMu8Leg2Matched = -1;
	   int mu1_Mu17TkMu8Leg2Matched = -1;
	   
	   int mu0_Mu17Mu8Leg1Matched = -1;
	   int mu1_Mu17Mu8Leg1Matched = -1;
	   int mu0_Mu17Mu8Leg2Matched = -1;
	   int mu1_Mu17Mu8Leg2Matched = -1;
	
	   /*std::cout <<"o1 "<< vSelectedTriggerObj_Mu17Mu8Leg1.size()   << std::endl;
	   std::cout <<"o2 "<< vSelectedTriggerObj_Mu17Mu8Leg2.size()   << std::endl;
	   std::cout <<"o3 "<< vSelectedTriggerObj_Mu17TkMu8Leg1.size() << std::endl;
	   std::cout <<"o4 "<< vSelectedTriggerObj_Mu17TkMu8Leg2.size() << std::endl;*/
	
	
	   for (unsigned int i=0; i < vSelectedTriggerObj_Mu17Mu8Leg1.size() ; i++)
	   {
	     if (base.GetDeltaR(vSelectedTriggerObj_Mu17Mu8Leg1.at(i).eta(),vSelectedTriggerObj_Mu17Mu8Leg1.at(i).phi(),vSelectedMuons.at(0).eta(),vSelectedMuons.at(0).phi()) < 0.3) 
	       mu0_Mu17Mu8Leg1Matched = i;
	     if (base.GetDeltaR(vSelectedTriggerObj_Mu17Mu8Leg1.at(i).eta(),vSelectedTriggerObj_Mu17Mu8Leg1.at(i).phi(),vSelectedMuons.at(1).eta(),vSelectedMuons.at(1).phi()) < 0.3)
	       mu1_Mu17Mu8Leg1Matched = i;
	   }
	   for (unsigned int i=0; i < vSelectedTriggerObj_Mu17Mu8Leg2.size() ; i++)
	   {
	     if (base.GetDeltaR(vSelectedTriggerObj_Mu17Mu8Leg2.at(i).eta(),vSelectedTriggerObj_Mu17Mu8Leg2.at(i).phi(),vSelectedMuons.at(0).eta(),vSelectedMuons.at(0).phi()) < 0.3) 
	       mu0_Mu17Mu8Leg2Matched = i;
	     if (base.GetDeltaR(vSelectedTriggerObj_Mu17Mu8Leg2.at(i).eta(),vSelectedTriggerObj_Mu17Mu8Leg2.at(i).phi(),vSelectedMuons.at(1).eta(),vSelectedMuons.at(1).phi()) < 0.3)
	       mu1_Mu17Mu8Leg2Matched = i;
	   }
	   for (unsigned int i=0; i < vSelectedTriggerObj_Mu17TkMu8Leg1.size() ; i++)
	   {
	     if (base.GetDeltaR(vSelectedTriggerObj_Mu17TkMu8Leg1.at(i).eta(),vSelectedTriggerObj_Mu17TkMu8Leg1.at(i).phi(),vSelectedMuons.at(0).eta(),vSelectedMuons.at(0).phi()) < 0.3) 
	       mu0_Mu17TkMu8Leg1Matched = i;
	     if (base.GetDeltaR(vSelectedTriggerObj_Mu17TkMu8Leg1.at(i).eta(),vSelectedTriggerObj_Mu17TkMu8Leg1.at(i).phi(),vSelectedMuons.at(1).eta(),vSelectedMuons.at(1).phi()) < 0.3)
	       mu1_Mu17TkMu8Leg1Matched = i;
	   }	  	   
	   for (unsigned int i=0; i < vSelectedTriggerObj_Mu17TkMu8Leg2.size() ; i++)
	   {
	     if (base.GetDeltaR(vSelectedTriggerObj_Mu17TkMu8Leg2.at(i).eta(),vSelectedTriggerObj_Mu17TkMu8Leg2.at(i).phi(),vSelectedMuons.at(0).eta(),vSelectedMuons.at(0).phi()) < 0.3) 
	       mu0_Mu17TkMu8Leg2Matched = i;
	     if (base.GetDeltaR(vSelectedTriggerObj_Mu17TkMu8Leg2.at(i).eta(),vSelectedTriggerObj_Mu17TkMu8Leg2.at(i).phi(),vSelectedMuons.at(1).eta(),vSelectedMuons.at(1).phi()) < 0.3)
	       mu1_Mu17TkMu8Leg2Matched = i;
	   }	  
		  
	  
	   //--------------------------------
	   
	   if ( mu0_IsoMu20Matched !=-1 )
	   {
	     theHistoManager->fillHisto("ProbeMuonPt", "", "", "", vSelectedMuons.at(1).pt(),1);
	     if ( mu1_IsoMu20Matched !=-1 ) theHistoManager->fillHisto("ProbeMuonPt_IsoMu20_Matched", "", "", "", vSelectedMuons.at(1).pt(),1);   
	     
	     //
	     //if ( mu0_Mu17TkMu8Leg1Matched !=-1 ) 
	     //{ 
	       theHistoManager->fillHisto("ProbeMuonPt_Mu17TkMu8Leg2", "", "", "", vSelectedMuons.at(1).pt(),1);   
	       if ( mu1_Mu17TkMu8Leg2Matched !=-1 ) theHistoManager->fillHisto("ProbeMuonPt_Mu17TkMu8Leg2_Matched", "", "", "", vSelectedMuons.at(1).pt(),1);   
	      //}
	     
	     theHistoManager->fillHisto("ProbeMuonPt_Mu17TkMu8Leg1", "", "", "", vSelectedMuons.at(1).pt(),1);   
	     if ( mu1_Mu17TkMu8Leg1Matched !=-1 ) theHistoManager->fillHisto("ProbeMuonPt_Mu17TkMu8Leg1_Matched", "", "", "", vSelectedMuons.at(1).pt(),1);   
	     
	     //
	     //if ( mu0_Mu17Mu8Leg1Matched !=-1 ) 
	     //{ 
	       theHistoManager->fillHisto("ProbeMuonPt_Mu17Mu8Leg2", "", "", "", vSelectedMuons.at(1).pt(),1);   
	       if ( mu1_Mu17Mu8Leg2Matched !=-1 ) theHistoManager->fillHisto("ProbeMuonPt_Mu17Mu8Leg2_Matched", "", "", "", vSelectedMuons.at(1).pt(),1);   
	      //}
	     
	     theHistoManager->fillHisto("ProbeMuonPt_Mu17Mu8Leg1", "", "", "", vSelectedMuons.at(1).pt(),1);   
	     if ( mu1_Mu17Mu8Leg1Matched !=-1 ) theHistoManager->fillHisto("ProbeMuonPt_Mu17Mu8Leg1_Matched", "", "", "", vSelectedMuons.at(1).pt(),1);   
	   
	   }
	   
	   
	   //--------------------------------	   

	   if ( trigCode_c == 1 || trigCode_c == 2 || trigCode_c == 3 || trigCode_c == 6 || trigCode_c == 7 || trigCode_c == 8 ) 
	   {
	     if ( trigCode_noDz == 20 || trigCode_noDz == 30 || trigCode_noDz == 70 || trigCode_noDz == 80 )
	     {  
	       _n_Mu17_Mu8_noDz++;
	       if ( trigCode_d == 2 || trigCode_d == 3 || trigCode_d == 7 || trigCode_d == 8 )  _n_Mu17_Mu8++;
	       }
	     
	     if ( trigCode_noDz == 50 || trigCode_noDz == 60 || trigCode_noDz == 70 || trigCode_noDz == 80 )
	     { 
	       _n_Mu17_TkMu8_noDz++;
	       if ( trigCode_d == 5 || trigCode_d == 6 || trigCode_d == 7 || trigCode_d == 8 )  _n_Mu17_TkMu8++;
	       }   
	    }
	 
        }//isOSSmu
		
	////////////////////////////////////////////////
	//
	////////////////////////////////////////////////
	
	if (isOSSel)
	{ 
	   std::cout <<"isOSSel "<< isOSSel <<" "<< vTriggerObj->size() <<" "<< trigCode << std::endl;
	
	   for (unsigned int i=0; i < vTriggerObj->size() ; i++)
           {  
	    	
              std::vector<std::string> pathNamesAll = vTriggerObj->at(i).pathNamesAll();
	      std::vector<std::string> filterLabels = vTriggerObj->at(i).filterLabels();
	           
             
	      for(int j=0;j<vTriggerObj->at(i).pathNamesAll_n();j++)
              { 	   	
           	//---------------------------------
		std::size_t ok1 = pathNamesAll.at(j).find("HLT_Ele23_WPLoose_Gsf_v");
		std::size_t ok2 = pathNamesAll.at(j).find("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v");
		   
	   	//if ( ( ok1!=std::string::npos || ok2!=std::string::npos ) && ( trigCode_c >= 5 && trigCode_c <= 8 ) ) 
		if ( ( ok1!=std::string::npos || ok2!=std::string::npos ) ) 		
		{  
		
		  vSelectedTriggerObj_Ele23.push_back(vTriggerObj->at(i));
		  break;
		  	  
		  /*
		  if ( base.GetDeltaR(vTriggerObj->at(i).eta(),vTriggerObj->at(i).phi(),vSelectedElectrons.at(0).eta(),vSelectedElectrons.at(0).phi()) < 0.3 )
		  vSelectedTriggerObj_Ele23_recoMatched.push_back(0);
		  else if ( base.GetDeltaR(vTriggerObj->at(i).eta(),vTriggerObj->at(i).phi(),vSelectedElectrons.at(1).eta(),vSelectedElectrons.at(1).phi()) < 0.3 )
		  vSelectedTriggerObj_Ele23_recoMatched.push_back(1);
		  else vSelectedTriggerObj_Ele23_recoMatched.push_back(-1);*/
		  
		 }
	       }
	      
	      
	       for(int k=0;k<vTriggerObj->at(i).filterLabels().size();k++)
               {
	       
		std::size_t ok3 = filterLabels.at(k).find("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter");
		std::size_t ok4 = filterLabels.at(k).find("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter");
		
		if ( ok3 != std::string::npos )
		  { vSelectedTriggerObj_Ele17Ele12Leg1.push_back(vTriggerObj->at(i)); }
		if ( ok4 != std::string::npos )
		  { vSelectedTriggerObj_Ele17Ele12Leg2.push_back(vTriggerObj->at(i)); }   
               }
	     
            }

	   //---------------------------------
	   int el0_Ele23Matched = -1;
	   int el1_Ele23Matched = -1;
	   for (unsigned int i=0; i < vSelectedTriggerObj_Ele23.size() ; i++)
	   {
	     if (base.GetDeltaR(vSelectedTriggerObj_Ele23.at(i).eta(),vSelectedTriggerObj_Ele23.at(i).phi(),vSelectedElectrons.at(0).eta(),vSelectedElectrons.at(0).phi()) < 0.3) 
	       el0_Ele23Matched = i;
	     if (base.GetDeltaR(vSelectedTriggerObj_Ele23.at(i).eta(),vSelectedTriggerObj_Ele23.at(i).phi(),vSelectedElectrons.at(1).eta(),vSelectedElectrons.at(1).phi()) < 0.3)
	       el1_Ele23Matched = i;
	   }
	   
	   
	   //---------------------------------
	   int el0_Ele17Ele12Leg1Matched = -1;
	   int el1_Ele17Ele12Leg1Matched = -1;
	   int el0_Ele17Ele12Leg2Matched = -1;
	   int el1_Ele17Ele12Leg2Matched = -1;
	   
	   for (unsigned int i=0; i < vSelectedTriggerObj_Ele17Ele12Leg1.size() ; i++)
	   {
	     if (base.GetDeltaR(vSelectedTriggerObj_Ele17Ele12Leg1.at(i).eta(),vSelectedTriggerObj_Ele17Ele12Leg1.at(i).phi(),vSelectedElectrons.at(0).eta(),vSelectedElectrons.at(0).phi()) < 0.3) 
	       el0_Ele17Ele12Leg1Matched = i;
	     if (base.GetDeltaR(vSelectedTriggerObj_Ele17Ele12Leg1.at(i).eta(),vSelectedTriggerObj_Ele17Ele12Leg1.at(i).phi(),vSelectedElectrons.at(1).eta(),vSelectedElectrons.at(1).phi()) < 0.3)
	       el1_Ele17Ele12Leg1Matched = i;
	   }
	   for (unsigned int i=0; i < vSelectedTriggerObj_Ele17Ele12Leg2.size() ; i++)
	   {
	     if (base.GetDeltaR(vSelectedTriggerObj_Ele17Ele12Leg2.at(i).eta(),vSelectedTriggerObj_Ele17Ele12Leg2.at(i).phi(),vSelectedElectrons.at(0).eta(),vSelectedElectrons.at(0).phi()) < 0.3) 
	       el0_Ele17Ele12Leg2Matched = i;
	     if (base.GetDeltaR(vSelectedTriggerObj_Ele17Ele12Leg2.at(i).eta(),vSelectedTriggerObj_Ele17Ele12Leg2.at(i).phi(),vSelectedElectrons.at(1).eta(),vSelectedElectrons.at(1).phi()) < 0.3)
	       el1_Ele17Ele12Leg2Matched = i;
	   } 
	   
	   //---------------------------------
	   if ( el0_Ele23Matched !=-1 )
	   {
	     theHistoManager->fillHisto("ProbeElPt", "", "", "", vSelectedElectrons.at(1).pt(),1);
	     if ( el1_Ele23Matched !=-1 )
	     theHistoManager->fillHisto("ProbeElPt_Ele23_Matched", "", "", "", vSelectedElectrons.at(1).pt(),1);  
	     
	     if ( el0_Ele17Ele12Leg1Matched !=-1 )
	     {
	       theHistoManager->fillHisto("ProbeElPtLeg2", "", "", "", vSelectedElectrons.at(1).pt(),1);
	       if ( el1_Ele17Ele12Leg2Matched !=-1 ) theHistoManager->fillHisto("ProbeElPtLeg2_Matched", "", "", "", vSelectedElectrons.at(1).pt(),1);  
              }
	      
	    theHistoManager->fillHisto("ProbeElPtLeg1", "", "", "", vSelectedElectrons.at(1).pt(),1);
	    if ( el1_Ele17Ele12Leg1Matched !=-1 )
	    theHistoManager->fillHisto("ProbeElPtLeg1_Matched", "", "", "", vSelectedElectrons.at(1).pt(),1);  
	      
	    }
	     	 
	   //---------------------------------
	   if ( trigCode_c >= 5 && trigCode_c <= 8 ) 
	   {
	     if ( trigCode_noDz == 10 || trigCode_noDz == 30 || trigCode_noDz == 60 || trigCode_noDz == 80 )
	     {
	      _n_Ele17_Ele12_noDz++;
	      if ( trigCode_d == 1 || trigCode_d == 3 || trigCode_d == 6 || trigCode_d == 8 )  _n_Ele17_Ele12++;
	     }
	   }
	 
        }//isOSSel
		
		     
    }
    

}
