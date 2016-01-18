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

    theHistoManager->addHisto("ProbeMuonPt", "", "", _sampleName.Data(), 100, 0, 200);
    theHistoManager->addHisto("ProbeMuonPt_IsoMu20_Matched", "", "", _sampleName.Data(), 100, 0, 200);
    
    theHistoManager->addHisto("ProbeElPt", "", "", _sampleName.Data(), 100, 0, 200);
    theHistoManager->addHisto("ProbeElPt_Ele23_Matched", "", "", _sampleName.Data(), 100, 0, 200);
    
    theHistoManager->addHisto("ProbeElPtLeg2", "", "", _sampleName.Data(), 100, 0, 200);
    theHistoManager->addHisto("ProbeElPtLeg2_Matched", "", "", _sampleName.Data(), 100, 0, 200);
}  

void TTbarHiggsTriggerAnalysis::writeHistograms()
{  
    _outputFile->cd();
    
    //
    
    TH1F* eff_IsoMu20 = theHistoManager->getHisto1D("ProbeMuonPt", "", "", _sampleName.Data());
    eff_IsoMu20->SetTitle("Eff_IsoMu20");
    eff_IsoMu20->Divide(theHistoManager->getHisto1D("ProbeMuonPt_IsoMu20_Matched", "", "", _sampleName.Data()));
    theHistoManager->addHisto1D(eff_IsoMu20);

    TH1F* eff_Ele23 = theHistoManager->getHisto1D("ProbeElPt", "", "", _sampleName.Data());
    eff_Ele23->SetTitle("Eff_Ele23");
    eff_Ele23->Divide(theHistoManager->getHisto1D("ProbeElPt_Ele23_Matched", "", "", _sampleName.Data()));
    theHistoManager->addHisto1D(eff_Ele23);

    TH1F* eff_EleLeg2 = theHistoManager->getHisto1D("ProbeElPtLeg2", "", "", _sampleName.Data());
    eff_EleLeg2->SetTitle("Eff_EleLeg2");
    eff_EleLeg2->Divide(theHistoManager->getHisto1D("ProbeElPtLeg2_Matched", "", "", _sampleName.Data()));
    theHistoManager->addHisto1D(eff_EleLeg2);

    //
    
    std::vector<TH1F*> the1DHisto =  theHistoManager->getHisto1D_list();
    std::vector<TH2F*> the2DHisto =  theHistoManager->getHisto2D_list();

    for(unsigned int i=0; i<the1DHisto.size(); i++)  the1DHisto[i]->Write();
    for(unsigned int i=0; i<the2DHisto.size(); i++)  the2DHisto[i]->Write();

    _outputFile->Close();
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
	

	////////////////////////////////////////////////
	//
	////////////////////////////////////////////////
	
        vSelectedMuons.clear();
	vSelectedElectrons.clear();
	vSelectedLeptons.clear();
	vSelectedTriggerObj_IsoMu20.clear();
	vSelectedTriggerObj_IsoMu20_recoMatched.clear();
	
	  
	for(unsigned int imuon=0; imuon < vMuon->size() ; imuon++)
        {             
	    if ( vMuon->at(imuon).lepMVA() > 0.65 && vMuon->at(imuon).isMediumMuon() == true )
	    { 
	      Lepton l; l.setLepton(&vMuon->at(imuon),imuon,1);
              vSelectedMuons.push_back(vMuon->at(imuon));    
	      vSelectedLeptons.push_back(l);	    
	     }
	}
	
	
        for(unsigned int ielectron=0; ielectron < vElectron->size() ; ielectron++)
        {   	 
	    if ( vElectron->at(ielectron).lepMVA() > 0.65 )
            { 
	      Lepton l; l.setLepton(&vElectron->at(ielectron),ielectron,1);
              vSelectedElectrons.push_back(vElectron->at(ielectron));
	      vSelectedLeptons.push_back(l);	     
              }
	}
	
	
	//std::cout <<"sizes " << vSelectedMuons.size() <<" " <<vSelectedElectrons.size() <<" " <<  vSelectedLeptons.size() << std::endl;
	
	////////////////////////////////////////////////
	//
	////////////////////////////////////////////////
	
	bool isOSSmu = false;
        bool isOSSel = false;
      
	bool nMu = ( vSelectedMuons.size() == 2 );
	bool nEl = ( vSelectedElectrons.size() == 2 );
	
	if (nMu && !nEl ) 
        { 
         if ( ( ( vSelectedLeptons.at(0).id() == -vSelectedLeptons.at(1).id() ) && ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M() - 91.188 < 10 ) ) )
         isOSSmu = true ;
        }
	
	if (nEl && !nMu ) 
        { 
         if ( ( ( vSelectedLeptons.at(0).id() == -vSelectedLeptons.at(1).id() ) && ( ( vSelectedLeptons.at(0).p4() + vSelectedLeptons.at(1).p4() ).M() - 91.188 < 10 ) ) )
         isOSSel = true ;
        }
	
	Base base;
	  
	int  trigCode = vEvent->at(0).ev_trigger_pass_byname_1();
	//std::cout << "trigCode " <<trigCode <<  std::endl;
	int  trigCode_c = trigCode /100 % 10;
	//std::cout << "trigCode_c " << trigCode_c << std::endl;
	int  trigCode_m = trigCode /1000 % 10;
	//std::cout << "trigCode_m " << trigCode_m << std::endl;
	
	////////////////////////////////////////////////
	//
	////////////////////////////////////////////////
	     
	if (isOSSmu)
	{  
	   //std::cout <<"isOSSmu " << isOSSmu <<" "<< vTriggerObj->size() << std::endl;
	   
	   for (unsigned int i=0; i < vTriggerObj->size() ; i++)
           {  std::cout <<"------------"<<i<< std::endl;
	   
              std::vector<std::string> pathNamesAll = vTriggerObj->at(i).pathNamesAll();
	      for(int j=0;j<vTriggerObj->at(i).pathNamesAll_n();j++)
              { 
	        //
	        //Single muon triggers
		//	   	
           	std::size_t ok1 = pathNamesAll.at(j).find("HLT_IsoMu20");
		std::size_t ok2 = pathNamesAll.at(j).find("HLT_IsoTkMu20");
		
		//std::cout << pathNamesAll.at(j) << std::endl;
		
	   	if ( (ok1!=std::string::npos || ok2!=std::string::npos) && (trigCode == 3 || trigCode == 8) )
		{  
		  std::cout << vTriggerObj->at(i).pT()<<" "<<vTriggerObj->at(i).eta() <<" "<<vTriggerObj->at(i).phi()<< std::endl;
		  vSelectedTriggerObj_IsoMu20.push_back(vTriggerObj->at(i));
		  
		  if ( base.GetDeltaR(vTriggerObj->at(i).eta(),vTriggerObj->at(i).phi(),vSelectedMuons.at(0).eta(),vSelectedMuons.at(0).phi()) < 0.3 )
		  vSelectedTriggerObj_IsoMu20_recoMatched.push_back(0);
		  else if ( base.GetDeltaR(vTriggerObj->at(i).eta(),vTriggerObj->at(i).phi(),vSelectedMuons.at(1).eta(),vSelectedMuons.at(1).phi()) < 0.3 )
		  vSelectedTriggerObj_IsoMu20_recoMatched.push_back(1);
		  else vSelectedTriggerObj_IsoMu20_recoMatched.push_back(-1);
		  
		}
		
		//
	        //Di-muon triggers
		//	   	
           	std::size_t ok3 = pathNamesAll.at(j).find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL");	
		
              }   
            }

	
	   int mu0_IsoMu20Matched = -1;
	   int mu1_IsoMu20Matched = -1;
	   
	   //std::cout << vSelectedTriggerObj_IsoMu20.size() << std::endl;
	   
	   for (unsigned int i=0; i < vSelectedTriggerObj_IsoMu20.size() ; i++)
	   {
	     if (base.GetDeltaR(vSelectedTriggerObj_IsoMu20.at(i).eta(),vSelectedTriggerObj_IsoMu20.at(i).phi(),vSelectedMuons.at(0).eta(),vSelectedMuons.at(0).phi()) < 0.3) 
	       mu0_IsoMu20Matched = i;
	     if (base.GetDeltaR(vSelectedTriggerObj_IsoMu20.at(i).eta(),vSelectedTriggerObj_IsoMu20.at(i).phi(),vSelectedMuons.at(1).eta(),vSelectedMuons.at(1).phi()) < 0.3)
	       mu1_IsoMu20Matched = i;
	   }
	    
	   //std::cout <<"IsoMu20Matched " << mu0_IsoMu20Matched << " " << mu1_IsoMu20Matched <<  std::endl;
	    
	   if ( mu0_IsoMu20Matched !=-1 )
	   {
	     std::cout <<"Tag muon found !" << std::endl;
	     theHistoManager->fillHisto("ProbeMuonPt", "", "", "", vSelectedMuons.at(1).pt(),1);
	     if ( mu1_IsoMu20Matched !=-1 ) theHistoManager->fillHisto("ProbeMuonPt_IsoMu20_Matched", "", "", "", vSelectedMuons.at(1).pt(),1);   
	   }
	 
        }//isOSSmu
		
	////////////////////////////////////////////////
	//
	////////////////////////////////////////////////
	
	if (isOSSel)
	{ 
	   std::cout <<"isOSSel " << isOSSel <<" "<< vTriggerObj->size()<< " " << trigCode << std::endl;
	   
	   for (unsigned int i=0; i < vTriggerObj->size() ; i++)
           {  std::cout <<"------------"<<i<< std::endl;
	   std::cout << vTriggerObj->at(i).pT()<<" "<<vTriggerObj->at(i).eta() <<" "<<vTriggerObj->at(i).phi()<< std::endl;
		
              std::vector<std::string> pathNamesAll = vTriggerObj->at(i).pathNamesAll();
	      for(int j=0;j<vTriggerObj->at(i).pathNamesAll_n();j++)
              { 	   	
           	//---------------------------------
		std::size_t ok1 = pathNamesAll.at(j).find("HLT_Ele23_WPLoose_Gsf");
		std::size_t ok2 = pathNamesAll.at(j).find("HLT_Ele23_CaloIdL_TrackIdL_IsoVL");
		
		std::cout <<"ok "<< pathNamesAll.at(j) << ok1 <<" " << ok2 << std::endl;
		
	   	if ( (ok1!=std::string::npos || ok2!=std::string::npos) && trigCode_m >= 1 && ( trigCode_c >= 5 && trigCode_c <= 8 ) ) 
		{  
		
		  vSelectedTriggerObj_Ele23.push_back(vTriggerObj->at(i));
		  
		  std::cout <<" pushing back "<< std::endl;
		  
		  if ( base.GetDeltaR(vTriggerObj->at(i).eta(),vTriggerObj->at(i).phi(),vSelectedElectrons.at(0).eta(),vSelectedElectrons.at(0).phi()) < 0.3 )
		  vSelectedTriggerObj_Ele23_recoMatched.push_back(0);
		  else if ( base.GetDeltaR(vTriggerObj->at(i).eta(),vTriggerObj->at(i).phi(),vSelectedElectrons.at(1).eta(),vSelectedElectrons.at(1).phi()) < 0.3 )
		  vSelectedTriggerObj_Ele23_recoMatched.push_back(1);
		  else vSelectedTriggerObj_Ele23_recoMatched.push_back(-1);
		  
		}
		
		//---------------------------------
		std::size_t ok3 = pathNamesAll.at(j).find("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter");
		std::size_t ok4 = pathNamesAll.at(j).find("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter");
		
		if ( ok3 != std::string::npos )
		  vSelectedTriggerObj_Ele17Ele12Leg1.push_back(vTriggerObj->at(i));
		if ( ok4 != std::string::npos )
		  vSelectedTriggerObj_Ele17Ele12Leg2.push_back(vTriggerObj->at(i));
		
              }   
            }


std::cout << "vSelectedTriggerObj_Ele23.size() " << vSelectedTriggerObj_Ele23.size() << std::endl;

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
	   
	std::cout << "el0_Ele23Matched " << el0_Ele23Matched << " " << el1_Ele23Matched << std::endl;
   
	   
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
	     std::cout <<"Tag electron found !" << std::endl;
	     
	     theHistoManager->fillHisto("ProbeElPt", "", "", "", vSelectedElectrons.at(1).pt(),1);
	     if ( el1_Ele23Matched !=-1 ) theHistoManager->fillHisto("ProbeElPt_Ele23_Matched", "", "", "", vSelectedElectrons.at(1).pt(),1);  
	     
	     if ( el0_Ele17Ele12Leg1Matched != -1 )
	     { 
	        theHistoManager->fillHisto("ProbeElPtLeg2", "", "", "", vSelectedElectrons.at(1).pt(),1);
	        if ( el1_Ele17Ele12Leg2Matched !=-1 ) theHistoManager->fillHisto("ProbeElPtLeg2_Matched", "", "", "", vSelectedElectrons.at(1).pt(),1);  
	     }	   
	   }
	   
	   
	   
	 
        }//isOSSel
		
		     

    }
    

}
