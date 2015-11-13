#include "../include/TTbarDileptonAnalysis.h"


TTbarDileptonAnalysis::TTbarDileptonAnalysis(){


}



TTbarDileptonAnalysis::TTbarDileptonAnalysis(TString inputfilename, TTree *tree, TString theSampleName){
    
    if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(inputfilename.Data());
      if (!f || !f->IsOpen()) {
        f = new TFile(inputfilename.Data());
      }
      f->GetObject("Nt",tree);
    }
    Init(tree);
    
    theHistoManager = new HistoManager();
    
    sampleName = theSampleName;
    
    
}

void TTbarDileptonAnalysis::createHistograms(){
  
  
    outputfile = new TFile("output.root", "recreate");
    
    
   theHistoManager->addHisto("CutFlow",      "noSel", "emu", sampleName.Data(), 10, 0, 10);
    
    
   theHistoManager->addHisto("MuonPt",      "noSel", "emu", sampleName.Data(), 100, 0, 200);
   theHistoManager->addHisto("MuonEta",     "noSel", "emu", sampleName.Data(), 100, -3, 3);
   
   theHistoManager->addHisto("ElectronPt",  "noSel", "emu", sampleName.Data(), 100, 0, 200);
   theHistoManager->addHisto("ElectronEta", "noSel", "emu", sampleName.Data(), 100, -3, 3);
   
   theHistoManager->addHisto("JetPt",       "noSel", "emu", sampleName.Data(), 100, 0, 200);
   theHistoManager->addHisto("JetEta",      "noSel", "emu", sampleName.Data(), 100, -3, 3);
}

void TTbarDileptonAnalysis::writeHistograms(){
  
  outputfile->cd();
  
  
  std::vector<TH1F*> the1DHisto =  theHistoManager->getHisto1D_list();
  std::vector<TH2F*> the2DHisto =  theHistoManager->getHisto2D_list();
  
  for(unsigned int i=0; i<the1DHisto.size(); i++)  the1DHisto[i]->Write();
  for(unsigned int i=0; i<the2DHisto.size(); i++)  the2DHisto[i]->Write();
  
  
  
  
}



void TTbarDileptonAnalysis::Init(TTree *tree)
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

void TTbarDileptonAnalysis::Loop(){

  if (fChain == 0) return;
  
  
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   
  
  
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = fChain->LoadTree(jentry);
      if (ientry < 0) break;
      
      if(jentry%10000 == 0) std::cout << "number of processed events " << jentry << std::endl;
      
      if(jentry > 1000000) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      float theweight = vEvent->at(0).mc_weight();
    
      
      std::vector<Muon>     vSelectedMuons;
      std::vector<Electron> vSelectedElectrons;
      std::vector<Jet>      vSelectedNonBTagJets;
      std::vector<Jet>      vSelectedBTagJets;
      
      //-------------------------------------------
      //selection muons
      //-------------------------------------------
      for(unsigned int imuon=0; imuon < vMuon->size() ; imuon++){
        if( !vMuon->at(imuon).passPtEta() ) continue;
        if(  vMuon->at(imuon).pt() < 20 )   continue;
	if( !vMuon->at(imuon).isLoose() )   continue;
	
	theHistoManager->fillHisto("MuonPt",  "noSel", "emu", sampleName.Data(),  vMuon->at(imuon).pt(), theweight);
	theHistoManager->fillHisto("MuonEta", "noSel", "emu", sampleName.Data(),  vMuon->at(imuon).eta(), theweight);
	
	vSelectedMuons.push_back(vMuon->at(imuon));
	
     }     
     
      //-------------------------------------------
      //selection electrons
      //-------------------------------------------
     for(unsigned int ielectron=0; ielectron < vElectron->size() ; ielectron++){
        if( !vElectron->at(ielectron).passPtEta() ) continue;
        if(  vElectron->at(ielectron).pt() < 20 )   continue;
	if( !vElectron->at(ielectron).isLoose() )   continue;
	
	theHistoManager->fillHisto("ElectronPt",  "noSel", "emu", sampleName.Data(),  vElectron->at(ielectron).pt(), theweight);
	theHistoManager->fillHisto("ElectronEta", "noSel", "emu", sampleName.Data(),  vElectron->at(ielectron).eta(), theweight);
	
	vSelectedElectrons.push_back(vElectron->at(ielectron));
	
     }  
     
       
      //-------------------------------------------
      //selection jets
      //-------------------------------------------
     for(unsigned int ijet=0; ijet < vJet->size() ; ijet++){
        if(  vJet->at(ijet).pt() < 30   )         continue;
        if(  fabs(vJet->at(ijet).eta()) > 2.4 )   continue;
	
	theHistoManager->fillHisto("JetPt",  "noSel", "emu", sampleName.Data(),  vJet->at(ijet).pt(), theweight);
	theHistoManager->fillHisto("JetEta", "noSel", "emu", sampleName.Data(),  vJet->at(ijet).eta(), theweight);
	
	if(vJet->at(ijet).CSVv2() >= 0.244 ) vSelectedBTagJets.push_back(vJet->at(ijet));
	else                                 vSelectedNonBTagJets.push_back(vJet->at(ijet));
	
     }
	
      //-------------------------------------------
      //apply the event selection
      //-------------------------------------------  
      theHistoManager->fillHisto("CutFlow", "noSel", "emu", sampleName.Data(),  0, theweight);
      
      //at least 2 leptons
      if(vSelectedMuons.size() == 1 && vSelectedElectrons.size() == 1){
      
        theHistoManager->fillHisto("CutFlow", "noSel", "emu", sampleName.Data(),  1, theweight);
	
        //at least 2 jets
        if( (vSelectedBTagJets.size() + vSelectedNonBTagJets.size()) >= 2){
          theHistoManager->fillHisto("CutFlow", "noSel", "emu", sampleName.Data(),  2, theweight);
	  
	  if(vSelectedBTagJets.size() >=1){
	  
            theHistoManager->fillHisto("CutFlow", "noSel", "emu", sampleName.Data(),  3, theweight);
	    
	  }//NBjet selection
	  
	  
	}//end NJet selection
      }//end nlepton selection
      
      
   }
}





