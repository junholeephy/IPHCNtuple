#include "../include/TTbarHiggsTFAnalysis.h"
#include "TSystem.h"

TTbarHiggsTFAnalysis::TTbarHiggsTFAnalysis() 
{
 
}


TTbarHiggsTFAnalysis::TTbarHiggsTFAnalysis(TString inputFileName, TChain *tree, TString theSampleName, TString treeName)
{    
//   gSystem->Load("libNtuple.so");

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
   
   sampleName = theSampleName;

   outputfile = new TFile("output.root", "recreate");
}

void TTbarHiggsTFAnalysis::createHistograms()
{    
 
    outputfile->cd();

    theHistoManager->addHisto("MuonPt", "noSel", "emu", sampleName.Data(), 100, 0, 200);
   
}

void TTbarHiggsTFAnalysis::writeHistograms()
{  
    outputfile->cd();

    std::vector<TH1F*> the1DHisto =  theHistoManager->getHisto1D_list();
    std::vector<TH2F*> the2DHisto =  theHistoManager->getHisto2D_list();

    for(unsigned int i=0; i<the1DHisto.size(); i++)  the1DHisto[i]->Write();
    for(unsigned int i=0; i<the2DHisto.size(); i++)  the2DHisto[i]->Write();

    outputfile->Close();
}


void TTbarHiggsTFAnalysis::Init(TChain *tree)
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

void TTbarHiggsTFAnalysis::Loop()
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

      
        for(unsigned int imuon=0; imuon < vMuon->size() ; imuon++)
        {
           std::cout<< vMuon->at(imuon).lepMVA() << std::endl;

        }     

    }

}
