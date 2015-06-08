#include "../include/NtupleProducer.h"

#include "Riostream.h"
#include "TSystem.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <assert.h>

Tree *ntP;
TChain *ch;
Ntuple *nt;
std::vector<int> *evdebug;

unsigned int idx;

int main(int argc, char *argv[])
{
   if( argc < 3 )
     {
	std::cout << "NtupleProducer usage:" << std::endl;
	std::cout << "--file: input filename" << std::endl;
	std::cout << "--tree: TTree name" << std::endl;
	std::cout << "--nmax: max number of events" << std::endl;
	exit(1);
     }
   
   const char *fname_str = "output.root";
   const char *stream_str = "FlatTree/tree";
   int nmax = -1;
   
   for(int i=0;i<argc;i++)
     {
	if( ! strcmp(argv[i],"--file") ) fname_str = argv[i+1];
	if( ! strcmp(argv[i],"--tree") ) stream_str = argv[i+1];
	if( ! strcmp(argv[i],"--nmax") ) nmax = atoi(argv[i+1]);
     }   
   
   const char *fname  = fname_str;
   const char *stream = stream_str;
   
   std::cout << "--file=" << fname << std::endl;
   std::cout << "--tree=" << stream << std::endl;
   
   Tree tree(0,const_cast<char*>(fname),stream);
   ntP = &tree;
   
   ch = tree.fChain;
   Long64_t nentries = ch->GetEntries();
   ntP->registerInputBranches(ch);
   
   nt = new Ntuple();

   nt->Init();
   nt->createVar();
   nt->setBranchAddress();

   Electron el;
   Muon mu;
   Jet jet;
   Event ev;
   Truth truth;

   evdebug = new std::vector<int>();
//   evdebug->push_back(120);
   
   for(Long64_t i=0;i<nentries;i++)
     {
	if( i > nmax && nmax >= 0 ) break; 
	
	ch->GetEntry(i);

	nt->clearVar();	
	
/*	bool isHtoWW = (abs(ntP->mc_truth_h0W1_id) == 24 &&
			abs(ntP->mc_truth_h0W2_id) == 24);
	bool isHtoZZ = (abs(ntP->mc_truth_h0Z1_id) == 23 &&
			abs(ntP->mc_truth_h0Z2_id) == 23);
	bool isHtoTT = (abs(ntP->mc_truth_h0tau1_id) == 15 &&
			abs(ntP->mc_truth_h0tau2_id) == 15);	
*/	
//	if( !(isHtoWW || isHtoZZ || isHtoTT) ) continue;
	
	// event
	ev.init();
	ev.read();
	
//	if( isHtoWW ) ev._tth_channel = 0;
//	else if( isHtoZZ ) ev._tth_channel = 1;
//	else if( isHtoTT ) ev._tth_channel = 2;
	
	nt->NtEvent->push_back(ev);
	
	// muons
	for(int j=0;j<ntP->mu_n;j++)
	  {
	     idx = j;
	     
	     mu.init();
	     mu.read();
	     mu.sel();
		  
	     nt->NtMuon->push_back(mu);
	  }					

	// electrons
	for(int j=0;j<ntP->el_n;j++)
	  {
	     idx = j;
	     
	     el.init();
	     el.read();
	     el.sel();
		  
	     nt->NtElectron->push_back(el);
	  }					
	
	// jets
	for(int j=0;j<ntP->jet_n;j++)
	  {
	     idx = j;
	     
	     jet.init();
	     jet.read();
	     jet.sel();
	     
	     nt->NtJet->push_back(jet);
	  }					

	// truth
	truth.init();
	truth.read();
	
	nt->NtTruth->push_back(truth);
	
	nt->fill();
     }   

   delete nt;
}
