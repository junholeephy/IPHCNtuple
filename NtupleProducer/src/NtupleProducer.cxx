#include "../include/NtupleProducer.h"

#include "Riostream.h"
#include "TSystem.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <assert.h>

Tree                 *ntP;
TChain                *ch;
Ntuple                *nt;
std::vector<int> *evdebug;

unsigned int          idx;

int main(int argc, char *argv[])
{
    if( argc < 5 )
    {
        std::cout << "NtupleProducer usage:"        << std::endl;
        std::cout << "--file: input filename"       << std::endl;
        std::cout << "--outfile: output filename"   << std::endl;
        std::cout << "--tree: TTree name"           << std::endl;
        std::cout << "--nmax: max number of events" << std::endl;
        std::cout << "--isdata: data or MC ?"       << std::endl;
        exit(1);
    }

    const char *fname_str     = "input.txt";
    const char *fname_out_str = "output.root";
    const char *stream_str    = "FlatTree/tree";
    int        nmax           = -1;
    bool       isdata         = 0;

    for(int i=0;i<argc;i++)
    {
        if( ! strcmp(argv[i],"--file") )    fname_str      = argv[i+1];
	if( ! strcmp(argv[i],"--outfile") ) fname_out_str  = argv[i+1];	
        if( ! strcmp(argv[i],"--tree") )    stream_str     = argv[i+1];
        if( ! strcmp(argv[i],"--nmax") )    nmax           = atoi(argv[i+1]);
	if( ! strcmp(argv[i],"--isdata") )  isdata         = (bool) atoi(argv[i+1]);
	
    }   

    const char *fname = fname_str;
    const char *stream = stream_str;
    const char *fname_out = fname_out_str;

    std::cout << "--file="   << fname      << std::endl;
    std::cout << "--outfile="<< fname_out  << std::endl;
    std::cout << "--tree="   << stream     << std::endl;
    std::cout << "--nmax="   << nmax       << std::endl;
    std::cout << "--isdata=" << isdata     << std::endl;

    Tree tree(0,const_cast<char*>(fname),stream);
    ntP = &tree;

    ch = tree.fChain;
    Long64_t nentries = ch->GetEntries();
    ntP->registerInputBranches(ch);

    nt = new Ntuple(fname_out);

    nt->Init();
    std::cout << "Initialization DONE" << std::endl;
    nt->createVar();
    std::cout << "Creation DONE" << std::endl;
    nt->setBranchAddress();
    std::cout << "SetBranchAddress DONE" << std::endl;

    Event     ev;
    Electron  el;
    Muon      mu;
    Tau      tau;
    Jet      jet;
    Truth  truth;
    GenJet genjet;
    TriggerObj trigObj;
    
    evdebug = new std::vector<int>();
    //   evdebug->push_back(120);

    int nlep = 0;
    int njet = 0;
    
   
    for(Long64_t i=0;i<nentries;i++)
    {
        if( i > nmax && nmax >= 0 ) break; 

        //std::cout << "i:" << i << std::endl;
        ch->GetEntry(i);
   
 
        nt->clearVar();	

        //	if( !(isHtoWW || isHtoZZ || isHtoTT) ) continue;

        // event
        ev.init();

        ev.read(isdata); 

        //	if( isHtoWW ) ev._tth_channel = 0;
        //	else if( isHtoZZ ) ev._tth_channel = 1;
        //	else if( isHtoTT ) ev._tth_channel = 2;

        nt->NtEvent->push_back(ev);
	
         // electrons
        for(int j=0;j<ntP->el_n;j++)
        {
            idx = j;

            el.init();
            el.read();
            if( el.sel() ) nlep++;
           
            if( el.sel()) nt->NtElectron->push_back(el);
        }
   
        // muons
        for(int j=0;j<ntP->mu_n;j++)
        {
            idx = j;

            mu.init();
            mu.read();
         
            if( mu.sel()) nt->NtMuon->push_back(mu);
        }
  
        int x_tau = 0;
        // taus 
        for(int j=0;j<ntP->tau_n;j++)
        {
            idx = j;

            tau.init();
            tau.read();
         
            if (tau.sel()) nt->NtTau->push_back(tau);
        }	

        // jets
        for(int j=0;j<ntP->jet_n;j++)
        {
            idx = j;

            jet.init();
            jet.read(isdata);
            
            if (jet.sel()) nt->NtJet->push_back(jet);
        }
	
	//trigger objects
	for(int j=0;j<ntP->triggerobject_n;j++)
        {
            idx = j;
            
            trigObj.init();
            trigObj.read();
          
            if (trigObj.sel()) nt->NtTriggerObj->push_back(trigObj);
        }
	
	if (!isdata)
	{
	  // genjets
          for(int j=0;j<ntP->genJet_n;j++)
          {
              idx = j;
	      
              genjet.init();
              genjet.read();
	      
              if (genjet.sel()) nt->NtGenJet->push_back(genjet);
          }

          // truth
          truth.init();
          truth.read();
          truth.readMultiLepton();

          nt->NtTruth->push_back(truth);
	}

        nt->fill();
    }  

    delete evdebug;
    delete nt;
}
