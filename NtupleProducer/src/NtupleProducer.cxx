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
    if( argc < 3 )
    {
        std::cout << "NtupleProducer usage:"        << std::endl;
        std::cout << "--file: input filename"       << std::endl;
        std::cout << "--tree: TTree name"           << std::endl;
        std::cout << "--nmax: max number of events" << std::endl;
        exit(1);
    }

    const char *fname_str  = "output.root";
    const char *stream_str = "FlatTree/tree";
    int        nmax        = -1;

    for(int i=0;i<argc;i++)
    {
        if( ! strcmp(argv[i],"--file") ) fname_str  = argv[i+1];
        if( ! strcmp(argv[i],"--tree") ) stream_str = argv[i+1];
        if( ! strcmp(argv[i],"--nmax") ) nmax       = atoi(argv[i+1]);
    }   

    const char *fname  = fname_str;
    const char *stream = stream_str;

    std::cout << "--file=" << fname  << std::endl;
    std::cout << "--tree=" << stream << std::endl;
    std::cout << "--nmax=" << nmax   << std::endl;

    Tree tree(0,const_cast<char*>(fname),stream);
    ntP = &tree;

    ch = tree.fChain;
    Long64_t nentries = ch->GetEntries();
    ntP->registerInputBranches(ch);

    nt = new Ntuple();

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

    evdebug = new std::vector<int>();
    //   evdebug->push_back(120);

    int sync_muon = 0;
    int sync_el   = 0;
    int sync_jet  = 0;
    int sync_tau  = 0;

    for(Long64_t i=0;i<nentries;i++)
    {
        if( i > nmax && nmax >= 0 ) break; 

        //std::cout << "i:" << i << std::endl;
        ch->GetEntry(i);
        nt->clearVar();	

        /*	
         	bool isHtoWW = (abs(ntP->mc_truth_h0W1_id) == 24 &&
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

        int x_mu = 0;
        // muons
        for(int j=0;j<ntP->mu_n;j++)
        {
            idx = j;

            mu.init();
            mu.read();

            if( mu.sel() && (x_mu == 0) )
            {
                sync_muon = sync_muon + 1;
                x_mu = 1;
            }

            //std::cout << "sync_muon: " << sync_muon << std::endl;

            if ( mu.sel() ) nt->NtMuon->push_back(mu);
            //if ( x_mu == 1 ) break;
        }

        //std::cout << "Muons DONE" << std::endl;

        int x_el = 0;
        // electrons
        for(int j=0;j<ntP->el_n;j++)
        {
            idx = j;

            el.init();
            el.read();
            
            if( el.sel() && (x_el == 0) )
            {
                sync_el = sync_el + 1;
                x_el = 1;
            }

            //std::cout << "sync_el: " << sync_el << std::endl;

            if ( el.sel() ) nt->NtElectron->push_back(el);
            //if ( x_el == 1 ) break;
        }

        //std::cout << "Electrons DONE" << std::endl;

        int x_tau = 0;
        // taus 
        for(int j=0;j<ntP->tau_n;j++)
        {
            idx = j;

            tau.init();
            tau.read();
            //tau.sel();

            if( tau.sel() && (x_tau == 0) )
            {
                sync_tau = sync_tau + 1;
                x_tau = 1;
            }

            //std::cout << "sync_tau: " << sync_tau << std::endll;

            if ( tau.sel() ) nt->NtTau->push_back(tau);
            //if ( x_tau == 1 ) break;
        }		

        //std::cout << "Taus DONE" << std::endl;

        int x_jet = 0;
        // jets
        for(int j=0;j<ntP->jet_n;j++)
        {
            idx = j;

            jet.init();
            jet.read();
            //jet.sel();

            //if ( jet.sel() && (x_jet == 0) )
            //{
            //    sync_jet = sync_jet + 1;
            //    x_jet = 1;
            //}

            if( jet.sel() ) nt->NtJet->push_back(jet);
            //if ( x_jet == 1 ) break;
        }

        //std::cout << "Jets DONE" << std::endl;

        // truth
        truth.init();
        //truth.read();

        //nt->NtTruth->push_back(truth);

        //std::cout << "Truth DONE (kinda) " << std::endl;

        nt->fill();

    }   

    delete nt;
}
