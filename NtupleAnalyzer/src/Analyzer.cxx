#include "../include/Hist.h"
#include "../include/TTbarDileptonAnalysis.h"
#include "../include/TTbarHiggsMultileptonAnalysis.h"
#include "../include/TTbarHiggsTFAnalysis.h"

#include "TChain.h"

#include <iostream>

int main(int argc, char *argv[])
{
   if( argc < 3 )                                                                                                          
     {
	std::cout << "NtupleAnalyzer usage:"        << std::endl;
	std::cout << "--file: input filename"       << std::endl;
	std::cout << "--tree: TTree name"           << std::endl;
	std::cout << "--nmax: max number of events" << std::endl;
	exit(1);
     }

   const char *fname_str  = "list.txt";
   const char *stream_str = "Nt";
   int nmax = -1;

   for(int i=0;i<argc;i++)
     {
	if( ! strcmp(argv[i],"--file") ) fname_str  = argv[i+1];
	if( ! strcmp(argv[i],"--tree") ) stream_str = argv[i+1];
	if( ! strcmp(argv[i],"--nmax") ) nmax = atoi(argv[i+1]);
     }
   
   const char *fname  = fname_str;                                                                                         
   const char *stream = stream_str;                                                                                        
   
   std::cout << "--file=" << fname  << std::endl;                                                                          
   std::cout << "--tree=" << stream << std::endl;                                                                          
   std::cout << "--nmax=" << nmax   << std::endl;

   TChain *thetree = 0;
   
   //TTH MEM analysis
      TTbarHiggsMultileptonAnalysis *TTHanalysis = new TTbarHiggsMultileptonAnalysis(fname,thetree,"TTbarHiggs",stream);
   // TTHanalysis->InitLHCO(1,1); // to print LHCO files
      TTHanalysis->createHistograms();
      TTHanalysis->Loop();
      TTHanalysis->writeHistograms();

   ////TTH Transfer Function analysis
   //TTbarHiggsTFAnalysis *TTHTFanalysis = new TTbarHiggsTFAnalysis(fname,thetree,"TTbarHiggs",stream);
   //TTHTFanalysis->createHistograms();
   //TTHTFanalysis->Loop();
   //TTHTFanalysis->writeHistograms();



}
