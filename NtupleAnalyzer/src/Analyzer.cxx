#include "../include/Hist.h"
#include "../include/TTbarDileptonAnalysis.h"
#include "../include/TTbarHiggsMultileptonAnalysis.h"
#include "TTree.h"

#include <iostream>

int main(int argc, char *argv[])
{
    /*if( argc < 3 )
      {
      std::cout << "Usage: ./Analyzer [input file] [tool] [evc]" << std::endl;
      exit(1);
      }   */
    TTree * thetree = 0;

    //  TTbarDileptonAnalysis * thettbaranalysis = new TTbarDileptonAnalysis("../../NtupleProducer/test/output.root", thetree, "TTbar");
    //  thettbaranalysis->createHistograms();
    //  thettbaranalysis->Loop();
    //  thettbaranalysis->writeHistograms();

    //std::cout << "Jusqu'ici, tout va bien 0.0" << std::endl;
    TTbarHiggsMultileptonAnalysis * TTHanalysis = new TTbarHiggsMultileptonAnalysis("../../NtupleProducer/test/output.root", thetree, "TTbarHiggs");
    //std::cout << "Jusqu'ici, tout va bien 0.1" << std::endl;
    TTHanalysis->createHistograms();
    //std::cout << "Jusqu'ici, tout va bien 0.2" << std::endl;
    TTHanalysis->Loop();
    //std::cout << "Jusqu'ici, tout va bien 0.3" << std::endl;
    TTHanalysis->writeHistograms();
    //std::cout << "Jusqu'ici, tout va bien 0.4" << std::endl;

}
