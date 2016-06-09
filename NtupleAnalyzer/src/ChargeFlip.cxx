#include "../include/ChargeFlip.h"

// fill the histograms (done once)
void fillQFhistos(TFile* fileFR)
{

    h_QF_wgt = (TH2D*)fileFR->Get("chargeMisId");

    for( int iEta=1; iEta<3; iEta++)
    {
        for( int iPt=1; iPt<4; iPt++ ) std::cout << "h_QF_wgt[iPt][iEta] with [iPt]: "<< iPt << " [iEta]: " << iEta << " weight: " << h_QF_wgt->GetBinContent(iPt,iEta) << std::endl;
    }
    
    return;
}

double get_QF_wgt_2l( std::vector<double> leptonsPts,
                      std::vector<double> leptonsEtas)
{

    double weight = 0;

    for( int ileptons=0; ileptons<int(leptonsPts.size()); ileptons++ )
    {
        weight = 1;
        //
    }

    return weight;
}

