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

    double weight    = 1;
    double weight_QF = 1;

    for( int ilepton=0; ilepton<int(leptonsPts.size()); ilepton++ )
    {

        int iPt, iEta;
        double leptonPt  = leptonsPts[ilepton];
        double leptonEta = fabs( leptonsEtas[ilepton] );

        if      ( leptonPt < 25 ) iPt = 1;
        else if ( leptonPt < 50 ) iPt = 2;
        else                      iPt = 3;

        if      ( leptonEta < 1.479 ) iEta = 1;
        else if ( leptonEta < 2.5   ) iEta = 2;

        weight_QF = h_QF_wgt->GetBinContent(iPt,iEta); 
        weight    = weight * weight_QF;
 
        //std::cout << "weight from charge flip : " << weight << std::endl;

        //
    }

    return weight;
}

