#include "../include/FakeRate.h"

// fill the histograms (done once)
void fillFRhistos(TFile* fileFR)
{

    h_FR_wgt_el = (TH2D*)fileFR->Get("FR_mva075_el_data_comb");
    h_FR_wgt_mu = (TH2D*)fileFR->Get("FR_mva075_mu_data_comb");

    std::cout << "Fake Rate Electrons ==========" << std::endl;

    for( int iEta=1; iEta<3; iEta++)
    {
        for( int iPt=1; iPt<6; iPt++ ) std::cout << "h_FR_wgt_el[iPt][iEta] with [iPt]: "<< iPt << " [iEta]: " << iEta << " weight: " << h_FR_wgt_el->GetBinContent(iPt,iEta) << std::endl;
    }

    std::cout << "Fake Rate Muons    ==========" << std::endl;

    for( int iEta=1; iEta<3; iEta++)
    {
        for( int iPt=1; iPt<6; iPt++) std::cout << "h_FR_wgt_mu[iPt][iEta] with [iPt]: "<< iPt << " [iEta]: " << iEta << " weight: " << h_FR_wgt_mu->GetBinContent(iPt,iEta) << std::endl;
    }

    return;
}

double get_FR_wgt_2l( std::vector<double> leptonsPts,
                      std::vector<double> leptonsEtas,
                      std::vector<int>    leptonsId)
{
    double weight    = 1;
    double weight_FR = 1;

    if( leptonsPts.size() == 1)
    {
        int iPt, iEta;
        double leptonPt  = leptonsPts[0];
        double leptonEta = fabs( leptonsEtas[0] );
        int    leptonId  = abs(  leptonsId[0]   );

        if      ( leptonPt < 25 ) iPt = 1;
        else if ( leptonPt < 30 ) iPt = 2;
        else if ( leptonPt < 40 ) iPt = 3;
        else if ( leptonPt < 40 ) iPt = 4;
        else                      iPt = 5;

        if      ( leptonEta < 1.479 ) iEta = 1;
        else if ( leptonEta < 2.5   ) iEta = 2;
        
        if(leptonId == 11)
        {
            double f1 = h_FR_wgt_el->GetBinContent(iPt,iEta);
            weight_FR = f1 / (1-f1);
        }
        else if (leptonId == 13)
        {
            double f1 = h_FR_wgt_mu->GetBinContent(iPt,iEta);
            weight_FR = f1/ (1-f1);
        }

        weight    = weight * weight_FR;

        //std::cout << "weight from charge flip : " << weight << std::endl;
    }
    else
    {
        int iPt0, iEta0, iPt1, iEta1;
        double leptonPt0  = leptonsPts[0];
        double leptonPt1  = leptonsPts[1];
        double leptonEta0 = fabs( leptonsEtas[0] );
        double leptonEta1 = fabs( leptonsEtas[1] );
        int    leptonId0  = abs(  leptonsId[0]   );
        int    leptonId1  = abs(  leptonsId[1]   );




    }

    return weight;
}

double get_FR_wgt_3l( std::vector<double> leptonsPts,
                      std::vector<double> leptonsEtas,
                      std::vector<int>    leptonsId)
{

    double weight;

    for( int ileptons=0; ileptons<int(leptonsPts.size()); ileptons++ )
    {
        weight = 1;
        //
    }

    return weight;
}




