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

    double weight;

    for( int ileptons=0; ileptons<int(leptonsPts.size()); ileptons++ )
    {
        weight = 1;
        //
    }

    return weight;
}

