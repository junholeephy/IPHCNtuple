#include "../include/FakeRate.h"

// fill the histograms (done once)
void fillFRhistos(TFile* fileFR)
{

    h_FR_wgt_el = (TH2D*)fileFR->Get("FR_mva075_el_data_comb");
    h_FR_wgt_mu = (TH2D*)fileFR->Get("FR_mva075_mu_data_comb");

    //std::cout << "Fake Rate Electrons ==========" << std::endl;

    //for( int iEta=1; iEta<3; iEta++)
    //{
    //    for( int iPt=1; iPt<6; iPt++ ) std::cout << "h_FR_wgt_el[iPt][iEta] with [iPt]: "<< iPt << " [iEta]: " << iEta << " weight: " << h_FR_wgt_el->GetBinContent(iPt,iEta) << std::endl;
    //}

    //std::cout << "Fake Rate Muons    ==========" << std::endl;

    //for( int iEta=1; iEta<3; iEta++)
    //{
    //    for( int iPt=1; iPt<6; iPt++) std::cout << "h_FR_wgt_mu[iPt][iEta] with [iPt]: "<< iPt << " [iEta]: " << iEta << " weight: " << h_FR_wgt_mu->GetBinContent(iPt,iEta) << std::endl;
    //}

    return;
}

double get_FR_wgt_2l( std::vector<double> leptonsPts,
                      std::vector<double> leptonsEtas,
                      std::vector<int>    leptonsId)
{
    double weight    = 1;
    double weight_FR = 1;

    int iPt, iEta;
    double leptonPt, leptonEta, leptonId;

    for(int i=0; i<leptonsPts.size(); i++)
    {
        double leptonPt  = leptonsPts[i];
        double leptonEta = fabs( leptonsEtas[i] );
        int    leptonId  = abs(  leptonsId[i]   );

        int x = h_FR_wgt_el->GetXaxis()->FindBin(leptonPt);
        int y = h_FR_wgt_el->GetYaxis()->FindBin(leptonEta);

        if(leptonId == 11)
        {
            double f2 = h_FR_wgt_el->GetBinContent(x,y);
            weight_FR = f2 / (1-f2);
            //std::cout << "electron - pt : " << leptonPt << " eta : " << leptonEta << " weigt: " << weight_FR << std::endl;
        }
        else if (leptonId == 13)
        {
            double f2 = h_FR_wgt_mu->GetBinContent(x,y);
            weight_FR = f2 / (1-f2);
            //std::cout << "muon - pt : " << leptonPt << " eta : " << leptonEta << " weight: " << weight_FR << std::endl;
        }

        weight    = weight * weight_FR;

        //std::cout << "weight: " << weight << std::endl;
    }

    return weight;
}

double get_FR_wgt_3l( std::vector<double> leptonsPts,
                      std::vector<double> leptonsEtas,
                      std::vector<int>    leptonsIds)
{
    double weight    = 1;
    double weight_FR = 1;

    int iPt, iEta;
    double leptonPt, leptonEta, leptonId;

    for(int i=0; i<leptonsPts.size(); i++)
    {   
        double leptonPt   = leptonsPts[i];
        double leptonEta  = fabs( leptonsEtas[i] );
        int    leptonId   = abs(  leptonsIds[i]   );

        int x = h_FR_wgt_el->GetXaxis()->FindBin(leptonPt);
        int y = h_FR_wgt_el->GetYaxis()->FindBin(leptonEta);

        if(leptonId == 11)
        {
            double f2 = h_FR_wgt_el->GetBinContent(x,y);
            weight_FR = f2 / (1-f2);
            //std::cout << "electron - pt : " << leptonPt << " eta : " << leptonEta << " weigt: " << weight_FR << std::endl;
        }
        else if (leptonId == 13)
        {
            double f2 = h_FR_wgt_mu->GetBinContent(x,y);
            weight_FR = f2 / (1-f2);
            //std::cout << "muon - pt : " << leptonPt << " eta : " << leptonEta << " weight: " << weight_FR << std::endl;
        }

        weight    = weight * weight_FR;

        //std::cout << "weight: " << weight << std::endl;

    }

    return weight;
}




