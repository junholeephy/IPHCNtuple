#include "../include/ChargeFlip.h"

// fill the histograms (done once)
void fillQFhistos(TFile* fileFR)
{

    h_QF_wgt = (TH2D*)fileFR->Get("chargeMisId");

    //for( int iEta=1; iEta<3; iEta++)
    //{
    //    for( int iPt=1; iPt<4; iPt++ ) std::cout << "h_QF_wgt[iPt][iEta] with [iPt]: "<< iPt << " [iEta]: " << iEta << " weight: " << h_QF_wgt->GetBinContent(iPt,iEta) << std::endl;
    //}
    
    return;
}

double get_QF_wgt_2l( std::vector<double> leptonsPts,
                      std::vector<double> leptonsEtas)
{

    double weight    = 1;
    double weight_QF = 0;

    for( int ilepton=0; ilepton<int(leptonsPts.size()); ilepton++ )
    {

        double leptonPt  = leptonsPts[ilepton];
        double leptonEta = fabs( leptonsEtas[ilepton] );

        int x = h_QF_wgt->GetXaxis()->FindBin(leptonPt);
        int y = h_QF_wgt->GetYaxis()->FindBin(leptonEta);

        weight_QF = weight_QF + h_QF_wgt->GetBinContent(x,y); 
 
        //std::cout << "weight from charge flip : " << weight << std::endl;
    }

    weight    = weight * weight_QF;
    return weight;
}

