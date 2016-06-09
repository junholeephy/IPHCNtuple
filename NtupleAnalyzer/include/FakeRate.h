#ifndef FAKERATE_H
#define FAKERATE_H


void fillFRhistos(TFile *fileFR);

double get_FR_wgt_2l( std::vector<double> leptonsPts,
                      std::vector<double> leptonsEtas,
                      std::vector<int>    leptonsId );

double get_FR_wgt_3l( std::vector<double> leptonsPts,
                      std::vector<double> leptonsEtas,
                      std::vector<int>    leptonsId );

TH2D* h_FR_wgt_el;
TH2D* h_FR_wgt_mu;

#endif
