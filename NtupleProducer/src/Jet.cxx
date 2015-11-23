#include "include/NtupleProducer.h"
#include <iostream>
#include <iomanip>

ClassImp(Jet)

Jet::Jet()
{
}

Jet::~Jet()
{
}

void Jet::read(bool isdata)
{
    _ID = idx;

    // general informations
    if( CHECK(ntP->jet_E)     )       _E       = ntP->jet_E->at(idx);
    if( CHECK(ntP->jet_pt)    )       _pt      = ntP->jet_pt->at(idx);
    if( CHECK(ntP->jet_eta)   )       _eta     = ntP->jet_eta->at(idx);   
    if( CHECK(ntP->jet_phi)   )       _phi     = ntP->jet_phi->at(idx);   
    if( CHECK(ntP->jet_m)     )       _m       = ntP->jet_m->at(idx);

    // selection variables
    if( CHECK(ntP->jet_looseJetID) )  _isLoose = ntP->jet_looseJetID->at(idx);

    // other variables
    if( CHECK(ntP->jet_ntrk) )        _ntrk  = ntP->jet_ntrk->at(idx);
    if( CHECK(ntP->jet_CSVv2) )       _CSVv2 = ntP->jet_CSVv2->at(idx);


    // Gen Jets variables

    if (!isdata)
    {
       _jet_partonFlavour    = ntP->jet_partonFlavour->at(idx);
       _jet_hadronFlavour    = ntP->jet_hadronFlavour->at(idx);
       _jet_genJet_pt        = ntP->jet_genJet_pt->at(idx);
       _jet_genParton_pt     = ntP->jet_genParton_pt->at(idx);
       _jet_genParton_id     = ntP->jet_genParton_id ->at(idx);}

}

void Jet::init()
{
    //general information
    _E       = -666;
    _pt      = -666;
    _eta     = -666;
    _phi     = -666;
    _m       = -666;

    // selection variables
    _isLoose = 0;

    // other variables
    _ntrk    = -666;
    _CSVv2   = -666.;
 
    // Gen Jet variables
    _jet_partonFlavour    = -666.;
    _jet_hadronFlavour    = -666.;
    _jet_genJet_pt        = -666.;
    _jet_genParton_pt     = -666.;
    _jet_genParton_id     = -666.;
}

bool Jet::sel()
{
    // selection
    bool pass_pt      = (_pt        > 25.);
    bool pass_eta     = (fabs(_eta) < 2.4);
    bool pass_isLoose = (_isLoose        );
    bool pass_jetId   = 0;

    float pileupJetId = ntP->jet_pileupJetId->at(idx);

    if( fabs(_eta) >= 0 && fabs(_eta) < 2.5 )
        if( pileupJetId > -0.63 ) pass_jetId = 1;
        else if( fabs(_eta) >= 2.5 && fabs(_eta) < 2.75 )
            if( pileupJetId > -0.60 ) pass_jetId = 1;
            else if( fabs(_eta) >= 2.75 && fabs(_eta) < 3.0 )
                if( pileupJetId > -0.55 ) pass_jetId = 1;
                else if( fabs(_eta) >= 3.0 && fabs(_eta) < 5.2 )
                    if( pileupJetId > -0.45 ) pass_jetId = 1;

    bool pass_muOverlap = 1;
    int nMuon = nt->NtMuon->size();
    for(int im=0;im<nMuon;im++)
    {
        float dr = GetDeltaR(_eta,_phi,nt->NtMuon->at(im).eta(),nt->NtMuon->at(im).phi());
        if( dr < 0.4 && nt->NtMuon->at(im).pt() > 10. ) pass_muOverlap = 0; //&& nt->NtMuon->at(im).isTight() ) pass_muOverlap = 0;
    }  

    bool pass_elOverlap = 1;
    int nElectron = nt->NtElectron->size();
    for(int ie=0;ie<nElectron;ie++)
    {
        float dr = GetDeltaR(_eta,_phi,nt->NtElectron->at(ie).eta(),nt->NtElectron->at(ie).phi());
        if( dr < 0.4 && nt->NtElectron->at(ie).pt() > 10. ) pass_elOverlap = 0; //&& nt->NtElectron->at(ie).isTight() ) pass_elOverlap = 0;
    }

    bool pass_tauOverlap = 1;
    int nTau = nt->NtTau->size();
    for(int it=0;it<nTau;it++)
    {
        float dr = GetDeltaR(_eta,_phi,nt->NtTau->at(it).eta(),nt->NtTau->at(it).phi());
        if( dr < 0.4 && nt->NtTau->at(it).pt() > 10. ) pass_tauOverlap = 0; //&& nt->NtTau->at(it).isTight() ) pass_tauOverlap = 0;
    }

    bool isSelectionJet = ( pass_pt         &&
            pass_eta        &&
            pass_isLoose    &&
            pass_muOverlap  &&
            pass_elOverlap  &&
            pass_tauOverlap );

    cout<<std::setiosflags(ios::fixed)<<setprecision(5);
    // synchronization printout
    /*if( isSelectionJet ) std::cout  << nt->NtEvent->at(0).id()                                      << std::setw(10)
      << _pt                                                          << std::setw(10)
      << _eta                                                         << std::setw(10)
      << _phi                                                         << std::setw(10)
    //<< _isLoose                                                     << std::setw(10)
    << _E                                                           << std::setw(10)
    << _CSVv2                                                       << std::setw(10)
    << _metpt                                                       << std::setw(10)
    << _metphi                                                      << std::endl;
    */
    return isSelectionJet;
}
