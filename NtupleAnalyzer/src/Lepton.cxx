#include "../include/Lepton.h"

Lepton::Lepton()
{
    _pt         =  0;
    _eta        =  0;
    _phi        =  0;
    _E          =  0;

    _p4.SetPtEtaPhiM(0,0,0,0);

    _id         =  0;

    _idx        = -1;
    _isElectron =  0; 
    _isMuon     =  0;
    
    _isFakeableTTH = false;
    _isTightTTH    = false;
    _lepMVA_TTH    =  0.;

    _charge     =  0;
}

Lepton::~Lepton()
{
}
