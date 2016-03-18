#include "TH2F.h"
#include "TFile.h"

bool SortingLeptonPt( Lepton l1, Lepton l2)
{
    if( l1.pt() > l2.pt() ) return true;
    else return false;
}
