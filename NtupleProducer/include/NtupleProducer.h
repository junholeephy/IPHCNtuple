#ifndef NTUPLEPRODUCER_H
#define NTUPLEPRODUCER_H

#include "Tree.h"
//#include "SKYNtuple/D3PDfast.h"
#include "Electron.h"
#include "Muon.h"
#include "Event.h"
#include "Jet.h"
#include "Truth.h"
#include "Ntuple.h"

//#include "TrigRootAnalysis/TrigDecisionToolD3PD.h"
//#include "TrigRootAnalysis/TrigConfigSvcD3PD.h"

extern Tree *ntP;
extern TChain *ch;
extern Ntuple *nt;
extern std::vector<int> *evdebug;
/*extern SKY::D3PDfast *ntPfast;
extern int isData;
extern int dbug;
extern int nocut;
extern D3PD::TrigDecisionToolD3PD *m_tdt;
*/
#endif
