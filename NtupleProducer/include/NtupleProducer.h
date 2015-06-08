#ifndef NTUPLEPRODUCER_H
#define NTUPLEPRODUCER_H

#include "Tree.h"
#include "Electron.h"
#include "Muon.h"
#include "Event.h"
#include "Jet.h"
#include "Truth.h"
#include "Ntuple.h"

extern Tree *ntP;
extern TChain *ch;
extern Ntuple *nt;
extern std::vector<int> *evdebug;

#endif
