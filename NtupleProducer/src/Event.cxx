#include "include/NtupleProducer.h"

ClassImp(Event)
    
Event::Event()
{
}

Event::~Event()
{
}

void Event::read()
{
   _id   = ntP->ev_id;
   _run  = ntP->ev_run;
   _lumi = ntP->ev_lumi;
   _rho  = ntP->ev_rho;
   
   _metpt = ntP->met_pt;
   _metphi = ntP->met_phi;
   _metsumet = ntP->met_sumet;
   
   _mc_weight  = ntP->mc_weight;
}

void Event::init()
{
   _id        = -666;
   _run       = -666;
   _lumi      = -666;
   _rho       = -666;

   _metpt        = -666;
   _metphi       = -666;
   _metsumet      = -666;
   
   _mc_weight = -666;
   _tth_channel = -666;
}
