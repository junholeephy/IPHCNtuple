#ifndef EVENT_H
#define EVENT_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Base.h"

class Event : public Base
{
 public:
   Event();
   virtual ~Event();
   
   int id()         {return _id;};
   int run()         {return _run;};
   int lumi()         {return _lumi;};

   float rho()         {return _rho;};
   
   float metpt()         {return _metpt;};
   float metphi()         {return _metphi;};
   float metsumet()         {return _metsumet;};
   
   float mc_weight()         {return _mc_weight;};
   
   int tth_channel()         {return _tth_channel;};
   
   void read();
   void init();
	
 protected:

   int _id;
   int _run;
   int _lumi;
   
   float _rho;

   float _metpt;
   float _metphi;
   float _metsumet;

   float _mc_weight;

 public:
   int _tth_channel;
   
   ClassDef(Event,1)
};

#endif
