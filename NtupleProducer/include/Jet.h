#ifndef JET_H
#define JET_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Base.h"

class Jet : public Base
{
 public:
   Jet();
   virtual ~Jet();

   static bool sortPtPredicate(Jet lhs, Jet rhs)
     {return (lhs.pt() > rhs.pt());};
   
   int ID()    {return _ID;};
   
   void sel();
   
   // kinematics
   float E()         {return _E;};
   float pt()        {return _pt;};
   float eta()       {return _eta;};
   float phi()       {return _phi;};
   float m()         {return _m;};
   
   int ntrk()         {return _ntrk;};
   
   float CSV()         {return _CSV;};
   float CSVv2()         {return _CSVv2;};
   
   bool isTight()         {return _isTight;};
   
   void read();
   void init();
	
 protected:

   int _ID;
	
   float _E;
   float _pt;
   float _eta;
   float _phi;
   float _m;
   
   int _ntrk;
   
   float _CSV;
   float _CSVv2;

   bool _isTight;

   ClassDef(Jet,1)
};

#endif
