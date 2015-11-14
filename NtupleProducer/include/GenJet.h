#ifndef GENJET_H
#define GENJET_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Base.h"

class GenJet : public Base
{
    public:

        GenJet();
        virtual ~GenJet();
        
	int   ID()      {return _ID;};

        float pt()	{return _genJet_pt;};
        float eta()	{return _genJet_eta;};
        float phi()	{return _genJet_phi;};
        float E()	{return _genJet_E;};
        
        void read(); 
        void init();
        bool sel();

    protected:

        int    _ID;
	
        float  _genJet_pt;
        float  _genJet_eta;
        float  _genJet_phi;
        float  _genJet_E;

        ClassDef(GenJet,1)
};

#endif
