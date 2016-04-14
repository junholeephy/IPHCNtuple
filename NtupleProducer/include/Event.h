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

        int id()                            {return _id;};
        int run()                           {return _run;};
        int lumi()                          {return _lumi;};

        float rho()                         {return _rho;};

        float metpt()                       {return _metpt;};
        float metphi()                      {return _metphi;};
        float metpx()                       {return _metpt*cos(_metphi);};
        float metpy()                       {return _metpt*sin(_metphi);};
        float metsumet()                    {return _metsumet;};

        double metcov00()		    {return _metcov00;};
        double metcov01()		    {return _metcov01;};
        double metcov10()		    {return _metcov10;};
        double metcov11()		    {return _metcov11;};

        float metNoHF_pt()                  {return _metNoHF_pt;};
        float metNoHF_phi()                 {return _metNoHF_phi;};
        float metNoHF_sumet()               {return _metNoHF_sumet;};

        int pv_n()                          {return _pv_n;};
        float pv_z()                        {return _pv_z;};
        float pv_zError()                   {return _pv_zError;};

        float weight_scale_muF0p5()         {return _weight_scale_muF0p5;};
        float weight_scale_muF2()           {return _weight_scale_muF2;};
        float weight_scale_muR0p5()         {return _weight_scale_muR0p5;};
        float weight_scale_muR2()           {return _weight_scale_muR2;};

        float mc_weight()                   {return _mc_weight;};
        float mc_ptHat()                    {return _mc_ptHat;};
        int mc_pu_trueNumInt()              {return _mc_pu_trueNumInt;};

        int ev_trigger_pass()               {return _trigger_pass;};
        int ev_trigger_pass_byname()        {return _trigger_pass_byname;};
        int ev_trigger_pass_byname_1()      {return _trigger_pass_byname_1;};
        int ev_trigger_pass_byname_1_noDz() {return _trigger_pass_byname_1_noDz;};

        int tth_channel()                   {return _tth_channel;};

        float disc_TT()                     {return _disc_TT;};

        void read(bool isdata);
        void init();

    protected:

        int   _id;
        int   _run;
        int   _lumi;

        float _rho;

        float _metpt;
        float _metphi;
        float _metsumet;
	double _metcov00;
        double _metcov01;
        double _metcov10;
        double _metcov11;

        float _metNoHF_pt;
        float _metNoHF_phi;
        float _metNoHF_sumet;

        int   _pv_n;
        float _pv_z;
        float _pv_zError;

        float _weight_scale_muF0p5;
        float _weight_scale_muF2;
        float _weight_scale_muR0p5;
        float _weight_scale_muR2;
	
        std::vector<float> _pdf_weights;

        float _mc_weight;
        float _mc_ptHat;
        int   _mc_pu_trueNumInt;

        int   _trigger_pass;
        int   _trigger_pass_byname;
        int   _trigger_pass_byname_1;
        int   _trigger_pass_byname_1_noDz;

        float _disc_TT;

    public:
        int _tth_channel;

        ClassDef(Event,1)
};

#endif
