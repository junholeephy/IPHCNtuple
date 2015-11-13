#ifndef ELECTRON_H
#define ELECTRON_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Base.h"

class Electron : public Base
{
    public:

        Electron();
        virtual ~Electron();

        static bool sortPtPredicate(Electron lhs, Electron rhs) {return (lhs.pt() > rhs.pt());};

        int ID()                             {return _ID;};

        void setFakeType(int faketype)       {_fakeType = faketype;};
        int fakeType()                       {return _fakeType;};
        bool sel();

        // General informations
        float E()                            {return _E;};
        float pt()                           {return _pt;};
        float eta()                          {return _eta;};
        float phi()                          {return _phi;};
        float m()                            {return _m;};
        int   charge()                       {return _charge;};
        int   id()                           {return _id;};

        // Id
        bool isLoose()                       {return _isLoose;};
        bool isTight()                       {return _isTight;};
        bool isLooseMVA()                    {return _isLooseMVA;};
        bool isTightMVA()                    {return _isTightMVA;};

        // Variables for Id
        float dxy()                          {return _dxy;};
        float dz()                           {return _dz;};
        float iso()                          {return _iso;};
        bool passPtEta()                     {return _passPtEta;};

        // MVA
        float lepMVA()                       {return _lepMVA;};
        float lepMVA_neuRelIso()             {return _lepMVA_neuRelIso;};
        float lepMVA_chRelIso()              {return _lepMVA_chRelIso;};
        float lepMVA_jetDR()                 {return _lepMVA_jetDR;};
        float lepMVA_jetPtRatio()            {return _lepMVA_jetPtRatio;};
        float lepMVA_jetBTagCSV()            {return _lepMVA_jetBTagCSV;};
        float lepMVA_sip3d()                 {return _lepMVA_sip3d;};
        float lepMVA_dxy()                   {return _lepMVA_dxy;};
        float lepMVA_dz()                    {return _lepMVA_dz;};
        float lepMVA_mvaId()                 {return _lepMVA_mvaId;};

        // Other variables
        bool passChargeFlip()                  {return _passChargeFlip;};
        bool hasMatchedConversion()            {return _hasMatchedConversion;};

        bool pass_isGsfCtfScPixChargeConsistent() {return _pass_isGsfCtfScPixChargeConsistent;};

        float deltaEtaSuperClusterTrackAtVtx() {return _deltaEtaSuperClusterTrackAtVtx;};
        float deltaPhiSuperClusterTrackAtVtx() {return _deltaPhiSuperClusterTrackAtVtx;};
        float see()                            {return _see;};
        float hadronicOverEm()                 {return _hadronicOverEm;};
        float scleta()                         {return _scleta;};

        float dB3D()                           {return _dB3D;};
        float edB3D()                          {return _edB3D;};

        float effectiveArea(int dr,float eta);

        void read();
        void init();

    protected:

        int _ID;

        int _fakeType;

        // General informations
        float _E;
        float _pt;
        float _eta;
        float _phi;
        float _m;
        int   _charge;
        int   _id;

        // Id
        bool _isLooseCBId;
        bool _isMediumCBId;
        bool _isLoose;
        bool _isTight;
        bool _isLooseMVA;
        bool _isTightMVA;

        // Variables for Id
        float _dxy;
        float _dz;
        float _miniIso;
        int   _nlosthits;
        float _sip3d;
        bool  _isMedium;
        bool  _passCV;
        bool  _isPCC;
        float _iso;
        bool  _passPtEta;
        float _ip3d;
        float _ip3dErr;

        // MVA
        float _lepMVA;
        float _lepMVA_neuRelIso;
        float _lepMVA_chRelIso;
        float _lepMVA_jetDR;
        float _lepMVA_jetPtRatio;
        float _lepMVA_jetBTagCSV;
        float _lepMVA_sip3d;
        float _lepMVA_dxy;
        float _lepMVA_dz;
        float _lepMVA_mvaId;

        // Other variables
        bool _passChargeFlip;
        bool _hasMatchedConversion;

        bool _pass_isGsfCtfScPixChargeConsistent;

        float _deltaEtaSuperClusterTrackAtVtx;
        float _deltaPhiSuperClusterTrackAtVtx;
        float _see;
        float _hadronicOverEm;
        float _scleta;

        float _dB3D;
        float _edB3D;

        ClassDef(Electron,1)
};

#endif
