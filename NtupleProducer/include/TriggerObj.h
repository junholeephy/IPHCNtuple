#ifndef TriggerObj_H
#define TriggerObj_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Base.h"

class TriggerObj : public Base
{
    public:

        TriggerObj();
        virtual ~TriggerObj();
                
        float pT()                               {return _pT;             };
        float eta()                              {return _eta;            };
        float phi()                              {return _phi;            };
        std::string collection()                 {return _collection;     };
        
	int filterIds_n()                        {return _filterIds_n;    };
	std::vector<int> filterIds()             {return _filterIds;      };
	int filterLabels_n()                     {return _filterLabels_n; };
	std::vector<std::string> filterLabels()  {return _filterLabels;   };
	
        int                      pathNamesAll_n()      {return  _pathNamesAll_n;     };
        std::vector<std::string> pathNamesAll()        {return _pathNamesAll;	     };
	std::vector<bool>        pathNamesAll_isL3()   {return _pathNamesAll_isL3;   };
	std::vector<bool>        pathNamesAll_isLF()   {return _pathNamesAll_isLF;   };
	std::vector<bool>        pathNamesAll_isBoth() {return _pathNamesAll_isBoth; };
	std::vector<bool>        pathNamesAll_isNone() {return _pathNamesAll_isNone; };

        int ID() {return _ID;};
        void read();
	void init();
	bool sel();
       
        
	// Other variables

    protected:
        
	int _ID;
	
        float _pT; 
	float _eta; 
	float _phi; 
	std::string _collection;
	
	int _filterIds_n; 
	std::vector<int> _filterIds;
	
	int _filterLabels_n;
	std::vector<std::string> _filterLabels;

	int _pathNamesAll_n;
        std::vector<std::string> _pathNamesAll;	
        std::vector<bool> _pathNamesAll_isL3;
	std::vector<bool> _pathNamesAll_isLF;
	std::vector<bool> _pathNamesAll_isBoth;
	std::vector<bool> _pathNamesAll_isNone;

        int _pathNamesAll_offset; // to decode flattrees
	int _filterLabels_offset; 
	int _filterIds_offset; 
	
        ClassDef(TriggerObj,1)
};

#endif
