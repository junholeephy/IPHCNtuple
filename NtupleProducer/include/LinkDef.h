#include "Base.h"
#include "Electron.h"
#include "Muon.h"
#include "Event.h"
#include "Jet.h"
#include "Truth.h"

#ifdef __CINT__

#pragma link off all global;
#pragma link off all class;
#pragma link off all function;

#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedef;

//#pragma link C++ class std::pair<TString,double>+;
//#pragma link C++ class std::pair<unsigned int,double>+;
#pragma link C++ class Electron+;
#pragma link C++ class std::vector<Electron>+;
#pragma link C++ class Muon+;
#pragma link C++ class std::vector<Muon>+;
#pragma link C++ class Base+;
#pragma link C++ class std::vector<Base>+;
#pragma link C++ class Event+;
#pragma link C++ class std::vector<Event>+;
#pragma link C++ class Jet+;
#pragma link C++ class std::vector<Jet>+;
#pragma link C++ class Truth+;
#pragma link C++ class std::vector<Truth>+;

#endif
