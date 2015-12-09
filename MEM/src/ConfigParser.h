
#ifndef CONFIGPARSER_H
#define CONFIGPARSER_H

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

class ConfigParser
{
  public:
  ConfigParser();
  ~ConfigParser();

  void GetConfigFromFile(string);
  void LoadHypotheses(int*, string**, int**, int**); 
  void LoadIntegrationRange(double*, double*, double*);
  void LoadJetChoice(string*);

  ifstream fconf;

  int nHyp;
  string* sHyp;
  int* Hyp;

  int doTTLL, doTTHfl, doTTHsl, doTTW, doTTWJJ;
  int nPointsHypTTLL, nPointsHypTTHsl, nPointsHypTTHfl, nPointsHypTTW, nPointsHypTTWJJ;
  double valJetTFfracmin, valJetTFfracmax, valNeutMaxE;
  string valJetChoice;

  void ReadOptionValue(string*, int*);
  void ReadOptionValue(string*, double*);
  void ReadOptionValue(string*, string*);

  private:
};

ConfigParser::ConfigParser(){
}

void ConfigParser::GetConfigFromFile(string InputFile){

  fconf.open(InputFile.c_str());
  string line;
  string option;

  getline(fconf, line); 
  getline(fconf, line);
  getline(fconf, line);
  ReadOptionValue(&option, &doTTLL);
  ReadOptionValue(&option, &nPointsHypTTLL);
  ReadOptionValue(&option, &doTTHsl);
  ReadOptionValue(&option, &nPointsHypTTHsl);
  ReadOptionValue(&option, &doTTHfl);
  ReadOptionValue(&option, &nPointsHypTTHfl);
  ReadOptionValue(&option, &doTTW);
  ReadOptionValue(&option, &nPointsHypTTW);
  ReadOptionValue(&option, &doTTWJJ);
  ReadOptionValue(&option, &nPointsHypTTWJJ);

  getline(fconf, line);
  getline(fconf, line);
  getline(fconf, line);
  ReadOptionValue(&option, &valJetTFfracmin);
  ReadOptionValue(&option, &valJetTFfracmax);
  ReadOptionValue(&option, &valNeutMaxE);

  getline(fconf, line);
  getline(fconf, line);
  getline(fconf, line);
  ReadOptionValue(&option, &valJetChoice);

  fconf.close();
  return;
}

ConfigParser::~ConfigParser(){
}

void ConfigParser::ReadOptionValue(string* option, int* value){
  string svalue;
  fconf >> *option;
  fconf >> svalue;
  cout << *option << " "<<svalue << endl;
  *value = atoi(svalue.c_str());
  return;
}

void ConfigParser::ReadOptionValue(string* option, double* value){
  string svalue;
  fconf >> *option;
  fconf >> svalue;
  cout << *option << " "<<svalue << endl;
  *value = atof(svalue.c_str());
  return;
}

void ConfigParser::ReadOptionValue(string* option, string* value){
  string svalue;
  fconf >> *option;
  fconf >> svalue;
  cout << *option << " "<<svalue << endl;
  *value = svalue;
  return;
}

void ConfigParser::LoadHypotheses(int* nhyp, string** shyp, int** hyp, int** nPointsHyp){

  (*nhyp) = 0;

  if (doTTLL) (*nhyp)++;
  if (doTTHsl) (*nhyp)++;
  if (doTTHfl) (*nhyp)++;
  if (doTTW) (*nhyp)++;
  if (doTTWJJ) (*nhyp)++;

  (*shyp) = new string[(*nhyp)];
  (*hyp) = new int[(*nhyp)];
  (*nPointsHyp) = new int[(*nhyp)]; 

  int ih=0;
  if (doTTLL){
    (*shyp)[ih] = "TTLL";
    (*hyp)[ih] = kMEM_TTLL_TopAntitopDecay;
    (*nPointsHyp)[ih] = nPointsHypTTLL;
    ih++;
  }
  if (doTTHsl){
    (*shyp)[ih] = "TTHsl";
    (*hyp)[ih] = kMEM_TTH_TopAntitopHiggsSemiLepDecay;
    (*nPointsHyp)[ih] = nPointsHypTTHsl;
    ih++;
  }
  if (doTTHfl){
    (*shyp)[ih] = "TTHfl";
    (*hyp)[ih] = kMEM_TTH_TopAntitopHiggsDecay;
    (*nPointsHyp)[ih] = nPointsHypTTHfl;
    ih++;
  }
  if (doTTW){
    (*shyp)[ih] = "TTW";
    (*hyp)[ih] = kMEM_TTW_TopAntitopDecay;
    (*nPointsHyp)[ih] = nPointsHypTTW;
    ih++;
  }
  if (doTTWJJ){
    (*shyp)[ih] = "TTWJJ";
    (*hyp)[ih] = kMEM_TTWJJ_TopAntitopDecay;
    (*nPointsHyp)[ih] = nPointsHypTTWJJ;
    ih++;
  }

  for (int ih=0; ih<(*nhyp); ih++) cout << "Will run hyp "<<(*shyp)[ih]<<" code "<<(*hyp)[ih]<<" with "<<(*nPointsHyp)[ih]<<" iterations"<<endl;

  return;
}

void ConfigParser::LoadIntegrationRange(double* jetTFfracmin, double* jetTFfracmax, double* neutMaxE){

  *jetTFfracmin = valJetTFfracmin;
  *jetTFfracmax = valJetTFfracmax;
  *neutMaxE = valNeutMaxE;

  cout << "Integration will use JetTFfracmin="<< (*jetTFfracmin)<<" JetTFfracmax="<< (*jetTFfracmax) <<" NeutMaxE="<< (*neutMaxE) <<endl;

  return;
}

void ConfigParser::LoadJetChoice(string* jetChoice){

  *jetChoice = valJetChoice;
  cout << "In 2j categories, choose 2 jets with option "<< (*jetChoice) << endl;

  return;
}



#endif
