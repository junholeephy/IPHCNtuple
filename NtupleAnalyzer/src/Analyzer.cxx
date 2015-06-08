#include "../include/Hist.h"
#include "../include/TransferFunc.h"

char *fin;
char *tool;
char *evc;
TTree *tr;

std::vector<Electron>             *v_Electron;
std::vector<Muon>             *v_Muon;
std::vector<Event>             *v_Event;
std::vector<Jet>             *v_Jet;
std::vector<Truth>             *v_Truth;

int main(int argc, char *argv[])
{
   if( argc < 3 )
     {
	std::cout << "Usage: ./Analyzer [input file] [tool] [evc]" << std::endl;
	exit(1);
     }   
   
   fin = argv[1];
   tool = argv[2];
   evc = argv[3];
   
   TChain f("Nt");

   std::ifstream infile;
   infile.open(fin);
   
   std::string ifile = "";
   while( getline(infile, ifile) )
     {
	std::cout << ifile.c_str() << std::endl;
	f.Add(ifile.c_str());
     }	
   
   infile.close();

   std::vector<Electron>             *v_Electron             = new std::vector<Electron>();
   std::vector<Muon>                 *v_Muon                 = new std::vector<Muon>();
   std::vector<Event>                 *v_Event                 = new std::vector<Event>();
   std::vector<Jet>                 *v_Jet                 = new std::vector<Jet>();
   std::vector<Truth>                 *v_Truth                 = new std::vector<Truth>();
   
   f.SetBranchAddress("Electron", &v_Electron);
   f.SetBranchAddress("Muon", &v_Muon);
   f.SetBranchAddress("Jet", &v_Jet);
   f.SetBranchAddress("Event", &v_Event);
   f.SetBranchAddress("Truth", &v_Truth);

   int nent = f.GetEntries();
   std::cout << "Number of input events in the file: " << nent << std::endl;

   std::cout << "Initialisation completed" << std::endl;
   
   if( strcmp(tool,"plot") == 0 )
     {		
	Hist hist;
	
	hist.init();
	hist.setElectron(v_Electron);
	hist.setMuon(v_Muon);
	hist.setEvent(v_Event);
	hist.setJet(v_Jet);
	
	if( strcmp(evc,"1") == 0 )
	  {
	     std::cout << "Running in data challenge mode" << std::endl;
	  }	     
	
	for(int i=0;i<nent;i++)
	  {	
	     f.GetEntry(i);
	     
	     std::string fcur = f.GetCurrentFile()->GetName();
	     if( strcmp(evc,"1") == 0 )
	       {		 
		  bool passFINAL = hist.printout(true);
	       }
	     else
	       {	
		  bool passFINAL = hist.printout(false);
		  hist.fill();
	       }	     
	  }   
	
	hist.close();
     }   
   else if( strcmp(tool,"tran") == 0 )
     {		
	TransferFunc tran;
	
	tran.init();
	tran.setElectron(v_Electron);
	tran.setMuon(v_Muon);
	tran.setEvent(v_Event);
	tran.setJet(v_Jet);
	tran.setTruth(v_Truth);

	std::cout << "Transfer function mode" << std::endl;
	
	for(int i=0;i<nent;i++)
	  {	
	     f.GetEntry(i);
	     
	     std::string fcur = f.GetCurrentFile()->GetName();
	     
	     tran.run();
	  }   
	
	tran.close();
     }
   else
     {
	std::cout << "Select a proper tool" << std::endl;
	exit(1);
     }   

}
