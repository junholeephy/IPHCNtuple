#include "include/NtupleProducer.h"

ClassImp(Truth)
    
Truth::Truth()
{
}

Truth::~Truth()
{
}

void Truth::read()
{
   int UNINT = -666;

   int gen_n = ntP->gen_n;
   int gen_n_sel = 0;
   for(int i=0;i<gen_n;i++)
     {
	int status = ntP->gen_status->at(i);
	if( status != 1 && status != 3 ) continue;
	_gen_pt.push_back(ntP->gen_pt->at(i));
	_gen_eta.push_back(ntP->gen_eta->at(i));
	_gen_phi.push_back(ntP->gen_phi->at(i));
	_gen_m.push_back(ntP->gen_m->at(i));
	_gen_id.push_back(ntP->gen_id->at(i));
	_gen_status.push_back(status);

	int mom = ntP->gen_mother_index->at(i);
	
	if( mom >= ntP->gen_index->size() )
	  {
	     
	     std::cout << "Problem with gen length" << std::endl;
	     exit(1);
	  }
	
	int mother_id = ntP->gen_id->at(mom);

	_gen_mother_id.push_back(mother_id);
	
	gen_n_sel++;
     }
/*
   _mc_truth_n = 0;
   
   // top as ttbar
   int mc_truth_t1_id = (ntP->mc_truth_t1_id) ? ntP->mc_truth_t1_id : UNINT;
   if( mc_truth_t1_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_t1_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_t1_p4));
	_mc_truth_label.push_back(2);
	_mc_truth_n++;
     }   
   int mc_truth_t2_id = (ntP->mc_truth_t2_id) ? ntP->mc_truth_t2_id : UNINT;
   if( mc_truth_t2_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_t2_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_t2_p4));
	_mc_truth_label.push_back(3);
	_mc_truth_n++;
     }   
   
   int mc_truth_tW1_id = (ntP->mc_truth_tW1_id) ? ntP->mc_truth_tW1_id : UNINT;
   if( mc_truth_tW1_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_tW1_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_tW1_p4));
	_mc_truth_label.push_back(20);
	_mc_truth_n++;
     }   
   int mc_truth_tb1_id = (ntP->mc_truth_tb1_id) ? ntP->mc_truth_tb1_id : UNINT;
   if( mc_truth_tb1_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_tb1_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_tb1_p4));
	_mc_truth_label.push_back(21);
	_mc_truth_n++;
     }   
   int mc_truth_tWl1_id = (ntP->mc_truth_tWl1_id) ? ntP->mc_truth_tWl1_id : UNINT;
   if( mc_truth_tWl1_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_tWl1_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_tWl1_p4));
	_mc_truth_label.push_back(220);
	_mc_truth_n++;
     }   
   int mc_truth_tWq11_id = (ntP->mc_truth_tWq11_id) ? ntP->mc_truth_tWq11_id : UNINT;
   if( mc_truth_tWq11_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_tWq11_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_tWq11_p4));
	_mc_truth_label.push_back(222);
	_mc_truth_n++;
     }   
   int mc_truth_tWq21_id = (ntP->mc_truth_tWq21_id) ? ntP->mc_truth_tWq21_id : UNINT;
   if( mc_truth_tWq21_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_tWq21_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_tWq21_p4));
	_mc_truth_label.push_back(223);
	_mc_truth_n++;
     }   

   int mc_truth_tW2_id = (ntP->mc_truth_tW2_id) ? ntP->mc_truth_tW2_id : UNINT;
   if( mc_truth_tW2_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_tW2_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_tW2_p4));
	_mc_truth_label.push_back(30);
	_mc_truth_n++;
     }   
   int mc_truth_tb2_id = (ntP->mc_truth_tb2_id) ? ntP->mc_truth_tb2_id : UNINT;
   if( mc_truth_tb2_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_tb2_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_tb2_p4));
	_mc_truth_label.push_back(31);
	_mc_truth_n++;
     }   
   int mc_truth_tWl2_id = (ntP->mc_truth_tWl2_id) ? ntP->mc_truth_tWl2_id : UNINT;
   if( mc_truth_tWl2_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_tWl2_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_tWl2_p4));
	_mc_truth_label.push_back(330);
	_mc_truth_n++;
     }   
   int mc_truth_tWq12_id = (ntP->mc_truth_tWq12_id) ? ntP->mc_truth_tWq12_id : UNINT;
   if( mc_truth_tWq12_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_tWq12_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_tWq12_p4));
	_mc_truth_label.push_back(332);
	_mc_truth_n++;
     }   
   int mc_truth_tWq22_id = (ntP->mc_truth_tWq22_id) ? ntP->mc_truth_tWq22_id : UNINT;
   if( mc_truth_tWq22_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_tWq22_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_tWq22_p4));
	_mc_truth_label.push_back(333);
	_mc_truth_n++;
     }

   // higgs
   int mc_truth_h0_id = ntP->mc_truth_h0_id;
   if( mc_truth_h0_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_h0_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_h0_p4));
	_mc_truth_label.push_back(1);
	_mc_truth_n++;
     }

   // WW
   int mc_truth_h0W1_id = ntP->mc_truth_h0W1_id;
   if( mc_truth_h0W1_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_h0W1_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_h0W1_p4));
	_mc_truth_label.push_back(12);
	_mc_truth_n++;
     }
   int mc_truth_h0Wl1_id = ntP->mc_truth_h0Wl1_id;
   if( mc_truth_h0Wl1_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_h0Wl1_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_h0Wl1_p4));
	_mc_truth_label.push_back(120);
	_mc_truth_n++;
     }
   int mc_truth_h0Wq11_id = ntP->mc_truth_h0Wq11_id;
   if( mc_truth_h0Wq11_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_h0Wq11_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_h0Wq11_p4));
	_mc_truth_label.push_back(122);
	_mc_truth_n++;
     }
   int mc_truth_h0Wq21_id = ntP->mc_truth_h0Wq21_id;
   if( mc_truth_h0Wq21_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_h0Wq21_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_h0Wq21_p4));
	_mc_truth_label.push_back(123);
	_mc_truth_n++;
     }
   int mc_truth_h0W2_id = ntP->mc_truth_h0W2_id;
   if( mc_truth_h0W2_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_h0W2_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_h0W2_p4));
	_mc_truth_label.push_back(13);
	_mc_truth_n++;
     }
   int mc_truth_h0Wl2_id = ntP->mc_truth_h0Wl2_id;
   if( mc_truth_h0Wl2_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_h0Wl2_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_h0Wl2_p4));
	_mc_truth_label.push_back(130);
	_mc_truth_n++;
     }
   int mc_truth_h0Wq12_id = ntP->mc_truth_h0Wq12_id;
   if( mc_truth_h0Wq12_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_h0Wq12_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_h0Wq12_p4));
	_mc_truth_label.push_back(132);
	_mc_truth_n++;
     }
   int mc_truth_h0Wq22_id = ntP->mc_truth_h0Wq22_id;
   if( mc_truth_h0Wq22_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_h0Wq22_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_h0Wq22_p4));
	_mc_truth_label.push_back(133);
	_mc_truth_n++;
     }

   // ZZ
   int mc_truth_h0Z1_id = ntP->mc_truth_h0Z1_id;
   if( mc_truth_h0Z1_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_h0Z1_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_h0Z1_p4));
	_mc_truth_label.push_back(14);
	_mc_truth_n++;
     }
   int mc_truth_h0Zl11_id = ntP->mc_truth_h0Zl11_id;
   if( mc_truth_h0Zl11_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_h0Zl11_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_h0Zl11_p4));
	_mc_truth_label.push_back(140);
	_mc_truth_n++;
     }
   int mc_truth_h0Zl21_id = ntP->mc_truth_h0Zl21_id;
   if( mc_truth_h0Zl21_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_h0Zl21_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_h0Zl21_p4));
	_mc_truth_label.push_back(141);
	_mc_truth_n++;
     }
   int mc_truth_h0Zq11_id = ntP->mc_truth_h0Zq11_id;
   if( mc_truth_h0Zq11_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_h0Zq11_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_h0Zq11_p4));
	_mc_truth_label.push_back(142);
	_mc_truth_n++;
     }
   int mc_truth_h0Zq21_id = ntP->mc_truth_h0Zq21_id;
   if( mc_truth_h0Zq21_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_h0Zq21_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_h0Zq21_p4));
	_mc_truth_label.push_back(143);
	_mc_truth_n++;
     }

   int mc_truth_h0Z2_id = ntP->mc_truth_h0Z2_id;
   if( mc_truth_h0Z2_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_h0Z2_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_h0Z2_p4));
	_mc_truth_label.push_back(15);
	_mc_truth_n++;
     }
   int mc_truth_h0Zl12_id = ntP->mc_truth_h0Zl12_id;
   if( mc_truth_h0Zl12_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_h0Zl12_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_h0Zl12_p4));
	_mc_truth_label.push_back(150);
	_mc_truth_n++;
     }
   int mc_truth_h0Zl22_id = ntP->mc_truth_h0Zl22_id;
   if( mc_truth_h0Zl22_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_h0Zl22_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_h0Zl22_p4));
	_mc_truth_label.push_back(151);
	_mc_truth_n++;
     }
   int mc_truth_h0Zq12_id = ntP->mc_truth_h0Zq12_id;
   if( mc_truth_h0Zq12_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_h0Zq12_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_h0Zq12_p4));
	_mc_truth_label.push_back(152);
	_mc_truth_n++;
     }
   int mc_truth_h0Zq22_id = ntP->mc_truth_h0Zq22_id;
   if( mc_truth_h0Zq22_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_h0Zq22_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_h0Zq22_p4));
	_mc_truth_label.push_back(153);
	_mc_truth_n++;
     }

   // tautau
   int mc_truth_h0tau1_id = ntP->mc_truth_h0tau1_id;
   if( mc_truth_h0tau1_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_h0tau1_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_h0tau1_p4));
	_mc_truth_label.push_back(16);
	_mc_truth_n++;
     }
   int mc_truth_h0tau2_id = ntP->mc_truth_h0tau2_id;
   if( mc_truth_h0tau2_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_h0tau2_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_h0tau2_p4));
	_mc_truth_label.push_back(17);
	_mc_truth_n++;
     }

   // bbbar
   int mc_truth_h0b1_id = ntP->mc_truth_h0b1_id;
   if( mc_truth_h0b1_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_h0b1_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_h0b1_p4));
	_mc_truth_label.push_back(18);
	_mc_truth_n++;
     }
   int mc_truth_h0b2_id = ntP->mc_truth_h0b2_id;
   if( mc_truth_h0b2_id != UNINT )
     {
	_mc_truth_id.push_back(mc_truth_h0b2_id);
	_mc_truth_p4.push_back(*(ntP->mc_truth_h0b2_p4));
	_mc_truth_label.push_back(19);
	_mc_truth_n++;
     }*/
   
}

void Truth::init()
{
   _mc_truth_n = 0;
   
   _mc_truth_id.clear();
   _mc_truth_p4.clear();
   _mc_truth_label.clear();

   _gen_pt.clear();
   _gen_eta.clear();
   _gen_phi.clear();
   _gen_m.clear();
   _gen_id.clear();
   _gen_status.clear();
   _gen_mother_id.clear();
}

