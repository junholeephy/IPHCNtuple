#ifndef TRUTH_H
#define TRUTH_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Base.h"

class Truth : public Base
{
 public:
   Truth();
   virtual ~Truth();

   int mc_truth_n()         {return _mc_truth_n;};
   std::vector<int> mc_truth_id()         {return _mc_truth_id;};
   std::vector<int> mc_truth_label()         {return _mc_truth_label;};
   std::vector<TLorentzVector> mc_truth_p4()         {return _mc_truth_p4;};

   int gen_n() {return _gen_n;};
   std::vector<float> gen_pt() {return _gen_pt;};
   std::vector<float> gen_eta() {return _gen_eta;};
   std::vector<float> gen_phi() {return _gen_phi;};
   std::vector<float> gen_m() {return _gen_m;};
   std::vector<int> gen_id() {return _gen_id;};
   std::vector<int> gen_status() {return _gen_status;};
   std::vector<int> gen_mother_id() {return _gen_mother_id;};
   
   void read();
   void init();
	
 protected:

   int _mc_truth_n;
   std::vector<int> _mc_truth_id;
   std::vector<int> _mc_truth_label;
   std::vector<TLorentzVector> _mc_truth_p4;
   
   int _gen_n;
   std::vector<float> _gen_pt;
   std::vector<float> _gen_eta;
   std::vector<float> _gen_phi;
   std::vector<float> _gen_m;
   std::vector<int> _gen_id;
   std::vector<int> _gen_status;
   std::vector<int> _gen_mother_id;
   
   ClassDef(Truth,1)
};

#endif
