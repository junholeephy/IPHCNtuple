{
   gROOT->SetBatch(1);
   
   TFile *f = TFile::Open("../../FlatTreeProducer/test/output.root");
   TTree *tr = (TTree*)f->Get("FlatTree/tree");
   TLorentzVector *mc_truth_tth_t1_p4 = 0;
   tr->SetBranchAddress("mc_truth_tth_t1_p4",&mc_truth_tth_t1_p4);
   for(int i=0;i<tr->GetEntries();i++)
     {
	tr->GetEntry(i);
	std::cout << mc_truth_tth_t1_p4->Pt() << std::endl;
     }   
   
   gApplication->Terminate();
}
